//! Profiles, profile-profile alignment and the shared progressive machinery.
//!
//! A *profile* is a block of already-aligned rows summarised column by column
//! as weighted residue counts plus a gap weight. Aligning two profiles is the
//! same Gotoh recurrence as [`crate::pairwise`], with the substitution score of
//! two columns replaced by the weighted average over all residue pairs, and
//! with **position-specific gap penalties** in the style of ClustalW (Thompson,
//! Higgins & Gibson 1994, *Nucleic Acids Res* 22:4673-4680, section "Gap
//! penalties"):
//!
//! * a column that already contains gaps is cheap to gap again;
//! * a column within eight residues of an existing gap is expensive;
//! * for protein, a column inside a run of five hydrophilic residues is cheap.
//!
//! Sequences are weighted so that a clade of near-identical sequences does not
//! outvote a lone divergent one; both ClustalW's tree-derived weights and the
//! position-based weights of Henikoff & Henikoff (1994, *J Mol Biol*
//! 243:574-578) are available.

use tolviewer_core::alphabet::is_gap;
use tolviewer_core::{Alphabet, Error, Result, GAP};

use crate::matrix::SubstMatrix;
use crate::pairwise::{OP_DELETE, OP_INSERT, OP_MATCH};
use crate::tree::GuideTree;
use crate::Progress;

/// Number of residue slots in a profile column (`A`..`Z`).
const NLET: usize = 26;
/// Index used for anything that is neither a gap nor `A`..`Z`.
const UNKNOWN: usize = 23; // 'X'

const NEG: f32 = -1e18;

/// ClustalW's default hydrophilic residue set.
const HYDROPHILIC: &[u8] = b"GPSNDQEKR";
/// Length of a hydrophilic run that triggers the reduced gap penalty.
const HYDRO_RUN: usize = 5;
/// Window in which an existing gap raises the penalty for opening a new one.
const GAP_PROXIMITY: usize = 8;

/// Map a residue byte to its profile slot; gaps map to `None`.
#[inline]
fn slot(c: u8) -> Option<usize> {
    if is_gap(c) {
        return None;
    }
    let u = c.to_ascii_uppercase();
    if u.is_ascii_uppercase() {
        Some((u - b'A') as usize)
    } else {
        Some(UNKNOWN)
    }
}

/// Everything the profile aligner needs that does not depend on the profiles.
#[derive(Debug, Clone)]
pub(crate) struct AlignCtx {
    pub mat: &'static SubstMatrix,
    /// Gap penalties, already rescaled to the matrix in use.
    pub gap_open: f32,
    pub gap_extend: f32,
    pub terminal_factor: f32,
    pub alphabet: Alphabet,
    /// Largest DP matrix, in cells, before the aligner falls back to a band.
    pub cell_budget: u64,
    /// Use the FFT segment finder to restrict the DP (MAFFT-style).
    pub use_fft: bool,
}

impl AlignCtx {
    /// Gap penalties expressed in "BLOSUM62 units" rescaled onto the matrix
    /// actually in use, so `gap_open = 10` means the same thing whether the
    /// user picked BLOSUM62 (mean diagonal 5.8) or the IUB matrix (1.9).
    pub(crate) fn new(
        mat: &'static SubstMatrix,
        gap_open: f32,
        gap_extend: f32,
        terminal_factor: f32,
        alphabet: Alphabet,
    ) -> Self {
        let scale = mat.mean_diagonal() / SubstMatrix::blosum62().mean_diagonal();
        AlignCtx {
            mat,
            gap_open: gap_open * scale,
            gap_extend: gap_extend * scale,
            terminal_factor: terminal_factor.max(0.0),
            alphabet,
            cell_budget: 64_000_000,
            use_fft: false,
        }
    }
}

/// A block of aligned rows plus the column statistics the DP needs.
#[derive(Debug, Clone)]
pub(crate) struct Profile {
    /// Gapped rows, all `width` long.
    pub rows: Vec<Vec<u8>>,
    /// Index of each row in the caller's original sequence list.
    pub ids: Vec<usize>,
    /// Weight of each row.
    pub weights: Vec<f32>,
    pub width: usize,
    /// Total weight of all rows.
    total: f32,
    /// Per column, the residues present as `(slot, weight)`.
    sparse: Vec<Vec<(u8, f32)>>,
    /// Position-specific multipliers for the gap-open and gap-extend penalties,
    /// one per column.
    gop_mult: Vec<f32>,
    gep_mult: Vec<f32>,
}

impl Profile {
    /// Build a profile from aligned rows. Rows shorter than the widest are
    /// padded with gaps, so a ragged block cannot produce a ragged profile.
    pub(crate) fn new(
        rows: Vec<Vec<u8>>,
        ids: Vec<usize>,
        weights: Vec<f32>,
        ctx: &AlignCtx,
    ) -> Profile {
        let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
        let mut rows = rows;
        for r in &mut rows {
            r.resize(width, GAP);
        }
        let total: f32 = weights.iter().sum::<f32>().max(f32::MIN_POSITIVE);

        let mut dense = vec![[0.0f32; NLET]; width];
        let mut gap_weight = vec![0.0f32; width];
        for (row, &w) in rows.iter().zip(weights.iter()) {
            for (c, &ch) in row.iter().enumerate() {
                match slot(ch) {
                    Some(s) => dense[c][s] += w,
                    None => gap_weight[c] += w,
                }
            }
        }
        let sparse: Vec<Vec<(u8, f32)>> = dense
            .iter()
            .map(|col| {
                col.iter()
                    .enumerate()
                    .filter(|(_, &v)| v > 0.0)
                    .map(|(i, &v)| (i as u8, v))
                    .collect()
            })
            .collect();

        let (gop_mult, gep_mult) = position_specific_penalties(&rows, &gap_weight, total, ctx);

        Profile { rows, ids, weights, width, total, sparse, gop_mult, gep_mult }
    }

    /// A profile holding a single ungapped sequence.
    pub(crate) fn single(seq: &[u8], id: usize, weight: f32, ctx: &AlignCtx) -> Profile {
        Profile::new(vec![seq.to_vec()], vec![id], vec![weight], ctx)
    }

    pub(crate) fn len(&self) -> usize {
        self.rows.len()
    }

    /// Consensus residue of each column, used by the FFT segment finder.
    pub(crate) fn consensus(&self) -> Vec<u8> {
        (0..self.width)
            .map(|c| {
                self.sparse[c]
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|&(s, _)| b'A' + s)
                    .unwrap_or(GAP)
            })
            .collect()
    }

    /// Weighted residue frequency of column `c` as a fraction of total weight,
    /// used by the FFT vectoriser.
    pub(crate) fn column_freq(&self, c: usize) -> impl Iterator<Item = (u8, f32)> + '_ {
        let inv = 1.0 / self.total;
        self.sparse[c].iter().map(move |&(s, w)| (s, w * inv))
    }
}

/// ClustalW-style position-specific gap penalty multipliers.
fn position_specific_penalties(
    rows: &[Vec<u8>],
    gap_weight: &[f32],
    total: f32,
    ctx: &AlignCtx,
) -> (Vec<f32>, Vec<f32>) {
    let width = gap_weight.len();
    let mut gop = vec![1.0f32; width];
    let mut gep = vec![1.0f32; width];
    if width == 0 {
        return (gop, gep);
    }

    // Distance from each column to the nearest column that contains a gap.
    let has_gap: Vec<bool> = gap_weight.iter().map(|&g| g > 0.0).collect();
    let mut dist = vec![usize::MAX; width];
    let mut d = usize::MAX;
    for c in 0..width {
        if has_gap[c] {
            d = 0;
        } else if d != usize::MAX {
            d = d.saturating_add(1);
        }
        dist[c] = d;
    }
    d = usize::MAX;
    for c in (0..width).rev() {
        if has_gap[c] {
            d = 0;
        } else if d != usize::MAX {
            d = d.saturating_add(1);
        }
        dist[c] = dist[c].min(d);
    }

    // Columns inside a run of at least HYDRO_RUN hydrophilic residues in any
    // single sequence (protein only).
    let mut hydro = vec![false; width];
    if !ctx.alphabet.is_nucleotide() {
        for row in rows {
            let mut run = 0usize;
            for (c, &ch) in row.iter().enumerate() {
                if HYDROPHILIC.contains(&ch.to_ascii_uppercase()) {
                    run += 1;
                } else {
                    run = 0;
                }
                if run >= HYDRO_RUN {
                    for h in hydro.iter_mut().take(c + 1).skip(c + 1 - run) {
                        *h = true;
                    }
                }
            }
        }
    }

    for c in 0..width {
        if has_gap[c] {
            // Thompson et al. 1994: "if there is a gap in the column, GOP is
            // reduced by 0.3 x (number of sequences without a gap)/(number of
            // sequences)" and the extension penalty is halved.
            let without = ((total - gap_weight[c]) / total).clamp(0.0, 1.0);
            gop[c] = 0.3 * without;
            gep[c] = 0.5;
        } else if hydro[c] {
            // "... reduced by one third within a stretch of five hydrophilic
            // residues".
            gop[c] = 1.0 / 3.0;
        } else if dist[c] != usize::MAX && dist[c] <= GAP_PROXIMITY {
            // "... increased close to existing gaps", by the published factor
            // 2 + (8 - d) * 2 / 8.
            gop[c] = 2.0 + (GAP_PROXIMITY - dist[c]) as f32 * 2.0 / GAP_PROXIMITY as f32;
        }
    }
    (gop, gep)
}

/// The band of columns the DP is allowed to visit, one inclusive range per row.
#[derive(Debug, Clone)]
pub(crate) struct Band {
    pub lo: Vec<u32>,
    pub hi: Vec<u32>,
}

impl Band {
    /// No restriction at all.
    pub(crate) fn full(n: usize, m: usize) -> Band {
        Band { lo: vec![0; n + 1], hi: vec![m as u32; n + 1] }
    }

    /// A band of half-width `w` around the main diagonal of an `n x m` matrix.
    pub(crate) fn diagonal(n: usize, m: usize, w: usize) -> Band {
        let mut lo = Vec::with_capacity(n + 1);
        let mut hi = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let centre = if n == 0 { 0.0 } else { i as f64 * m as f64 / n as f64 };
            lo.push((centre - w as f64).max(0.0) as u32);
            hi.push(((centre + w as f64).min(m as f64)) as u32);
        }
        // Guarantee the corners and monotonicity.
        lo[0] = 0;
        hi[n] = m as u32;
        for i in 1..=n {
            lo[i] = lo[i].max(lo[i - 1]);
            hi[i] = hi[i].max(hi[i - 1]);
        }
        for i in (0..n).rev() {
            hi[i] = hi[i].min(hi[i + 1]);
            lo[i] = lo[i].min(lo[i + 1]);
        }
        Band { lo, hi }
    }

    fn cells(&self) -> u64 {
        self.lo
            .iter()
            .zip(self.hi.iter())
            .map(|(&l, &h)| (h as u64).saturating_sub(l as u64) + 1)
            .sum()
    }

    /// Make the band legal: monotone, containing both corners, non-empty.
    pub(crate) fn sanitise(&mut self, n: usize, m: usize) {
        self.lo.resize(n + 1, 0);
        self.hi.resize(n + 1, m as u32);
        for i in 0..=n {
            self.hi[i] = self.hi[i].min(m as u32);
            self.lo[i] = self.lo[i].min(self.hi[i]);
        }
        self.lo[0] = 0;
        self.hi[n] = m as u32;
        for i in 1..=n {
            self.lo[i] = self.lo[i].max(self.lo[i - 1]);
            self.hi[i] = self.hi[i].max(self.hi[i - 1]);
        }
        for i in (0..n).rev() {
            self.hi[i] = self.hi[i].min(self.hi[i + 1]);
            self.lo[i] = self.lo[i].min(self.lo[i + 1]);
        }
        for i in 0..=n {
            if self.lo[i] > self.hi[i] {
                self.hi[i] = self.lo[i];
            }
        }
    }
}

/// Per-column substitution scores of profile `a` against each residue slot.
fn score_rows(a: &Profile, mat: &SubstMatrix) -> Vec<[f32; NLET]> {
    let inv = 1.0 / a.total;
    a.sparse
        .iter()
        .map(|col| {
            let mut out = [0.0f32; NLET];
            for &(x, w) in col {
                let f = w * inv;
                for (y, o) in out.iter_mut().enumerate() {
                    *o += f * mat.score(b'A' + x, b'A' + y as u8);
                }
            }
            out
        })
        .collect()
}

/// Gap penalty arrays for the DP, indexed by insertion position.
fn gap_arrays(p: &Profile, ctx: &AlignCtx) -> (Vec<f32>, Vec<f32>) {
    let w = p.width;
    let mut open = Vec::with_capacity(w + 1);
    let mut ext = Vec::with_capacity(w + 1);
    for i in 0..=w {
        // Position `i` sits before column `i`; use that column's multipliers,
        // falling back to the last column at the right-hand edge.
        let c = i.min(w.saturating_sub(1));
        let (gm, em) = if w == 0 { (1.0, 1.0) } else { (p.gop_mult[c], p.gep_mult[c]) };
        let terminal = i == 0 || i == w;
        let f = if terminal { ctx.terminal_factor } else { 1.0 };
        open.push(ctx.gap_open * gm * f);
        ext.push(ctx.gap_extend * em * f);
    }
    (open, ext)
}

/// Profile-profile alignment. Returns the traceback operations, where a
/// `DELETE` consumes a column of `a` and an `INSERT` a column of `b`.
pub(crate) fn align_profiles(a: &Profile, b: &Profile, ctx: &AlignCtx, band: &Band) -> Vec<u8> {
    let n = a.width;
    let m = b.width;
    if n == 0 || m == 0 {
        let mut ops = Vec::with_capacity(n + m);
        ops.extend(std::iter::repeat_n(OP_DELETE, n));
        ops.extend(std::iter::repeat_n(OP_INSERT, m));
        return ops;
    }

    let a_sc = score_rows(a, ctx.mat);
    let inv_b = 1.0 / b.total;
    let (vopen, vext) = gap_arrays(b, ctx); // gaps inserted into b, indexed by b column
    let (hopen, hext) = gap_arrays(a, ctx); // gaps inserted into a, indexed by a column

    // Flat traceback storage over the band.
    let mut off = Vec::with_capacity(n + 2);
    let mut acc: usize = 0;
    for i in 0..=n {
        off.push(acc);
        acc += (band.hi[i] - band.lo[i]) as usize + 1;
    }
    off.push(acc);
    let mut tb = vec![0u8; acc];

    let mut prev_m = vec![NEG; m + 1];
    let mut prev_x = vec![NEG; m + 1];
    let mut prev_y = vec![NEG; m + 1];
    let mut cur_m = vec![NEG; m + 1];
    let mut cur_x = vec![NEG; m + 1];
    let mut cur_y = vec![NEG; m + 1];

    // Row 0.
    {
        let (lo, hi) = (band.lo[0] as usize, band.hi[0] as usize);
        prev_m[0] = 0.0;
        for j in lo.max(1)..=hi {
            let opened = prev_m[j - 1] - hopen[0] - hext[0];
            let extended = prev_y[j - 1] - hext[0];
            if extended > opened {
                prev_y[j] = extended;
                tb[off[0] + j - lo] |= 2 << 4;
            } else {
                prev_y[j] = opened;
                // predecessor state 0 (match lane); the bits are already zero
            }
        }
    }

    for i in 1..=n {
        let (lo, hi) = (band.lo[i] as usize, band.hi[i] as usize);
        let (plo, phi) = (band.lo[i - 1] as usize, band.hi[i - 1] as usize);
        let base = off[i];
        let sc = &a_sc[i - 1];
        let hop = hopen[i];
        let hex = hext[i];
        // Reset the working row, including one cell to the left of the band so
        // stale values from older rows cannot leak in.
        for j in lo.saturating_sub(1)..=hi {
            cur_m[j] = NEG;
            cur_x[j] = NEG;
            cur_y[j] = NEG;
        }
        for j in lo..=hi {
            let idx = base + j - lo;
            // Match: from (i-1, j-1).
            if j >= 1 && j > plo && j - 1 <= phi {
                let (best, st) = argmax3(prev_m[j - 1], prev_x[j - 1], prev_y[j - 1]);
                let mut s = 0.0f32;
                for &(y, w) in &b.sparse[j - 1] {
                    s += sc[y as usize] * w * inv_b;
                }
                cur_m[j] = best + s;
                tb[idx] |= st;
            }
            // Gap in b: from (i-1, j).
            if j >= plo && j <= phi {
                let opened = prev_m[j].max(prev_y[j]) - vopen[j] - vext[j];
                let ost = if prev_m[j] >= prev_y[j] { 0 } else { 2 };
                let extended = prev_x[j] - vext[j];
                if extended > opened {
                    cur_x[j] = extended;
                    tb[idx] |= 1 << 2;
                } else {
                    cur_x[j] = opened;
                    tb[idx] |= ost << 2;
                }
            }
            // Gap in a: from (i, j-1).
            if j >= 1 && j > lo.saturating_sub(1) {
                let opened = cur_m[j - 1].max(cur_x[j - 1]) - hop - hex;
                let ost = if cur_m[j - 1] >= cur_x[j - 1] { 0 } else { 1 };
                let extended = cur_y[j - 1] - hex;
                if extended > opened {
                    cur_y[j] = extended;
                    tb[idx] |= 2 << 4;
                } else {
                    cur_y[j] = opened;
                    tb[idx] |= ost << 4;
                }
            }
        }
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_x, &mut cur_x);
        std::mem::swap(&mut prev_y, &mut cur_y);
    }

    let (_, mut state) = argmax3(prev_m[m], prev_x[m], prev_y[m]);
    let mut ops = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        // On the top or left edge only one move is possible, whatever lane the
        // traceback thinks it is in; taking it here keeps the indices in range.
        if i == 0 {
            ops.push(OP_INSERT);
            j -= 1;
            continue;
        }
        if j == 0 {
            ops.push(OP_DELETE);
            i -= 1;
            continue;
        }
        let lo = band.lo[i] as usize;
        if j < lo || j > band.hi[i] as usize {
            while i > 0 {
                ops.push(OP_DELETE);
                i -= 1;
            }
            while j > 0 {
                ops.push(OP_INSERT);
                j -= 1;
            }
            break;
        }
        let idx = off[i] + j - lo;
        match state {
            0 => {
                ops.push(OP_MATCH);
                state = tb[idx] & 0b11;
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(OP_DELETE);
                state = (tb[idx] >> 2) & 0b11;
                i -= 1;
            }
            _ => {
                ops.push(OP_INSERT);
                state = (tb[idx] >> 4) & 0b11;
                j -= 1;
            }
        }
    }
    ops.reverse();
    ops
}

#[inline]
fn argmax3(m: f32, x: f32, y: f32) -> (f32, u8) {
    let (mut v, mut s) = (m, 0u8);
    if x > v {
        v = x;
        s = 1;
    }
    if y > v {
        v = y;
        s = 2;
    }
    (v, s)
}

/// Align two profiles, choosing an unrestricted or a banded DP depending on
/// the size of the problem and on whether the FFT segment finder is enabled.
pub(crate) fn align_profiles_auto(a: &Profile, b: &Profile, ctx: &AlignCtx) -> Vec<u8> {
    let (n, m) = (a.width, b.width);
    if n == 0 || m == 0 {
        return align_profiles(a, b, ctx, &Band::full(n, m));
    }
    // Two lone sequences too big for a quadratic traceback: hand them to the
    // linear-space pairwise aligner, which is exact and needs O(n + m) memory.
    // This is the path a user aligning two organellar genomes takes.
    if a.len() == 1
        && b.len() == 1
        && (n as u64 + 1) * (m as u64 + 1) > crate::pairwise::LINEAR_SPACE_CELLS
        && !a.rows[0].iter().any(|&c| is_gap(c))
        && !b.rows[0].iter().any(|&c| is_gap(c))
    {
        return crate::pairwise::global_ops(
            &a.rows[0],
            &b.rows[0],
            ctx.mat,
            ctx.gap_open,
            ctx.gap_extend,
            ctx.terminal_factor,
        );
    }

    let mut band = if ctx.use_fft { crate::mafft::maybe_band(a, b, ctx) } else { None }
        .unwrap_or_else(|| Band::full(n, m));
    band.sanitise(n, m);

    if band.cells() > ctx.cell_budget {
        // Still too big: fall back to a band around the diagonal sized to the
        // budget. Approximate, and the only alternative to refusing outright.
        let w = (ctx.cell_budget / (2 * (n as u64 + 1))).max(64) as usize;
        band = Band::diagonal(n, m, w);
        band.sanitise(n, m);
    }
    align_profiles(a, b, ctx, &band)
}

/// Merge two profiles along a traceback.
pub(crate) fn merge(a: &Profile, b: &Profile, ops: &[u8], ctx: &AlignCtx) -> Profile {
    let width = ops.len();
    let mut rows: Vec<Vec<u8>> = Vec::with_capacity(a.len() + b.len());
    for row in &a.rows {
        let mut out = Vec::with_capacity(width);
        let mut i = 0usize;
        for &op in ops {
            match op {
                OP_INSERT => out.push(GAP),
                _ => {
                    out.push(row[i]);
                    i += 1;
                }
            }
        }
        rows.push(out);
    }
    for row in &b.rows {
        let mut out = Vec::with_capacity(width);
        let mut j = 0usize;
        for &op in ops {
            match op {
                OP_DELETE => out.push(GAP),
                _ => {
                    out.push(row[j]);
                    j += 1;
                }
            }
        }
        rows.push(out);
    }
    let ids: Vec<usize> = a.ids.iter().chain(b.ids.iter()).copied().collect();
    let weights: Vec<f32> = a.weights.iter().chain(b.weights.iter()).copied().collect();
    Profile::new(rows, ids, weights, ctx)
}

// ---------------------------------------------------------------------------
// Sequence weighting
// ---------------------------------------------------------------------------

/// Position-based sequence weights (Henikoff & Henikoff 1994).
///
/// In each column, a residue type seen `n` times among `k` distinct types
/// contributes `1 / (k * n)` to every sequence carrying it. Gaps count as a
/// type of their own, which keeps sequences that are mostly gaps from getting a
/// weight of zero. Weights are normalised to a mean of 1.
pub(crate) fn henikoff_weights(rows: &[Vec<u8>]) -> Vec<f32> {
    let n = rows.len();
    if n == 0 {
        return Vec::new();
    }
    let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut w = vec![0.0f32; n];
    for c in 0..width {
        let mut counts = [0u32; NLET + 1];
        for row in rows.iter() {
            let s = row.get(c).copied().map_or(NLET, |ch| slot(ch).unwrap_or(NLET));
            counts[s] += 1;
        }
        let distinct = counts.iter().filter(|&&x| x > 0).count();
        if distinct <= 1 {
            continue; // an invariant column carries no information about weights
        }
        for (i, row) in rows.iter().enumerate() {
            let s = row.get(c).copied().map_or(NLET, |ch| slot(ch).unwrap_or(NLET));
            w[i] += 1.0 / (distinct as f32 * counts[s] as f32);
        }
    }
    normalise(&mut w);
    w
}

/// ClustalW's tree-derived weights: each leaf's weight is the sum, over the
/// branches on its path to the root, of the branch length divided by the number
/// of leaves below that branch (Thompson et al. 1994, "Sequence weighting").
pub(crate) fn tree_weights(tree: &GuideTree, n: usize) -> Vec<f32> {
    let mut w = vec![0.0f32; n];
    fn walk(t: &GuideTree, acc: f32, w: &mut [f32]) {
        match t {
            GuideTree::Leaf(i) => {
                if let Some(slot) = w.get_mut(*i) {
                    *slot = acc;
                }
            }
            GuideTree::Node { left, right, left_len, right_len } => {
                let lc = left.leaf_count().max(1) as f32;
                let rc = right.leaf_count().max(1) as f32;
                walk(left, acc + left_len / lc, w);
                walk(right, acc + right_len / rc, w);
            }
        }
    }
    walk(tree, 0.0, &mut w);
    normalise(&mut w);
    w
}

fn normalise(w: &mut [f32]) {
    let n = w.len();
    if n == 0 {
        return;
    }
    let sum: f32 = w.iter().sum();
    if sum <= 0.0 || !sum.is_finite() {
        w.iter_mut().for_each(|x| *x = 1.0);
        return;
    }
    let k = n as f32 / sum;
    for x in w.iter_mut() {
        *x = (*x * k).max(1e-4);
    }
}

// ---------------------------------------------------------------------------
// Progressive alignment
// ---------------------------------------------------------------------------

/// Align `seqs` progressively down `tree`. Returns one gapped row per input
/// sequence, in the caller's original order.
pub(crate) fn progressive(
    seqs: &[Vec<u8>],
    tree: &GuideTree,
    weights: &[f32],
    ctx: &AlignCtx,
    progress: &dyn Progress,
    message: &str,
) -> Result<Vec<Vec<u8>>> {
    let profile = progressive_profile(seqs, tree, weights, ctx, progress, message)?;
    Ok(to_rows(&profile, seqs.len()))
}

/// Rows of a profile, put back in the caller's sequence order.
pub(crate) fn to_rows(p: &Profile, n: usize) -> Vec<Vec<u8>> {
    let mut out = vec![Vec::new(); n];
    for (row, &id) in p.rows.iter().zip(p.ids.iter()) {
        if let Some(slot) = out.get_mut(id) {
            *slot = row.clone();
        }
    }
    // Anything the tree never mentioned becomes an all-gap row so the result
    // stays rectangular.
    let width = p.width;
    for row in out.iter_mut() {
        if row.is_empty() && width > 0 {
            *row = vec![GAP; width];
        }
    }
    out
}

/// Progressive alignment, returning the root profile.
pub(crate) fn progressive_profile(
    seqs: &[Vec<u8>],
    tree: &GuideTree,
    weights: &[f32],
    ctx: &AlignCtx,
    progress: &dyn Progress,
    message: &str,
) -> Result<Profile> {
    let mut cache = SubtreeCache::default();
    progressive_cached(seqs, tree, weights, ctx, progress, message, &mut cache, &mut 0)
}

/// Alignments of already-computed subtrees, so that a second progressive pass
/// on a re-estimated tree can skip every subtree whose topology survived.
///
/// Only MUSCLE uses this; a single-pass engine passes a cache with `enabled`
/// false, because storing every intermediate profile costs O(n^2 * columns)
/// and would be pure waste.
#[derive(Debug, Default)]
pub(crate) struct SubtreeCache {
    entries: std::collections::HashMap<String, Profile>,
    /// Whether to store anything at all.
    pub enabled: bool,
    /// Cells stored so far, against [`SubtreeCache::BUDGET`].
    stored: u64,
}

impl SubtreeCache {
    /// Stop caching past this many residue cells (~256 MB of rows plus their
    /// column summaries). Beyond it MUSCLE simply recomputes the subtrees.
    const BUDGET: u64 = 256_000_000;

    /// A cache that will actually store subtree alignments.
    pub(crate) fn enabled() -> Self {
        SubtreeCache { enabled: true, ..Default::default() }
    }

    fn get(&self, key: &str) -> Option<&Profile> {
        if self.enabled {
            self.entries.get(key)
        } else {
            None
        }
    }

    fn insert(&mut self, key: String, p: &Profile) {
        if !self.enabled {
            return;
        }
        let cells = (p.len() as u64) * (p.width as u64);
        if self.stored + cells > Self::BUDGET {
            return;
        }
        self.stored += cells;
        self.entries.insert(key, p.clone());
    }
}

pub(crate) fn subtree_key(t: &GuideTree) -> String {
    let mut s = String::new();
    fn walk(t: &GuideTree, s: &mut String) {
        match t {
            GuideTree::Leaf(i) => s.push_str(&format!("{i},")),
            GuideTree::Node { left, right, .. } => {
                s.push('(');
                walk(left, s);
                walk(right, s);
                s.push(')');
            }
        }
    }
    walk(t, &mut s);
    s
}

/// Progressive alignment that reuses cached subtree alignments where the
/// topology is unchanged (MUSCLE stage 2, Edgar 2004, *Nucleic Acids Res*
/// 32:1792-1797, "Progressive alignment").
pub(crate) fn progressive_cached(
    seqs: &[Vec<u8>],
    tree: &GuideTree,
    weights: &[f32],
    ctx: &AlignCtx,
    progress: &dyn Progress,
    message: &str,
    cache: &mut SubtreeCache,
    reused: &mut usize,
) -> Result<Profile> {
    let total_nodes = tree.leaf_count().max(1);
    let mut done = 0usize;

    // Iterative post-order traversal; recursion would risk the stack on the
    // caterpillar trees that near-tied distances produce.
    enum Step<'a> {
        Visit(&'a GuideTree),
        Combine(&'a GuideTree),
    }
    let mut work = vec![Step::Visit(tree)];
    let mut stack: Vec<Profile> = Vec::new();

    while let Some(step) = work.pop() {
        match step {
            Step::Visit(node) => match node {
                GuideTree::Leaf(i) => {
                    let seq = seqs.get(*i).cloned().unwrap_or_default();
                    let w = weights.get(*i).copied().unwrap_or(1.0);
                    stack.push(Profile::single(&seq, *i, w, ctx));
                }
                GuideTree::Node { left, right, .. } => {
                    match cache.enabled.then(|| subtree_key(node)).and_then(|k| cache.get(&k)) {
                        Some(p) => {
                            *reused += 1;
                            stack.push(p.clone());
                        }
                        None => {
                            work.push(Step::Combine(node));
                            work.push(Step::Visit(right));
                            work.push(Step::Visit(left));
                        }
                    }
                }
            },
            Step::Combine(node) => {
                let b = stack.pop().ok_or_else(|| Error::algorithm("guide tree walk underflow"))?;
                let a = stack.pop().ok_or_else(|| Error::algorithm("guide tree walk underflow"))?;
                let ops = align_profiles_auto(&a, &b, ctx);
                let merged = merge(&a, &b, &ops, ctx);
                cache.insert(subtree_key(node), merged.clone());
                stack.push(merged);
                done += 1;
                if !progress.tick(done as f32 / total_nodes as f32, message) {
                    return Err(Error::Cancelled);
                }
            }
        }
    }
    stack.pop().ok_or_else(|| Error::algorithm("guide tree produced no alignment"))
}

// ---------------------------------------------------------------------------
// Sum-of-pairs score
// ---------------------------------------------------------------------------

/// Weighted sum-of-pairs score of an alignment.
///
/// This is the single objective function used by every accept/reject decision
/// in the crate, so refinement cannot be fooled by two different scorers
/// disagreeing. Substitutions are summed column by column from the weighted
/// residue counts; gaps are charged an affine penalty per sequence pair, with
/// a gap that lies outside a sequence's own first-to-last residue span treated
/// as terminal and scaled by `terminal_factor`.
///
/// Cost is O(width * (distinct residues^2 + 64)), i.e. independent of the
/// number of sequences beyond one pass to build the counts.
pub(crate) fn sp_score(rows: &[Vec<u8>], weights: &[f32], ctx: &AlignCtx) -> f32 {
    let k = rows.len();
    if k < 2 {
        return 0.0;
    }
    let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    if width == 0 {
        return 0.0;
    }
    let w: Vec<f32> = (0..k).map(|i| weights.get(i).copied().unwrap_or(1.0)).collect();

    // First and last residue column of each row, for terminal gap detection.
    let mut lead = vec![width; k];
    let mut tail = vec![0usize; k];
    for (i, row) in rows.iter().enumerate() {
        for (c, &ch) in row.iter().enumerate() {
            if !is_gap(ch) {
                if lead[i] == width {
                    lead[i] = c;
                }
                tail[i] = c;
            }
        }
    }

    // Penalty table over the eight (gap at c-1, gap at c, terminal) classes.
    let mut pen = [[0.0f32; 8]; 8];
    for (ci, row) in pen.iter_mut().enumerate() {
        for (cj, cell) in row.iter_mut().enumerate() {
            let (p0, n0, t0) = (ci & 1 != 0, ci & 2 != 0, ci & 4 != 0);
            let (p1, n1, t1) = (cj & 1 != 0, cj & 2 != 0, cj & 4 != 0);
            if n0 == n1 {
                continue; // both gapped or both aligned: no gap in this pair
            }
            let terminal = if n0 { t0 } else { t1 };
            let f = if terminal { ctx.terminal_factor } else { 1.0 };
            let continuing = p0 == n0 && p1 == n1;
            *cell = ctx.gap_extend * f + if continuing { 0.0 } else { ctx.gap_open * f };
        }
    }

    let mut score = 0.0f32;
    let mut counts = [0.0f32; NLET];
    let mut classes = [0.0f32; 8];
    let mut present: Vec<usize> = Vec::new();
    for c in 0..width {
        counts.fill(0.0);
        classes.fill(0.0);
        present.clear();
        // Sum of w_i^2 * S(x_i, x_i) over rows with a residue here, so the
        // self-pairs can be removed from the closed form below.
        let mut self_pairs = 0.0f32;
        for (i, row) in rows.iter().enumerate() {
            let ch = row.get(c).copied().unwrap_or(GAP);
            if let Some(s) = slot(ch) {
                if counts[s] == 0.0 {
                    present.push(s);
                }
                counts[s] += w[i];
                self_pairs += w[i] * w[i] * ctx.mat.score(b'A' + s as u8, b'A' + s as u8);
            }
            let prev_gap = if c == 0 {
                // Before the first column everything is "outside", so a leading
                // gap is a continuation rather than a fresh opening.
                is_gap(ch)
            } else {
                is_gap(row.get(c - 1).copied().unwrap_or(GAP))
            };
            let now_gap = is_gap(ch);
            let terminal = c < lead[i] || c > tail[i];
            let code = (prev_gap as usize) | ((now_gap as usize) << 1) | ((terminal as usize) << 2);
            classes[code] += w[i];
        }

        // sum_{i != j} w_i w_j S(x_i, x_j) = sum_x sum_y n_x n_y S - self pairs;
        // halve it for unordered pairs.
        let mut full = 0.0f32;
        for (a, &x) in present.iter().enumerate() {
            let cx = counts[x];
            full += cx * cx * ctx.mat.score(b'A' + x as u8, b'A' + x as u8);
            for &y in &present[a + 1..] {
                full += 2.0 * cx * counts[y] * ctx.mat.score(b'A' + x as u8, b'A' + y as u8);
            }
        }
        score += 0.5 * (full - self_pairs);

        for ci in 0..8 {
            if classes[ci] == 0.0 {
                continue;
            }
            for cj in (ci + 1)..8 {
                if classes[cj] == 0.0 {
                    continue;
                }
                score -= classes[ci] * classes[cj] * pen[ci][cj];
            }
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NoProgress;

    fn ctx() -> AlignCtx {
        AlignCtx::new(SubstMatrix::identity(), 4.0, 0.5, 1.0, Alphabet::Dna)
    }

    fn rows(v: &[&[u8]]) -> Vec<Vec<u8>> {
        v.iter().map(|r| r.to_vec()).collect()
    }

    #[test]
    fn identical_profiles_align_without_gaps() {
        let c = ctx();
        let a = Profile::new(rows(&[b"ACGTACGT"]), vec![0], vec![1.0], &c);
        let b = Profile::new(rows(&[b"ACGTACGT"]), vec![1], vec![1.0], &c);
        let ops = align_profiles_auto(&a, &b, &c);
        assert_eq!(ops.len(), 8);
        assert!(ops.iter().all(|&o| o == OP_MATCH));
    }

    #[test]
    fn merge_produces_a_rectangle() {
        let c = ctx();
        let a = Profile::new(rows(&[b"ACGTACGT", b"ACGTACGT"]), vec![0, 1], vec![1.0, 1.0], &c);
        let b = Profile::new(rows(&[b"ACGTTTACGT"]), vec![2], vec![1.0], &c);
        let ops = align_profiles_auto(&a, &b, &c);
        let m = merge(&a, &b, &ops, &c);
        assert_eq!(m.len(), 3);
        assert!(m.rows.iter().all(|r| r.len() == m.width));
        assert_eq!(m.width, ops.len());
    }

    #[test]
    fn existing_gaps_attract_new_gaps() {
        // Profile A already has a gap column; inserting into it should be
        // preferred over cutting a fresh gap elsewhere.
        let c = AlignCtx::new(SubstMatrix::identity(), 8.0, 0.5, 1.0, Alphabet::Dna);
        let a = Profile::new(
            rows(&[b"AAA-AAA", b"AAAGAAA", b"AAA-AAA", b"AAA-AAA"]),
            vec![0, 1, 2, 3],
            vec![1.0; 4],
            &c,
        );
        assert!(a.gop_mult[3] < 1.0, "gapped column must be cheap: {:?}", a.gop_mult);
        assert!(a.gop_mult[1] > 1.0, "column near a gap must be dear: {:?}", a.gop_mult);
    }

    #[test]
    fn hydrophilic_stretches_lower_the_penalty_for_protein() {
        let c = AlignCtx::new(SubstMatrix::blosum62(), 10.0, 0.5, 1.0, Alphabet::Protein);
        let p = Profile::new(rows(&[b"WWWWKRNDQEWWWW"]), vec![0], vec![1.0], &c);
        // Columns 4..10 are a run of hydrophilic residues.
        assert!(p.gop_mult[6] < 0.5, "{:?}", p.gop_mult);
        assert_eq!(p.gop_mult[0], 1.0);
    }

    #[test]
    fn henikoff_weights_downweight_duplicates() {
        let r = rows(&[b"AAAA", b"AAAA", b"AAAA", b"CCCC"]);
        let w = henikoff_weights(&r);
        assert_eq!(w.len(), 4);
        assert!(w[3] > w[0], "the odd one out must weigh more: {w:?}");
        assert!((w.iter().sum::<f32>() - 4.0).abs() < 1e-3);
    }

    #[test]
    fn tree_weights_follow_branch_lengths() {
        // ((0:0.1,1:0.1):0.4,2:0.5)
        let t = GuideTree::Node {
            left: Box::new(GuideTree::Node {
                left: Box::new(GuideTree::Leaf(0)),
                right: Box::new(GuideTree::Leaf(1)),
                left_len: 0.1,
                right_len: 0.1,
            }),
            right: Box::new(GuideTree::Leaf(2)),
            left_len: 0.4,
            right_len: 0.5,
        };
        let w = tree_weights(&t, 3);
        assert!(w[2] > w[0], "the lone deep leaf must weigh more: {w:?}");
        assert!((w[0] - w[1]).abs() < 1e-6);
    }

    #[test]
    fn sp_score_prefers_the_right_alignment() {
        let c = ctx();
        let good = rows(&[b"ACGTACGT", b"ACGTACGT"]);
        let bad = rows(&[b"ACGTACGT", b"TGCATGCA"]);
        let w = vec![1.0, 1.0];
        assert!(sp_score(&good, &w, &c) > sp_score(&bad, &w, &c));
    }

    #[test]
    fn sp_score_charges_gaps_affinely() {
        let c = ctx();
        let w = vec![1.0, 1.0];
        // One gap of two versus two gaps of one: the former must score higher.
        let one = rows(&[b"AA--AAAA", b"AAAAAAAA"]);
        let two = rows(&[b"AA-A-AAA", b"AAAAAAAA"]);
        assert!(sp_score(&one, &w, &c) > sp_score(&two, &w, &c));
    }

    #[test]
    fn progressive_alignment_of_three_sequences() {
        let c = ctx();
        let seqs: Vec<Vec<u8>> =
            vec![b"ACGTACGT".to_vec(), b"ACGTTACGT".to_vec(), b"ACGTACGT".to_vec()];
        let t = GuideTree::Node {
            left: Box::new(GuideTree::Node {
                left: Box::new(GuideTree::Leaf(0)),
                right: Box::new(GuideTree::Leaf(2)),
                left_len: 0.1,
                right_len: 0.1,
            }),
            right: Box::new(GuideTree::Leaf(1)),
            left_len: 0.2,
            right_len: 0.2,
        };
        let w = vec![1.0; 3];
        let out = progressive(&seqs, &t, &w, &c, &NoProgress, "aligning").unwrap();
        assert_eq!(out.len(), 3);
        let width = out[0].len();
        assert!(out.iter().all(|r| r.len() == width));
        assert_eq!(width, 9);
        for (i, r) in out.iter().enumerate() {
            let ungapped: Vec<u8> = r.iter().copied().filter(|&c| c != GAP).collect();
            assert_eq!(ungapped, seqs[i]);
        }
    }

    #[test]
    fn banded_and_full_agree_on_easy_input() {
        let c = ctx();
        let a = Profile::new(rows(&[b"ACGTACGTACGTACGT"]), vec![0], vec![1.0], &c);
        let b = Profile::new(rows(&[b"ACGTACGTACGTACGT"]), vec![1], vec![1.0], &c);
        let full = align_profiles(&a, &b, &c, &Band::full(a.width, b.width));
        let mut band = Band::diagonal(a.width, b.width, 3);
        band.sanitise(a.width, b.width);
        let banded = align_profiles(&a, &b, &c, &band);
        assert_eq!(full, banded);
    }

    #[test]
    fn empty_profiles_do_not_panic() {
        let c = ctx();
        let a = Profile::new(rows(&[b""]), vec![0], vec![1.0], &c);
        let b = Profile::new(rows(&[b"ACGT"]), vec![1], vec![1.0], &c);
        let ops = align_profiles_auto(&a, &b, &c);
        let m = merge(&a, &b, &ops, &c);
        assert_eq!(m.width, 4);
        assert_eq!(m.rows[0], b"----".to_vec());
    }
}
