//! Pairwise alignment with affine gaps.
//!
//! Global alignment follows Gotoh (1982, *J Mol Biol* 162:705-708): three
//! states per cell (aligned pair, gap in the second sequence, gap in the first)
//! so that a gap of length `L` costs `open + L * extend` instead of `L` times a
//! flat penalty.
//!
//! For long sequences the quadratic-memory traceback is replaced by the
//! linear-space divide and conquer of Myers & Miller (1988, *CABIOS* 4:11-17),
//! which is Hirschberg's algorithm extended to affine gaps. The switch happens
//! automatically at [`LINEAR_SPACE_CELLS`]; both paths return the same score
//! (checked by unit tests in this file), so callers never have to care which
//! one ran.
//!
//! ## Terminal gaps
//!
//! A gap run that touches the start or the end of the matrix is a *terminal*
//! gap. Its penalties are multiplied by `terminal_gap_factor`, so 0.0 gives
//! free end gaps (semi-global / "overlap" behaviour) and 1.0 gives a strict
//! global alignment. Terminal-ness depends only on the *global* row/column
//! index of the run, which is why the divide-and-conquer recursion slices the
//! penalty arrays rather than passing scalars: that keeps the linear-space path
//! exact for terminal gaps too, not just for the interior.

use std::ops::Range;

use crate::matrix::SubstMatrix;

/// Sentinel for "unreachable state". Small enough that adding a few thousand
/// penalties to it cannot reach `f32::MIN`, large enough that no real score
/// comes near it.
const NEG: f32 = -1e18;

/// Above this many DP cells, [`global`] switches to the linear-space
/// (Myers-Miller) path. 4M cells is ~4 MB of traceback bytes.
pub const LINEAR_SPACE_CELLS: u64 = 4_000_000;

/// Cells below which the recursion stops dividing and runs the quadratic DP.
const BASE_CELLS: u64 = 1 << 16;

/// Traceback operations, in alignment order.
pub(crate) const OP_MATCH: u8 = 0;
/// Consume one residue of `a` against a gap in `b`.
pub(crate) const OP_DELETE: u8 = 1;
/// Consume one residue of `b` against a gap in `a`.
pub(crate) const OP_INSERT: u8 = 2;

/// State indices used inside the DP and in the packed traceback.
const S_M: u8 = 0;
const S_X: u8 = 1;
const S_Y: u8 = 2;
const S_X0: u8 = 3;
/// Traceback marker for "a local alignment starts here".
const S_START: u8 = 3;

/// Position-dependent affine gap penalties.
///
/// `vopen[j]` / `vext[j]` apply to a gap in `b` (a vertical move) placed at
/// column `j`; `hopen[i]` / `hext[i]` apply to a gap in `a` (a horizontal move)
/// placed at row `i`. Making the penalties arrays rather than scalars is what
/// lets terminal gaps be charged differently, and it is the same shape the
/// profile aligner needs for ClustalW-style position-specific penalties.
#[derive(Debug, Clone)]
pub(crate) struct GapProfile {
    pub vopen: Vec<f32>,
    pub vext: Vec<f32>,
    pub hopen: Vec<f32>,
    pub hext: Vec<f32>,
}

impl GapProfile {
    /// Uniform penalties with terminal runs scaled by `terminal_factor`.
    pub(crate) fn uniform(n: usize, m: usize, open: f32, ext: f32, terminal_factor: f32) -> Self {
        let t = terminal_factor.max(0.0);
        let mut vopen = vec![open; m + 1];
        let mut vext = vec![ext; m + 1];
        let mut hopen = vec![open; n + 1];
        let mut hext = vec![ext; n + 1];
        // A gap in `b` sitting before the first or after the last residue of
        // `b` is a terminal gap; likewise for `a`.
        for idx in [0, m] {
            vopen[idx] = open * t;
            vext[idx] = ext * t;
        }
        for idx in [0, n] {
            hopen[idx] = open * t;
            hext[idx] = ext * t;
        }
        GapProfile { vopen, vext, hopen, hext }
    }
}

/// A rectangular subproblem, with its residues and penalties already sliced
/// (and reversed, for a backward pass).
struct Sub {
    a: Vec<u8>,
    b: Vec<u8>,
    vopen: Vec<f32>,
    vext: Vec<f32>,
    hopen: Vec<f32>,
    hext: Vec<f32>,
}

/// The four DP lanes for one row.
struct Lanes {
    m: Vec<f32>,
    x: Vec<f32>,
    x0: Vec<f32>,
    y: Vec<f32>,
}

impl Lanes {
    fn new(width: usize) -> Self {
        Lanes {
            m: vec![NEG; width],
            x: vec![NEG; width],
            x0: vec![NEG; width],
            y: vec![NEG; width],
        }
    }
    #[inline]
    fn best(&self, j: usize) -> f32 {
        self.m[j].max(self.x[j]).max(self.x0[j]).max(self.y[j])
    }
}

#[inline]
fn max2(a: f32, b: f32) -> f32 {
    if a >= b {
        a
    } else {
        b
    }
}

#[inline]
fn argmax3(a: f32, b: f32, c: f32, sa: u8, sb: u8, sc: u8) -> (f32, u8) {
    let (mut v, mut s) = (a, sa);
    if b > v {
        v = b;
        s = sb;
    }
    if c > v {
        v = c;
        s = sc;
    }
    (v, s)
}

#[inline]
fn argmax4(m: f32, x: f32, y: f32, x0: f32) -> (f32, u8) {
    let (mut v, mut s) = (m, S_M);
    if x > v {
        v = x;
        s = S_X;
    }
    if y > v {
        v = y;
        s = S_Y;
    }
    if x0 > v {
        v = x0;
        s = S_X0;
    }
    (v, s)
}

/// Score-only forward pass; returns the four lanes at the last row.
///
/// `enter_x` means the alignment is already inside a run of gaps in `b` when it
/// enters this subproblem, and that run's opening penalty has been charged by
/// the caller. That is the `x0` lane, which lives only in column 0 because a
/// subproblem is always entered at its top-left corner.
fn last_row(sub: &Sub, mat: &SubstMatrix, enter_x: bool) -> Lanes {
    let n = sub.a.len();
    let m = sub.b.len();
    let mut prev = Lanes::new(m + 1);
    let mut cur = Lanes::new(m + 1);

    if enter_x {
        prev.x0[0] = 0.0;
    } else {
        prev.m[0] = 0.0;
    }
    for j in 1..=m {
        let from = max2(max2(prev.m[j - 1], prev.x[j - 1]), prev.x0[j - 1]);
        prev.y[j] = max2(from - sub.hopen[0] - sub.hext[0], prev.y[j - 1] - sub.hext[0]);
    }

    for i in 1..=n {
        let ai = sub.a[i - 1];
        let hopen = sub.hopen[i];
        let hext = sub.hext[i];
        cur.m[0] = NEG;
        cur.y[0] = NEG;
        cur.x[0] =
            max2(max2(prev.m[0], prev.y[0]) - sub.vopen[0] - sub.vext[0], prev.x[0] - sub.vext[0]);
        cur.x0[0] = if enter_x { prev.x0[0] - sub.vext[0] } else { NEG };
        for j in 1..=m {
            let vopen = sub.vopen[j];
            let vext = sub.vext[j];
            let diag =
                max2(max2(prev.m[j - 1], prev.x[j - 1]), max2(prev.x0[j - 1], prev.y[j - 1]));
            cur.m[j] = diag + mat.score(ai, sub.b[j - 1]);
            cur.x[j] = max2(max2(prev.m[j], prev.y[j]) - vopen - vext, prev.x[j] - vext);
            cur.x0[j] = NEG;
            let left = max2(max2(cur.m[j - 1], cur.x[j - 1]), cur.x0[j - 1]);
            cur.y[j] = max2(left - hopen - hext, cur.y[j - 1] - hext);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev
}

/// Quadratic-memory DP with traceback. Returns the score and the operations.
fn full_dp(sub: &Sub, mat: &SubstMatrix, enter_x: bool, exit_x: bool) -> (f32, Vec<u8>) {
    let n = sub.a.len();
    let m = sub.b.len();
    let stride = m + 1;
    let mut tb = vec![0u8; stride * (n + 1)];
    let mut prev = Lanes::new(stride);
    let mut cur = Lanes::new(stride);

    if enter_x {
        prev.x0[0] = 0.0;
    } else {
        prev.m[0] = 0.0;
    }
    // The loop walks four lanes and the traceback array together, so an
    // iterator over one of them would not help.
    #[allow(clippy::needless_range_loop)]
    for j in 1..=m {
        let (from, st) = argmax3(prev.m[j - 1], prev.x[j - 1], prev.x0[j - 1], S_M, S_X, S_X0);
        let opened = from - sub.hopen[0] - sub.hext[0];
        let extended = prev.y[j - 1] - sub.hext[0];
        if extended > opened {
            prev.y[j] = extended;
            tb[j] |= S_Y << 4;
        } else {
            prev.y[j] = opened;
            tb[j] |= st << 4;
        }
    }

    for i in 1..=n {
        let ai = sub.a[i - 1];
        let hopen = sub.hopen[i];
        let hext = sub.hext[i];
        let row = i * stride;
        cur.m[0] = NEG;
        cur.y[0] = NEG;
        {
            let (from, st) =
                if prev.m[0] >= prev.y[0] { (prev.m[0], S_M) } else { (prev.y[0], S_Y) };
            let opened = from - sub.vopen[0] - sub.vext[0];
            let extended = prev.x[0] - sub.vext[0];
            if extended > opened {
                cur.x[0] = extended;
                tb[row] |= S_X << 2;
            } else {
                cur.x[0] = opened;
                tb[row] |= st << 2;
            }
            cur.x0[0] = if enter_x { prev.x0[0] - sub.vext[0] } else { NEG };
        }
        for j in 1..=m {
            let idx = row + j;
            let vopen = sub.vopen[j];
            let vext = sub.vext[j];
            let (diag, dst) = argmax4(prev.m[j - 1], prev.x[j - 1], prev.y[j - 1], prev.x0[j - 1]);
            cur.m[j] = diag + mat.score(ai, sub.b[j - 1]);
            tb[idx] |= dst;

            let (from, st) =
                if prev.m[j] >= prev.y[j] { (prev.m[j], S_M) } else { (prev.y[j], S_Y) };
            let opened = from - vopen - vext;
            let extended = prev.x[j] - vext;
            if extended > opened {
                cur.x[j] = extended;
                tb[idx] |= S_X << 2;
            } else {
                cur.x[j] = opened;
                tb[idx] |= st << 2;
            }
            cur.x0[j] = NEG;

            let (lfrom, lst) = argmax3(cur.m[j - 1], cur.x[j - 1], cur.x0[j - 1], S_M, S_X, S_X0);
            let lopened = lfrom - hopen - hext;
            let lextended = cur.y[j - 1] - hext;
            if lextended > lopened {
                cur.y[j] = lextended;
                tb[idx] |= S_Y << 4;
            } else {
                cur.y[j] = lopened;
                tb[idx] |= lst << 4;
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }

    let (score, mut state) = if exit_x {
        if prev.x[m] >= prev.x0[m] {
            (prev.x[m], S_X)
        } else {
            (prev.x0[m], S_X0)
        }
    } else {
        argmax4(prev.m[m], prev.x[m], prev.y[m], prev.x0[m])
    };

    let mut ops = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let idx = i * stride + j;
        match state {
            S_M => {
                ops.push(OP_MATCH);
                state = tb[idx] & 0b11;
                i -= 1;
                j -= 1;
            }
            S_X => {
                ops.push(OP_DELETE);
                state = (tb[idx] >> 2) & 0b11;
                i -= 1;
            }
            S_X0 => {
                ops.push(OP_DELETE);
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
    (score, ops)
}

/// Everything the recursion needs about the *global* problem.
struct Problem<'a> {
    a: &'a [u8],
    b: &'a [u8],
    mat: &'a SubstMatrix,
    gp: GapProfile,
}

impl Problem<'_> {
    fn forward(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Sub {
        Sub {
            a: self.a[r0..r1].to_vec(),
            b: self.b[c0..c1].to_vec(),
            vopen: self.gp.vopen[c0..=c1].to_vec(),
            vext: self.gp.vext[c0..=c1].to_vec(),
            hopen: self.gp.hopen[r0..=r1].to_vec(),
            hext: self.gp.hext[r0..=r1].to_vec(),
        }
    }

    fn backward(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Sub {
        Sub {
            a: self.a[r0..r1].iter().rev().copied().collect(),
            b: self.b[c0..c1].iter().rev().copied().collect(),
            vopen: self.gp.vopen[c0..=c1].iter().rev().copied().collect(),
            vext: self.gp.vext[c0..=c1].iter().rev().copied().collect(),
            hopen: self.gp.hopen[r0..=r1].iter().rev().copied().collect(),
            hext: self.gp.hext[r0..=r1].iter().rev().copied().collect(),
        }
    }
}

/// Myers-Miller divide and conquer. Appends the operations for
/// `a[r0..r1] x b[c0..c1]` to `out`.
fn linear_path(
    p: &Problem<'_>,
    r0: usize,
    r1: usize,
    c0: usize,
    c1: usize,
    enter_x: bool,
    exit_x: bool,
    out: &mut Vec<u8>,
) {
    let n = r1 - r0;
    let m = c1 - c0;
    if n <= 1 || (n as u64 + 1) * (m as u64 + 1) <= BASE_CELLS {
        let sub = p.forward(r0, r1, c0, c1);
        let (_, ops) = full_dp(&sub, p.mat, enter_x, exit_x);
        out.extend_from_slice(&ops);
        return;
    }

    let mid = n / 2;
    let (best_j, best_type2) = {
        let f = last_row(&p.forward(r0, r0 + mid, c0, c1), p.mat, enter_x);
        let r = last_row(&p.backward(r0 + mid, r1, c0, c1), p.mat, exit_x);
        let mut best = NEG;
        let mut best_j = 0usize;
        let mut type2 = false;
        for j in 0..=m {
            let rj = m - j;
            let t1 = f.best(j) + r.best(rj);
            // A run of gaps in `b` crossing the split row is charged an opening
            // penalty by both halves; refund one. The `x0` lanes carry no
            // opening penalty of their own, so nothing is refunded for those.
            let fx = max2(f.x[j], f.x0[j]);
            let refund = p.gp.vopen[c0 + j];
            let t2 = max2(fx + r.x[rj] + refund, fx + r.x0[rj]);
            if t1 > best {
                best = t1;
                best_j = j;
                type2 = false;
            }
            if t2 > best {
                best = t2;
                best_j = j;
                type2 = true;
            }
        }
        (best_j, type2)
    };

    let j = best_j;
    linear_path(p, r0, r0 + mid, c0, c0 + j, enter_x, best_type2, out);
    linear_path(p, r0 + mid, r1, c0 + j, c1, best_type2, exit_x, out);
}

/// Turn a traceback into the two gapped rows.
pub(crate) fn rows_from_ops(a: &[u8], b: &[u8], ops: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut ra = Vec::with_capacity(ops.len());
    let mut rb = Vec::with_capacity(ops.len());
    let (mut i, mut j) = (0usize, 0usize);
    for &op in ops {
        match op {
            OP_MATCH => {
                ra.push(a[i]);
                rb.push(b[j]);
                i += 1;
                j += 1;
            }
            OP_DELETE => {
                ra.push(a[i]);
                rb.push(tolviewer_core::GAP);
                i += 1;
            }
            _ => {
                ra.push(tolviewer_core::GAP);
                rb.push(b[j]);
                j += 1;
            }
        }
    }
    (ra, rb)
}

/// Score an alignment given as operations, using the same gap model as the DP.
fn score_ops(a: &[u8], b: &[u8], mat: &SubstMatrix, gp: &GapProfile, ops: &[u8]) -> f32 {
    let mut score = 0.0f32;
    let (mut i, mut j) = (0usize, 0usize);
    let mut prev = u8::MAX;
    for &op in ops {
        match op {
            OP_MATCH => {
                score += mat.score(a[i], b[j]);
                i += 1;
                j += 1;
            }
            OP_DELETE => {
                if prev != OP_DELETE {
                    score -= gp.vopen[j];
                }
                score -= gp.vext[j];
                i += 1;
            }
            _ => {
                if prev != OP_INSERT {
                    score -= gp.hopen[i];
                }
                score -= gp.hext[i];
                j += 1;
            }
        }
        prev = op;
    }
    score
}

/// Global alignment with affine gaps (Gotoh). Returns the two gapped rows and
/// the score. Terminal gaps are charged in full; use [`global_ends`] for free
/// or discounted end gaps.
///
/// `gap_open` and `gap_extend` are positive numbers; a gap of length `L` costs
/// `gap_open + L * gap_extend`. Sequences longer than [`LINEAR_SPACE_CELLS`]
/// cells switch to the linear-space algorithm automatically.
pub fn global(
    a: &[u8],
    b: &[u8],
    m: &SubstMatrix,
    gap_open: f32,
    gap_extend: f32,
) -> (Vec<u8>, Vec<u8>, f32) {
    global_ends(a, b, m, gap_open, gap_extend, 1.0)
}

/// Global alignment where gap runs touching either end of the matrix are
/// charged `terminal_gap_factor` times the normal penalty. `0.0` gives free end
/// gaps (semi-global alignment), `1.0` reproduces [`global`].
pub fn global_ends(
    a: &[u8],
    b: &[u8],
    mat: &SubstMatrix,
    gap_open: f32,
    gap_extend: f32,
    terminal_gap_factor: f32,
) -> (Vec<u8>, Vec<u8>, f32) {
    let n = a.len();
    let m = b.len();
    let gp = GapProfile::uniform(n, m, gap_open, gap_extend, terminal_gap_factor);
    if n == 0 || m == 0 {
        let mut ops = Vec::with_capacity(n + m);
        ops.extend(std::iter::repeat_n(OP_DELETE, n));
        ops.extend(std::iter::repeat_n(OP_INSERT, m));
        let s = score_ops(a, b, mat, &gp, &ops);
        let (ra, rb) = rows_from_ops(a, b, &ops);
        return (ra, rb, s);
    }
    let p = Problem { a, b, mat, gp };
    let ops = ops_for(&p);
    let s = score_ops(a, b, mat, &p.gp, &ops);
    let (ra, rb) = rows_from_ops(a, b, &ops);
    (ra, rb, s)
}

fn ops_for(p: &Problem<'_>) -> Vec<u8> {
    let (n, m) = (p.a.len(), p.b.len());
    let cells = (n as u64 + 1) * (m as u64 + 1);
    if cells <= LINEAR_SPACE_CELLS {
        full_dp(&p.forward(0, n, 0, m), p.mat, false, false).1
    } else {
        let mut ops = Vec::with_capacity(n + m);
        linear_path(p, 0, n, 0, m, false, false, &mut ops);
        ops
    }
}

/// Global alignment as a traceback, for callers that want to apply the same
/// edits to more than the two rows.
///
/// The profile aligner uses this when both "profiles" are single ungapped
/// sequences and the matrix is too big for its own quadratic traceback: the
/// linear-space path here needs O(n + m) memory instead of O(n * m), which is
/// what makes aligning a pair of organellar genomes possible at all.
pub(crate) fn global_ops(
    a: &[u8],
    b: &[u8],
    mat: &SubstMatrix,
    gap_open: f32,
    gap_extend: f32,
    terminal_gap_factor: f32,
) -> Vec<u8> {
    let (n, m) = (a.len(), b.len());
    if n == 0 || m == 0 {
        let mut ops = Vec::with_capacity(n + m);
        ops.extend(std::iter::repeat_n(OP_DELETE, n));
        ops.extend(std::iter::repeat_n(OP_INSERT, m));
        return ops;
    }
    let gp = GapProfile::uniform(n, m, gap_open, gap_extend, terminal_gap_factor);
    ops_for(&Problem { a, b, mat, gp })
}

/// Local alignment (Smith-Waterman) with affine gaps.
///
/// Returns the two gapped rows of the best local alignment, its score, and the
/// half-open residue ranges it covers in `a` and in `b`. When no positively
/// scoring alignment exists the rows and ranges are empty and the score is 0.
///
/// This path is always quadratic in memory: local alignment is used on short
/// pairs (seeded hits, "find this motif"), not on whole genomes.
pub fn local(
    a: &[u8],
    b: &[u8],
    mat: &SubstMatrix,
    gap_open: f32,
    gap_extend: f32,
) -> (Vec<u8>, Vec<u8>, f32, Range<usize>, Range<usize>) {
    let n = a.len();
    let m = b.len();
    if n == 0 || m == 0 {
        return (Vec::new(), Vec::new(), 0.0, 0..0, 0..0);
    }
    let stride = m + 1;
    let mut tb = vec![0u8; stride * (n + 1)];
    let mut prev_m = vec![0.0f32; stride];
    let mut prev_x = vec![NEG; stride];
    let mut prev_y = vec![NEG; stride];
    let mut cur_m = vec![0.0f32; stride];
    let mut cur_x = vec![NEG; stride];
    let mut cur_y = vec![NEG; stride];
    let mut best = 0.0f32;
    let mut best_at = (0usize, 0usize);

    for i in 1..=n {
        let ai = a[i - 1];
        let row = i * stride;
        cur_m[0] = 0.0;
        cur_x[0] = NEG;
        cur_y[0] = NEG;
        for j in 1..=m {
            let idx = row + j;
            let (diag, dst) = argmax3(prev_m[j - 1], prev_x[j - 1], prev_y[j - 1], S_M, S_X, S_Y);
            let v = max2(diag, 0.0) + mat.score(ai, b[j - 1]);
            if v <= 0.0 || diag <= 0.0 {
                cur_m[j] = max2(v, 0.0);
                tb[idx] |= S_START;
            } else {
                cur_m[j] = v;
                tb[idx] |= dst;
            }

            let opened = max2(prev_m[j], prev_y[j]) - gap_open - gap_extend;
            let ost = if prev_m[j] >= prev_y[j] { S_M } else { S_Y };
            let extended = prev_x[j] - gap_extend;
            if extended > opened {
                cur_x[j] = extended;
                tb[idx] |= S_X << 2;
            } else {
                cur_x[j] = opened;
                tb[idx] |= ost << 2;
            }

            let lopened = max2(cur_m[j - 1], cur_x[j - 1]) - gap_open - gap_extend;
            let lst = if cur_m[j - 1] >= cur_x[j - 1] { S_M } else { S_X };
            let lextended = cur_y[j - 1] - gap_extend;
            if lextended > lopened {
                cur_y[j] = lextended;
                tb[idx] |= S_Y << 4;
            } else {
                cur_y[j] = lopened;
                tb[idx] |= lst << 4;
            }

            if cur_m[j] > best {
                best = cur_m[j];
                best_at = (i, j);
            }
        }
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_x, &mut cur_x);
        std::mem::swap(&mut prev_y, &mut cur_y);
    }

    if best <= 0.0 {
        return (Vec::new(), Vec::new(), 0.0, 0..0, 0..0);
    }

    let (ei, ej) = best_at;
    let (mut i, mut j) = (ei, ej);
    let mut state = S_M;
    let mut ops: Vec<u8> = Vec::new();
    while i > 0 && j > 0 {
        let idx = i * stride + j;
        match state {
            S_M => {
                let pred = tb[idx] & 0b11;
                ops.push(OP_MATCH);
                i -= 1;
                j -= 1;
                if pred == S_START {
                    break;
                }
                state = pred;
            }
            S_X => {
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
    let range_a = i..ei;
    let range_b = j..ej;
    let (ra, rb) = rows_from_ops(&a[range_a.clone()], &b[range_b.clone()], &ops);
    (ra, rb, best, range_a, range_b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::SubstMatrix;

    fn s(v: &[u8]) -> Vec<u8> {
        v.to_vec()
    }

    #[test]
    fn identical_sequences_align_without_gaps() {
        let m = SubstMatrix::identity();
        let (a, b, sc) = global(b"ACGTACGT", b"ACGTACGT", m, 5.0, 1.0);
        assert_eq!(a, s(b"ACGTACGT"));
        assert_eq!(b, s(b"ACGTACGT"));
        assert_eq!(sc, 8.0);
    }

    #[test]
    fn single_internal_gap_is_placed() {
        let m = SubstMatrix::identity();
        let (a, b, sc) = global(b"ACGT", b"ACT", m, 1.0, 0.5);
        assert_eq!(a, s(b"ACGT"));
        assert_eq!(b.len(), 4);
        assert_eq!(b.iter().filter(|&&c| c == b'-').count(), 1);
        // 3 matches minus one gap of length 1 = 3 - (1.0 + 0.5).
        assert!((sc - 1.5).abs() < 1e-5, "score {sc}");
    }

    #[test]
    fn affine_penalties_prefer_one_long_gap() {
        let m = SubstMatrix::identity();
        let a = b"AAAACCCC";
        let b = b"AAAAGGGGCCCC";
        let (ra, rb, sc) = global(a, b, m, 10.0, 0.5);
        assert_eq!(ra.len(), rb.len());
        // 8 matches minus one gap of 4 = 8 - (10 + 4*0.5) = -4.
        assert!((sc - (-4.0)).abs() < 1e-5, "score {sc}");
        let gaps: Vec<usize> =
            ra.iter().enumerate().filter(|(_, &c)| c == b'-').map(|(i, _)| i).collect();
        assert_eq!(gaps.len(), 4);
        assert_eq!(gaps[3] - gaps[0], 3, "the four gaps must be contiguous");
    }

    #[test]
    fn gap_cost_is_affine_not_linear() {
        let m = SubstMatrix::identity();
        let (open, ext) = (4.0f32, 1.0f32);
        let (_, _, s4) = global(b"AA", b"AGGGGA", m, open, ext);
        let (_, _, s2) = global(b"AA", b"AGGA", m, open, ext);
        // Doubling the gap length adds only 2 * ext, not 2 * (open + ext).
        assert!(((s2 - s4) - 2.0 * ext).abs() < 1e-5, "{s2} {s4}");
    }

    #[test]
    fn terminal_gaps_can_be_made_free() {
        let m = SubstMatrix::identity();
        let a = b"ACGTACGT";
        let b = b"TTTTACGTACGTTTTT";
        let (_, _, strict) = global_ends(a, b, m, 10.0, 1.0, 1.0);
        let (ra, rb, free) = global_ends(a, b, m, 10.0, 1.0, 0.0);
        assert!(free > strict);
        assert!((free - 8.0).abs() < 1e-5, "free-end score {free}");
        assert_eq!(ra.len(), rb.len());
        assert_eq!(rb, b.to_vec());
    }

    /// Deterministic pseudo-random generator: no `rand` dependency.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn random_pair(rng: &mut Rng, len: usize, alphabet: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let a: Vec<u8> = (0..len).map(|_| alphabet[rng.below(alphabet.len())]).collect();
        let mut b = Vec::with_capacity(len);
        for &c in &a {
            match rng.below(100) {
                0..=7 => {}
                8..=15 => {
                    b.push(alphabet[rng.below(alphabet.len())]);
                    b.push(c);
                }
                16..=29 => b.push(alphabet[rng.below(alphabet.len())]),
                _ => b.push(c),
            }
        }
        (a, b)
    }

    /// The linear-space path must return exactly the same score as the full DP.
    #[test]
    fn hirschberg_matches_full_dp() {
        let mat = SubstMatrix::blosum62();
        let mut rng = Rng(0x2545F4914F6CDD1D);
        for trial in 0..25 {
            let len = 40 + trial * 13;
            let (a, b) = random_pair(&mut rng, len, b"ACDEFGHIKLMNPQRSTVWY");
            for &(open, ext) in &[(10.0f32, 0.5f32), (4.0, 2.0), (1.0, 1.0)] {
                let gp = GapProfile::uniform(a.len(), b.len(), open, ext, 1.0);
                let p = Problem { a: &a, b: &b, mat, gp };
                let (full_score, full_ops) =
                    full_dp(&p.forward(0, a.len(), 0, b.len()), mat, false, false);
                let mut lin_ops = Vec::new();
                linear_path(&p, 0, a.len(), 0, b.len(), false, false, &mut lin_ops);
                let lin_score = score_ops(&a, &b, mat, &p.gp, &lin_ops);
                let full_check = score_ops(&a, &b, mat, &p.gp, &full_ops);
                assert!(
                    (full_score - full_check).abs() < 1e-2,
                    "full DP score {full_score} disagrees with rescoring {full_check}"
                );
                assert!(
                    (full_score - lin_score).abs() < 1e-2,
                    "trial {trial} open {open}: full {full_score} vs linear {lin_score}"
                );
            }
        }
    }

    /// Same, with discounted end gaps, exercising the position-dependent
    /// terminal penalties through the recursion.
    #[test]
    fn hirschberg_matches_full_dp_with_free_ends() {
        let mat = SubstMatrix::identity();
        let mut rng = Rng(0x9E3779B97F4A7C15);
        for trial in 0..20 {
            let len = 30 + trial * 11;
            let (mut a, mut b) = random_pair(&mut rng, len, b"ACGT");
            let pad = trial % 7;
            a.splice(0..0, std::iter::repeat_n(b'A', pad));
            b.extend(std::iter::repeat_n(b'T', pad));
            for &tgf in &[0.0f32, 0.5] {
                let gp = GapProfile::uniform(a.len(), b.len(), 6.0, 0.5, tgf);
                let p = Problem { a: &a, b: &b, mat, gp };
                let (full_score, _) =
                    full_dp(&p.forward(0, a.len(), 0, b.len()), mat, false, false);
                let mut lin_ops = Vec::new();
                linear_path(&p, 0, a.len(), 0, b.len(), false, false, &mut lin_ops);
                let lin_score = score_ops(&a, &b, mat, &p.gp, &lin_ops);
                assert!(
                    (full_score - lin_score).abs() < 1e-2,
                    "trial {trial} tgf {tgf}: full {full_score} vs linear {lin_score}"
                );
            }
        }
    }

    #[test]
    fn global_handles_empty_input() {
        let m = SubstMatrix::identity();
        let (a, b, _) = global(b"", b"ACGT", m, 5.0, 1.0);
        assert_eq!(a, s(b"----"));
        assert_eq!(b, s(b"ACGT"));
        let (a, b, _) = global(b"", b"", m, 5.0, 1.0);
        assert!(a.is_empty() && b.is_empty());
    }

    #[test]
    fn local_finds_the_embedded_match() {
        let m = SubstMatrix::identity();
        let a = b"TTTTTACGTACGTTTTTT";
        let b = b"GGGGACGTACGTGGGG";
        let (ra, rb, score, range_a, range_b) = local(a, b, m, 5.0, 1.0);
        assert_eq!(score, 8.0);
        assert_eq!(ra, s(b"ACGTACGT"));
        assert_eq!(rb, s(b"ACGTACGT"));
        assert_eq!(&a[range_a], b"ACGTACGT");
        assert_eq!(&b[range_b], b"ACGTACGT");
    }

    #[test]
    fn local_on_unrelated_sequences_is_not_negative() {
        let m = SubstMatrix::blosum62();
        let (_, _, score, _, _) = local(b"WWWWWW", b"PPPPPP", m, 10.0, 1.0);
        assert!(score >= 0.0);
    }

    #[test]
    fn local_with_an_internal_gap() {
        let m = SubstMatrix::identity();
        // Cheap gaps, so bridging the TTT insertion beats reporting just one of
        // the two ACGT blocks: 8 matches - (0.5 + 3 * 0.1) = 7.2 > 4.
        let (ra, rb, score, _, _) = local(b"CCACGTTTTACGTCC", b"AAACGTACGTAA", m, 0.5, 0.1);
        assert!((score - 7.2).abs() < 1e-4, "score {score}");
        assert_eq!(ra.len(), rb.len());
        assert_eq!(ra.iter().filter(|&&c| c == b'-').count(), 0);
        assert_eq!(rb.iter().filter(|&&c| c == b'-').count(), 3);
    }
}
