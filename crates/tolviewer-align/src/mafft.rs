//! MAFFT-style engine (FFT-NS-2).
//!
//! Follows Katoh, Misawa, Kuma & Miyata (2002, *Nucleic Acids Res*
//! 30:3059-3066) and Katoh et al. (2005, *Nucleic Acids Res* 33:511-518):
//!
//! 1. residues are mapped to a two-component vector (normalised volume and
//!    polarity for protein, a two-component indicator for nucleotides);
//! 2. the two groups are correlated with an FFT, which turns "find the offsets
//!    at which these two sequences look alike" into one O(N log N) pass instead
//!    of an O(nm) scan;
//! 3. the strongest offsets are turned into ungapped homologous segments, the
//!    segments are chained, and the group-to-group DP is restricted to a band
//!    around the chain;
//! 4. the whole progressive alignment is done twice, the second time on a tree
//!    re-estimated from the first alignment - that is the "-2" in FFT-NS-2.
//!
//! Step 3 is a pure accelerator. If no confident segment is found the band is
//! dropped and the unrestricted DP runs, so correctness never depends on the
//! heuristic firing.

use tolviewer_core::{Alphabet, Result};

use crate::distance::{self, DistanceMethod};
use crate::fft;
use crate::profile::{self, AlignCtx, Band, Profile};
use crate::{tree, AlignParams, Progress};

/// Below this many DP cells the unrestricted DP is cheap enough that banding
/// is not worth its risk.
///
/// This is set deliberately high. Measured on the 200 x 1000 benchmark in
/// `tests/accuracy.rs`, banding every merge cut the run from 5.8 s to 3.2 s but
/// cost 0.14 of SP accuracy, because the anchors are drawn from the *consensus*
/// of a profile and that consensus is unreliable once the profile holds
/// hundreds of diverged sequences. Restricting the band to problems the exact
/// DP would struggle with - roughly 4000 x 4000 columns and up, i.e. organellar
/// genomes rather than gene families - keeps the speed where it is needed and
/// the accuracy everywhere else.
const MIN_FFT_CELLS: u64 = 16_000_000;
/// How many correlation peaks to turn into candidate diagonals.
const MAX_DIAGONALS: usize = 24;
/// Smallest half-width of the band kept around an anchor; larger profiles get
/// a proportionally wider one.
const MIN_ANCHOR_MARGIN: usize = 64;
/// Shortest segment accepted as an anchor.
const MIN_SEGMENT: usize = 12;
/// Segment score threshold, in multiples of the matrix's mean diagonal.
const SEGMENT_THRESHOLD: f32 = 5.0;
/// Per-residue score subtracted while hunting for segments, as a fraction of
/// the matrix's mean diagonal.
///
/// Without it the sweep below cannot find *maximal* segments on a matrix whose
/// mismatch scores are non-negative - the IUB nucleotide matrix scores every
/// mismatch 0, so the running sum never falls back to zero and a single
/// "segment" swallows the whole diagonal. Subtracting a background makes a
/// mismatch genuinely cost something, which is what the maximal-segment sweep
/// assumes.
const SEGMENT_BACKGROUND: f32 = 0.4;
/// A chain must cover at least this fraction of the shorter profile before the
/// band is trusted; otherwise the unrestricted DP runs.
const MIN_COVERAGE: f32 = 0.2;

/// Grantham (1974) polarity and volume for the twenty standard residues, in
/// the order A R N D C Q E G H I L K M F P S T W Y V.
#[rustfmt::skip]
const GRANTHAM: [(u8, f32, f32); 20] = [
    (b'A', 8.1, 31.0),  (b'R', 10.5, 124.0), (b'N', 11.6, 56.0), (b'D', 13.0, 54.0),
    (b'C', 5.5, 55.0),  (b'Q', 10.5, 85.0),  (b'E', 12.3, 83.0), (b'G', 9.0, 3.0),
    (b'H', 10.4, 96.0), (b'I', 5.2, 111.0),  (b'L', 4.9, 111.0), (b'K', 11.3, 119.0),
    (b'M', 5.7, 105.0), (b'F', 5.2, 132.0),  (b'P', 8.0, 32.5),  (b'S', 9.2, 32.0),
    (b'T', 8.6, 61.0),  (b'W', 5.4, 170.0),  (b'Y', 6.2, 136.0), (b'V', 5.9, 84.0),
];

/// Two-component vector per residue slot (`A`..`Z`).
///
/// Protein uses volume and polarity, each centred and scaled to unit variance
/// over the twenty standard residues, exactly as Katoh et al. describe.
/// Nucleotides use the two-component indicator A = (1,0), G = (0,1),
/// C = (-1,0), T/U = (0,-1), whose correlation is positive only for identical
/// bases.
fn residue_vectors(alphabet: Alphabet) -> [[f64; 2]; 26] {
    let mut v = [[0.0f64; 2]; 26];
    if alphabet.is_nucleotide() {
        let set = |v: &mut [[f64; 2]; 26], c: u8, a: f64, b: f64| {
            v[(c - b'A') as usize] = [a, b];
        };
        set(&mut v, b'A', 1.0, 0.0);
        set(&mut v, b'G', 0.0, 1.0);
        set(&mut v, b'C', -1.0, 0.0);
        set(&mut v, b'T', 0.0, -1.0);
        set(&mut v, b'U', 0.0, -1.0);
        return v;
    }
    let pol_mean = GRANTHAM.iter().map(|g| g.1 as f64).sum::<f64>() / 20.0;
    let vol_mean = GRANTHAM.iter().map(|g| g.2 as f64).sum::<f64>() / 20.0;
    let pol_sd = (GRANTHAM.iter().map(|g| (g.1 as f64 - pol_mean).powi(2)).sum::<f64>() / 20.0)
        .sqrt()
        .max(1e-9);
    let vol_sd = (GRANTHAM.iter().map(|g| (g.2 as f64 - vol_mean).powi(2)).sum::<f64>() / 20.0)
        .sqrt()
        .max(1e-9);
    for &(c, pol, vol) in &GRANTHAM {
        v[(c - b'A') as usize] =
            [(vol as f64 - vol_mean) / vol_sd, (pol as f64 - pol_mean) / pol_sd];
    }
    v
}

/// Column-by-column two-component signal for a profile.
fn signal(p: &Profile, vecs: &[[f64; 2]; 26]) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0f64; p.width];
    let mut y = vec![0.0f64; p.width];
    for c in 0..p.width {
        for (s, f) in p.column_freq(c) {
            let v = vecs[s as usize];
            x[c] += f as f64 * v[0];
            y[c] += f as f64 * v[1];
        }
    }
    (x, y)
}

/// One ungapped homologous segment on the diagonal `j - i = d`.
#[derive(Debug, Clone, Copy)]
struct Segment {
    i0: usize,
    i1: usize, // exclusive
    d: isize,
    score: f32,
}

impl Segment {
    fn j0(&self) -> usize {
        (self.i0 as isize + self.d) as usize
    }
    fn j1(&self) -> usize {
        (self.i1 as isize + self.d) as usize
    }
}

/// Diagonals worth examining, strongest correlation first.
fn candidate_diagonals(a: &Profile, b: &Profile, alphabet: Alphabet) -> Vec<isize> {
    let (n, m) = (a.width, b.width);
    let vecs = residue_vectors(alphabet);
    let (ax, ay) = signal(a, &vecs);
    let (bx, by) = signal(b, &vecs);
    let size = n + m + 1;
    let cx = fft::cross_correlation(&ax, &bx, size);
    let cy = fft::cross_correlation(&ay, &by, size);
    let big = cx.len();

    let mut scored: Vec<(f32, isize)> = Vec::with_capacity(n + m);
    for k in 0..big {
        // Lags above m are the negative side of the circular correlation.
        let d: isize = if k <= m { k as isize } else { k as isize - big as isize };
        if d <= -(n as isize) || d >= m as isize {
            continue;
        }
        let lo = (-d).max(0) as usize;
        let hi = (m as isize - d).min(n as isize).max(0) as usize;
        if hi <= lo || hi - lo < MIN_SEGMENT {
            continue;
        }
        let overlap = (hi - lo) as f64;
        let v = (cx[k] + cy[k]) / overlap.sqrt();
        scored.push((v as f32, d));
    }
    scored.sort_by(|p, q| q.0.partial_cmp(&p.0).unwrap_or(std::cmp::Ordering::Equal));

    // Keep the strongest peaks, spreading them out so twenty adjacent lags on
    // one diagonal do not crowd out a second, weaker homologous region.
    let mut chosen: Vec<isize> = Vec::new();
    for (_, d) in scored {
        if chosen.len() >= MAX_DIAGONALS {
            break;
        }
        if chosen.iter().any(|&c| (c - d).abs() < 4) {
            continue;
        }
        chosen.push(d);
    }
    chosen
}

/// Maximal high-scoring ungapped segments on one diagonal, found with the
/// running-sum sweep used for maximal-segment-pair detection.
fn segments_on(ca: &[u8], cb: &[u8], d: isize, ctx: &AlignCtx, out: &mut Vec<Segment>) {
    let n = ca.len();
    let m = cb.len();
    let lo = (-d).max(0) as usize;
    let hi = ((m as isize - d).min(n as isize)).max(0) as usize;
    if hi <= lo {
        return;
    }
    let threshold = SEGMENT_THRESHOLD * ctx.mat.mean_diagonal();
    let background = SEGMENT_BACKGROUND * ctx.mat.mean_diagonal();
    let mut sum = 0.0f32;
    let mut start = lo;
    let mut best = 0.0f32;
    let mut best_end = lo;
    for i in lo..hi {
        let s = ctx.mat.score(ca[i], cb[(i as isize + d) as usize]) - background;
        sum += s;
        if sum > best {
            best = sum;
            best_end = i + 1;
        }
        if sum <= 0.0 {
            if best >= threshold && best_end > start && best_end - start >= MIN_SEGMENT {
                out.push(Segment { i0: start, i1: best_end, d, score: best });
            }
            sum = 0.0;
            best = 0.0;
            start = i + 1;
            best_end = start;
        }
    }
    if best >= threshold && best_end > start && best_end - start >= MIN_SEGMENT {
        out.push(Segment { i0: start, i1: best_end, d, score: best });
    }
}

/// Highest-scoring chain of mutually compatible segments.
fn chain(mut segs: Vec<Segment>) -> Vec<Segment> {
    if segs.is_empty() {
        return segs;
    }
    segs.sort_by_key(|s| (s.i0, s.j0()));
    let k = segs.len();
    let mut best = vec![0.0f32; k];
    let mut from = vec![usize::MAX; k];
    let mut top = 0usize;
    for i in 0..k {
        best[i] = segs[i].score;
        for j in 0..i {
            if segs[j].i1 <= segs[i].i0 && segs[j].j1() <= segs[i].j0() {
                let cand = best[j] + segs[i].score;
                if cand > best[i] {
                    best[i] = cand;
                    from[i] = j;
                }
            }
        }
        if best[i] > best[top] {
            top = i;
        }
    }
    let mut out = Vec::new();
    let mut cur = top;
    while cur != usize::MAX {
        out.push(segs[cur]);
        cur = from[cur];
    }
    out.reverse();
    out
}

/// Band for a problem big enough to be worth restricting, or `None` to run the
/// exact DP. This is the entry point the profile aligner uses.
pub(crate) fn maybe_band(a: &Profile, b: &Profile, ctx: &AlignCtx) -> Option<Band> {
    let (n, m) = (a.width, b.width);
    if (n as u64) * (m as u64) < MIN_FFT_CELLS {
        return None;
    }
    fft_band(a, b, ctx)
}

/// Band restricting the group-to-group DP, or `None` when no confident
/// homologous segment was found and the DP should run unrestricted.
pub(crate) fn fft_band(a: &Profile, b: &Profile, ctx: &AlignCtx) -> Option<Band> {
    let (n, m) = (a.width, b.width);
    if n == 0 || m == 0 {
        return None;
    }
    let diagonals = candidate_diagonals(a, b, ctx.alphabet);
    if diagonals.is_empty() {
        return None;
    }
    let ca = a.consensus();
    let cb = b.consensus();
    let mut segs = Vec::new();
    for d in diagonals {
        segments_on(&ca, &cb, d, ctx, &mut segs);
    }
    let anchors = chain(segs);
    if anchors.is_empty() {
        return None;
    }
    // A chain that explains only a sliver of the shorter profile is not enough
    // to justify constraining the DP: fall back to the exact algorithm.
    let covered: usize = anchors.iter().map(|s| s.i1 - s.i0).sum();
    if (covered as f32) < MIN_COVERAGE * n.min(m) as f32 {
        return None;
    }

    // Between two anchors the alignment path cannot leave the rectangle they
    // span, so that rectangle is an exact bound; inside an anchor we keep a
    // margin because the anchor itself is a heuristic.
    let margin = MIN_ANCHOR_MARGIN.max(n.min(m) / 20);
    let mut lo = vec![0u32; n + 1];
    let mut hi = vec![m as u32; n + 1];
    let clampj = |j: isize| j.clamp(0, m as isize) as u32;

    let mut prev_i = 0usize;
    let mut prev_j = 0usize;
    for s in &anchors {
        for (l, h) in lo.iter_mut().zip(hi.iter_mut()).take(s.i0 + 1).skip(prev_i) {
            *l = clampj(prev_j as isize);
            *h = clampj(s.j0() as isize);
        }
        for i in s.i0..=s.i1.min(n) {
            lo[i] = clampj(i as isize + s.d - margin as isize);
            hi[i] = clampj(i as isize + s.d + margin as isize);
        }
        prev_i = s.i1.min(n);
        prev_j = s.j1().min(m);
    }
    for i in prev_i..=n {
        lo[i] = clampj(prev_j as isize);
        hi[i] = m as u32;
    }

    let mut band = Band { lo, hi };
    band.sanitise(n, m);
    Some(band)
}

/// Run the MAFFT-style engine.
pub(crate) fn align(
    seqs: &[Vec<u8>],
    params: &AlignParams,
    alphabet: Alphabet,
    ctx: &AlignCtx,
    progress: &dyn Progress,
) -> Result<Vec<Vec<u8>>> {
    let n = seqs.len();
    let mut ctx = ctx.clone();
    ctx.use_fft = true;

    // Pass 1: 6-mer distances, guide tree, progressive alignment.
    let d = distance::matrix(
        seqs,
        DistanceMethod::Kmer { k: distance::default_k(alphabet) },
        alphabet,
        progress,
    )?;
    let t = tree::build(&d, params.tree, progress)?;
    let w = profile::tree_weights(&t, n);
    let mut best = profile::progressive(seqs, &t, &w, &ctx, progress, "MAFFT: first pass")?;
    // Weights for the objective are fixed from the first pass so that the score
    // of two candidate alignments is directly comparable.
    let obj_weights = profile::henikoff_weights(&best);
    let mut best_score = profile::sp_score(&best, &obj_weights, &ctx);

    // Passes 2..: re-estimate the tree from the current alignment and align
    // again, keeping the result only when the objective improves.
    let passes = 1 + params.iterations.clamp(1, 4);
    for pass in 1..passes {
        let d2 = distance::matrix(&best, DistanceMethod::FromAlignment, alphabet, progress)?;
        let t2 = tree::build(&d2, params.tree, progress)?;
        let w2 = profile::tree_weights(&t2, n);
        let msg = format!("MAFFT: pass {}", pass + 1);
        let cand = profile::progressive(seqs, &t2, &w2, &ctx, progress, &msg)?;
        let cs = profile::sp_score(&cand, &obj_weights, &ctx);
        if cs > best_score {
            best = cand;
            best_score = cs;
        } else {
            break;
        }
    }
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::SubstMatrix;
    use crate::profile::Profile;

    fn ctx() -> AlignCtx {
        // The transition-weighted matrix scores mismatches negatively, so
        // unrelated flanks are not spuriously "alignable" and the tests below
        // exercise the band rather than the tie-breaking.
        AlignCtx::new(SubstMatrix::clustal_dna(), 10.0, 0.5, 1.0, Alphabet::Dna)
    }

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

    fn random_dna(rng: &mut Rng, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.below(4)]).collect()
    }

    #[test]
    fn nucleotide_vectors_correlate_only_for_identical_bases() {
        let v = residue_vectors(Alphabet::Dna);
        let dot = |a: u8, b: u8| {
            let (x, y) = (v[(a - b'A') as usize], v[(b - b'A') as usize]);
            x[0] * y[0] + x[1] * y[1]
        };
        for &c in b"ACGT" {
            assert!(dot(c, c) > 0.0);
        }
        assert_eq!(dot(b'A', b'G'), 0.0);
        assert!(dot(b'A', b'C') < 0.0);
    }

    #[test]
    fn protein_vectors_are_centred() {
        let v = residue_vectors(Alphabet::Protein);
        let sum0: f64 = GRANTHAM.iter().map(|g| v[(g.0 - b'A') as usize][0]).sum();
        let sum1: f64 = GRANTHAM.iter().map(|g| v[(g.0 - b'A') as usize][1]).sum();
        assert!(sum0.abs() < 1e-9, "{sum0}");
        assert!(sum1.abs() < 1e-9, "{sum1}");
        // Large and small residues sit on opposite sides of the volume axis.
        assert!(v[(b'W' - b'A') as usize][0] > 0.0);
        assert!(v[(b'G' - b'A') as usize][0] < 0.0);
    }

    #[test]
    fn fft_band_locates_a_shifted_homologous_region() {
        let c = ctx();
        let mut rng = Rng(12345);
        let core = random_dna(&mut rng, 900);
        let mut a = random_dna(&mut rng, 60);
        a.extend_from_slice(&core);
        let mut b = random_dna(&mut rng, 200);
        b.extend_from_slice(&core);
        let pa = Profile::new(vec![a.clone()], vec![0], vec![1.0], &c);
        let pb = Profile::new(vec![b.clone()], vec![1], vec![1.0], &c);
        let band = fft_band(&pa, &pb, &c).expect("a 900 bp identical core must be found");
        // The true path runs along j = i + 140 through the core.
        let i = 500usize;
        assert!(
            band.lo[i] as usize <= i + 140 && band.hi[i] as usize >= i + 140,
            "band at {i} is {}..{}, wanted {}",
            band.lo[i],
            band.hi[i],
            i + 140
        );
        // And it must be much cheaper than the full matrix.
        let cells: u64 =
            band.lo.iter().zip(band.hi.iter()).map(|(&l, &h)| (h - l) as u64 + 1).sum();
        assert!(cells * 4 < (a.len() as u64 + 1) * (b.len() as u64 + 1), "{cells} cells");
    }

    #[test]
    fn banded_result_matches_unrestricted_dp_on_a_clean_case() {
        let c = ctx();
        let mut rng = Rng(99);
        let core = random_dna(&mut rng, 900);
        let mut a = random_dna(&mut rng, 60);
        a.extend_from_slice(&core);
        let mut b = random_dna(&mut rng, 200);
        b.extend_from_slice(&core);
        let pa = Profile::new(vec![a.clone()], vec![0], vec![1.0], &c);
        let pb = Profile::new(vec![b.clone()], vec![1], vec![1.0], &c);
        let full = profile::align_profiles(&pa, &pb, &c, &Band::full(pa.width, pb.width));
        let mut band = fft_band(&pa, &pb, &c).expect("band");
        band.sanitise(pa.width, pb.width);
        let banded = profile::align_profiles(&pa, &pb, &c, &band);
        assert_eq!(full, banded);
    }

    #[test]
    fn unrelated_sequences_fall_back_to_the_full_dp() {
        let c = ctx();
        let mut rng = Rng(7);
        // Two independent random sequences share no long ungapped segment.
        let a = random_dna(&mut rng, 700);
        let b = random_dna(&mut rng, 700);
        let pa = Profile::new(vec![a], vec![0], vec![1.0], &c);
        let pb = Profile::new(vec![b], vec![1], vec![1.0], &c);
        // Either no band at all, or one that still contains the whole diagonal.
        if let Some(band) = fft_band(&pa, &pb, &c) {
            assert!(band.lo.len() == pa.width + 1);
        }
    }
}
