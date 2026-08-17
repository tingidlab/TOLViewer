//! Pairwise distances and distance matrices.
//!
//! Three families are provided:
//!
//! * observed divergence ([`p_distance`]) and the two corrections TOLViewer
//!   needs — Jukes & Cantor (1969) for nucleotides and Kimura's protein
//!   correction as used by ClustalW (Thompson et al. 1994, "Calculation of the
//!   guide tree");
//! * an alignment-free k-mer distance (Edgar 2004, *BMC Bioinformatics* 5:113,
//!   "k-mer distance measures"), which is what makes MUSCLE's first pass and
//!   MAFFT's guide tree cheap on large inputs;
//! * [`matrix`], which fills the whole lower triangle in parallel.

use rayon::prelude::*;
use tolviewer_core::alphabet::is_gap;
use tolviewer_core::{Alphabet, Error, Result};

use crate::matrix::SubstMatrix;
use crate::{pairwise, Progress};

/// Value returned by the corrections when the observed divergence is at or
/// beyond the point where the model has no information left.
const MAX_DISTANCE: f32 = 5.0;

/// How to fill a [`DistMatrix`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceMethod {
    /// Alignment-free shared-k-mer distance. Fast, approximate, and the only
    /// sane choice for thousands of sequences.
    Kmer { k: usize },
    /// Align every pair with [`pairwise::global`] and take the corrected
    /// divergence of that alignment. Accurate and quadratic in sequence length.
    PairwiseAlignment,
    /// The sequences are already aligned rows; take the corrected divergence
    /// column by column.
    FromAlignment,
}

/// A symmetric distance matrix, stored as its strict lower triangle.
#[derive(Debug, Clone, Default)]
pub struct DistMatrix {
    n: usize,
    data: Vec<f32>,
}

impl DistMatrix {
    /// An `n x n` matrix of zeros.
    pub fn zeros(n: usize) -> Self {
        DistMatrix { n, data: vec![0.0; n * n.saturating_sub(1) / 2] }
    }

    /// Number of taxa.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    #[inline]
    fn index(i: usize, j: usize) -> usize {
        // Strict lower triangle, row-major: (1,0), (2,0), (2,1), (3,0), ...
        i * (i - 1) / 2 + j
    }

    /// Distance between `i` and `j`; 0 on the diagonal. Out-of-range indices
    /// yield 0 rather than panicking, so a caller that lost track of `n`
    /// degrades instead of crashing.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f32 {
        if i == j || i >= self.n || j >= self.n {
            return 0.0;
        }
        let (i, j) = if i > j { (i, j) } else { (j, i) };
        self.data.get(Self::index(i, j)).copied().unwrap_or(0.0)
    }

    /// Set the distance between `i` and `j` (order does not matter).
    pub fn set(&mut self, i: usize, j: usize, v: f32) {
        if i == j || i >= self.n || j >= self.n {
            return;
        }
        let (i, j) = if i > j { (i, j) } else { (j, i) };
        let idx = Self::index(i, j);
        if let Some(slot) = self.data.get_mut(idx) {
            *slot = v;
        }
    }

    /// Build from a full square matrix, taking the lower triangle.
    pub fn from_square(rows: &[Vec<f32>]) -> Self {
        let n = rows.len();
        let mut d = DistMatrix::zeros(n);
        for (i, row) in rows.iter().enumerate() {
            for (j, &v) in row.iter().enumerate().take(i) {
                d.set(i, j, v);
            }
        }
        d
    }
}

/// Fraction of differing residues over positions where both sequences have a
/// residue. Comparison is case-insensitive. Returns 0.0 when the sequences do
/// not overlap at all.
pub fn p_distance(a: &[u8], b: &[u8]) -> f32 {
    let mut diff = 0usize;
    let mut total = 0usize;
    for (&x, &y) in a.iter().zip(b.iter()) {
        if is_gap(x) || is_gap(y) {
            continue;
        }
        total += 1;
        if !x.eq_ignore_ascii_case(&y) {
            diff += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        diff as f32 / total as f32
    }
}

/// Jukes-Cantor (1969) correction `d = -3/4 ln(1 - 4p/3)`, saturating at
/// [`MAX_DISTANCE`] once `p >= 0.75`, where the formula is undefined.
pub fn jukes_cantor(p: f32) -> f32 {
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 0.75 {
        return MAX_DISTANCE;
    }
    let d = -0.75 * (1.0 - 4.0 * p / 3.0).ln();
    d.clamp(0.0, MAX_DISTANCE)
}

/// Kimura's correction for protein distances, as used by ClustalW:
/// `d = -ln(1 - p - 0.2 p^2)`.
///
/// The closed form goes singular just below `p = 0.79`; ClustalW switches to a
/// table of empirical PAM values there. We instead hold the formula to
/// `p = 0.75` and interpolate linearly from that value up to
/// [`MAX_DISTANCE`] at `p = 1`, which keeps the function monotone and finite
/// without pretending to more precision than a 90 %-divergent pair supports.
pub fn kimura_protein(p: f32) -> f32 {
    if p <= 0.0 {
        return 0.0;
    }
    const KNEE: f32 = 0.75;
    if p < KNEE {
        let inner = 1.0 - p - 0.2 * p * p;
        if inner <= 0.0 {
            return MAX_DISTANCE;
        }
        return (-inner.ln()).clamp(0.0, MAX_DISTANCE);
    }
    let at_knee = -(1.0f32 - KNEE - 0.2 * KNEE * KNEE).ln();
    let t = ((p - KNEE) / (1.0 - KNEE)).clamp(0.0, 1.0);
    at_knee + t * (MAX_DISTANCE - at_knee)
}

/// Apply the correction that suits the alphabet.
pub(crate) fn correct(p: f32, alphabet: Alphabet) -> f32 {
    if alphabet.is_nucleotide() {
        jukes_cantor(p)
    } else {
        kimura_protein(p)
    }
}

/// Compressed protein alphabet used for k-mer counting.
///
/// Seven physico-chemical groups; MUSCLE uses a compressed alphabet of the
/// same flavour for exactly this purpose (Edgar 2004, "Compressed alphabets"),
/// because raw 20-letter k-mers are too sparse to detect remote similarity.
const PROTEIN_GROUPS: [&[u8]; 7] = [
    b"AGST", // small
    b"C",    // cysteine on its own
    b"DENQ", // acidic / amide
    b"FWY",  // aromatic
    b"HKR",  // basic
    b"ILMV", // aliphatic
    b"P",    // proline on its own
];

/// Map a residue to its alphabet code, or `None` for gaps and ambiguity codes.
fn residue_code(c: u8, alphabet: Alphabet) -> Option<u8> {
    let c = c.to_ascii_uppercase();
    if alphabet.is_nucleotide() {
        match c {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ => None,
        }
    } else {
        PROTEIN_GROUPS.iter().position(|g| g.contains(&c)).map(|p| p as u8)
    }
}

fn alphabet_size(alphabet: Alphabet) -> u64 {
    if alphabet.is_nucleotide() {
        4
    } else {
        PROTEIN_GROUPS.len() as u64
    }
}

/// Sorted `(kmer code, count)` pairs for one sequence.
#[derive(Debug, Clone, Default)]
pub(crate) struct KmerProfile {
    counts: Vec<(u32, u32)>,
    /// Number of k-mer positions that were countable.
    total: u32,
}

pub(crate) fn kmer_profile(seq: &[u8], k: usize, alphabet: Alphabet) -> KmerProfile {
    let base = alphabet_size(alphabet);
    if k == 0 || seq.len() < k || base.pow(k as u32) > u32::MAX as u64 {
        return KmerProfile::default();
    }
    let mut codes: Vec<u32> = Vec::with_capacity(seq.len());
    let mut code: u64 = 0;
    let mut run = 0usize; // consecutive codeable residues
    let modulus = base.pow(k as u32);
    for &c in seq {
        match residue_code(c, alphabet) {
            Some(v) => {
                code = (code * base + v as u64) % modulus;
                run += 1;
                if run >= k {
                    codes.push(code as u32);
                }
            }
            None => {
                run = 0;
                code = 0;
            }
        }
    }
    let total = codes.len() as u32;
    codes.sort_unstable();
    let mut counts: Vec<(u32, u32)> = Vec::new();
    for c in codes {
        match counts.last_mut() {
            Some(last) if last.0 == c => last.1 += 1,
            _ => counts.push((c, 1)),
        }
    }
    KmerProfile { counts, total }
}

/// Distance between two precomputed k-mer profiles.
pub(crate) fn kmer_distance_profiles(a: &KmerProfile, b: &KmerProfile) -> f32 {
    let denom = a.total.min(b.total);
    if denom == 0 {
        return 1.0;
    }
    let mut shared = 0u32;
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.counts.len() && j < b.counts.len() {
        match a.counts[i].0.cmp(&b.counts[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += a.counts[i].1.min(b.counts[j].1);
                i += 1;
                j += 1;
            }
        }
    }
    let f = shared as f32 / denom as f32;
    (1.0 - f).clamp(0.0, 1.0)
}

/// Fast alignment-free distance from shared k-mer counts.
///
/// `F` is the number of k-mers the two sequences share (counting
/// multiplicities) divided by the number of k-mer positions in the shorter
/// sequence; the distance is `1 - F`. Protein sequences are first mapped onto
/// a seven-letter compressed alphabet.
pub fn kmer_distance(a: &[u8], b: &[u8], k: usize, alphabet: Alphabet) -> f32 {
    let pa = kmer_profile(a, k, alphabet);
    let pb = kmer_profile(b, k, alphabet);
    kmer_distance_profiles(&pa, &pb)
}

/// A sensible k for the alphabet: 6 for nucleotides, 3 for the compressed
/// protein alphabet (7^3 = 343 buckets, dense enough to be informative on
/// sequences of a few hundred residues).
pub(crate) fn default_k(alphabet: Alphabet) -> usize {
    if alphabet.is_nucleotide() {
        6
    } else {
        3
    }
}

/// Distance of one pair under `method`.
fn pair_distance(
    a: &[u8],
    b: &[u8],
    method: DistanceMethod,
    alphabet: Alphabet,
    profiles: &[KmerProfile],
    i: usize,
    j: usize,
) -> f32 {
    match method {
        DistanceMethod::Kmer { .. } => kmer_distance_profiles(&profiles[i], &profiles[j]),
        DistanceMethod::FromAlignment => correct(p_distance(a, b), alphabet),
        DistanceMethod::PairwiseAlignment => {
            let mat = SubstMatrix::choose(crate::MatrixChoice::Auto, alphabet);
            let scale = mat.mean_diagonal() / SubstMatrix::blosum62().mean_diagonal();
            let (ga, gb, _) = pairwise::global(a, b, mat, 10.0 * scale, 0.5 * scale);
            correct(p_distance(&ga, &gb), alphabet)
        }
    }
}

/// Full pairwise distance matrix, lower triangle, computed in parallel.
///
/// `progress.tick` is called from the calling thread between chunks of work,
/// never from inside the parallel closure, so implementations do not have to
/// be re-entrant. Returning `false` from `tick` aborts with
/// [`Error::Cancelled`].
pub fn matrix(
    seqs: &[Vec<u8>],
    method: DistanceMethod,
    alphabet: Alphabet,
    progress: &dyn Progress,
) -> Result<DistMatrix> {
    let n = seqs.len();
    let mut out = DistMatrix::zeros(n);
    if n < 2 {
        return Ok(out);
    }

    let profiles: Vec<KmerProfile> = match method {
        DistanceMethod::Kmer { k } => {
            let k = if k == 0 { default_k(alphabet) } else { k };
            seqs.par_iter().map(|s| kmer_profile(s, k, alphabet)).collect()
        }
        _ => Vec::new(),
    };

    let pairs: Vec<(usize, usize)> = (1..n).flat_map(|i| (0..i).map(move |j| (i, j))).collect();
    let total = pairs.len();
    // Chunks large enough to amortise the rayon join, small enough that
    // cancellation feels immediate.
    let chunk = match method {
        DistanceMethod::PairwiseAlignment => 64,
        _ => 4096,
    }
    .min(total.max(1));

    let mut done = 0usize;
    for block in pairs.chunks(chunk) {
        let values: Vec<f32> = block
            .par_iter()
            .map(|&(i, j)| pair_distance(&seqs[i], &seqs[j], method, alphabet, &profiles, i, j))
            .collect();
        for (&(i, j), &v) in block.iter().zip(values.iter()) {
            out.set(i, j, v);
        }
        done += block.len();
        if !progress.tick(done as f32 / total as f32, "computing pairwise distances") {
            return Err(Error::Cancelled);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Counting(std::sync::atomic::AtomicUsize, usize);
    impl Progress for Counting {
        fn tick(&self, _f: f32, _m: &str) -> bool {
            self.0.fetch_add(1, std::sync::atomic::Ordering::SeqCst) < self.1
        }
    }

    #[test]
    fn p_distance_ignores_gaps() {
        assert_eq!(p_distance(b"ACGT", b"ACGT"), 0.0);
        assert_eq!(p_distance(b"ACGT", b"ACGA"), 0.25);
        assert_eq!(p_distance(b"AC-T", b"AC-A"), 1.0 / 3.0);
        assert_eq!(p_distance(b"acgt", b"ACGT"), 0.0);
        assert_eq!(p_distance(b"----", b"ACGT"), 0.0);
    }

    #[test]
    fn jukes_cantor_is_monotone_and_saturates() {
        assert_eq!(jukes_cantor(0.0), 0.0);
        assert!(jukes_cantor(0.1) > 0.1);
        assert!(jukes_cantor(0.5) > jukes_cantor(0.3));
        assert_eq!(jukes_cantor(0.75), MAX_DISTANCE);
        assert_eq!(jukes_cantor(0.9), MAX_DISTANCE);
        // Known value: p = 0.25 -> -0.75 ln(2/3) = 0.30409
        assert!((jukes_cantor(0.25) - 0.304_09).abs() < 1e-4);
    }

    #[test]
    fn kimura_protein_matches_the_closed_form_below_the_knee() {
        assert_eq!(kimura_protein(0.0), 0.0);
        let p = 0.3f32;
        let expected = -(1.0 - p - 0.2 * p * p).ln();
        assert!((kimura_protein(p) - expected).abs() < 1e-5);
        // Monotone across the knee and finite everywhere.
        let mut prev = 0.0;
        let mut x = 0.0;
        while x <= 1.0 {
            let d = kimura_protein(x);
            assert!(d.is_finite(), "kimura_protein({x}) = {d}");
            assert!(d >= prev - 1e-6, "not monotone at {x}");
            prev = d;
            x += 0.01;
        }
    }

    #[test]
    fn kmer_distance_is_zero_for_identical_and_high_for_unrelated() {
        let a = b"ACGTACGTTTGACCATTGACCA".to_vec();
        assert_eq!(kmer_distance(&a, &a, 4, Alphabet::Dna), 0.0);
        let b = b"GGGGGGGGGGGGGGGGGGGGGG".to_vec();
        assert!(kmer_distance(&a, &b, 4, Alphabet::Dna) > 0.9);
    }

    #[test]
    fn kmer_distance_uses_the_compressed_protein_alphabet() {
        // K and R are in the same group, so a K/R swap costs nothing.
        let a = b"MKVLWIPQKSTYAGGH".to_vec();
        let b = b"MRVLWIPQRSTYAGGH".to_vec();
        assert_eq!(kmer_distance(&a, &b, 3, Alphabet::Protein), 0.0);
        // ... but a swap across groups does.
        let c = b"MKVLWIPQKSTYAGGH".to_vec();
        let d = b"MKVLWIPQKSTYAGGP".to_vec();
        assert!(kmer_distance(&c, &d, 3, Alphabet::Protein) > 0.0);
    }

    #[test]
    fn kmer_distance_handles_short_sequences() {
        assert_eq!(kmer_distance(b"AC", b"ACGT", 6, Alphabet::Dna), 1.0);
        assert_eq!(kmer_distance(b"", b"", 6, Alphabet::Dna), 1.0);
    }

    #[test]
    fn matrix_is_symmetric_and_zero_on_the_diagonal() {
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTAAAA".to_vec(),
            b"ACGTACGTACGTAAAC".to_vec(),
            b"TTTTTTTTGGGGGGGG".to_vec(),
        ];
        let d = matrix(&seqs, DistanceMethod::Kmer { k: 4 }, Alphabet::Dna, &crate::NoProgress)
            .expect("no cancellation");
        assert_eq!(d.len(), 3);
        for i in 0..3 {
            assert_eq!(d.get(i, i), 0.0);
            for j in 0..3 {
                assert_eq!(d.get(i, j), d.get(j, i));
            }
        }
        assert!(d.get(0, 1) < d.get(0, 2));
    }

    #[test]
    fn matrix_respects_cancellation() {
        let seqs: Vec<Vec<u8>> =
            (0..40).map(|i| format!("ACGTACGTACGTACGT{i}").into_bytes()).collect();
        let p = Counting(std::sync::atomic::AtomicUsize::new(0), 0);
        let e = matrix(&seqs, DistanceMethod::PairwiseAlignment, Alphabet::Dna, &p);
        assert!(matches!(e, Err(Error::Cancelled)));
    }

    #[test]
    fn empty_and_single_sequence_matrices() {
        let d =
            matrix(&[], DistanceMethod::Kmer { k: 3 }, Alphabet::Dna, &crate::NoProgress).unwrap();
        assert!(d.is_empty());
        let d = matrix(
            &[b"ACGT".to_vec()],
            DistanceMethod::Kmer { k: 3 },
            Alphabet::Dna,
            &crate::NoProgress,
        )
        .unwrap();
        assert_eq!(d.len(), 1);
        assert_eq!(d.get(0, 0), 0.0);
    }
}
