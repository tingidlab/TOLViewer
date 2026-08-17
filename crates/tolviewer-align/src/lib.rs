//! Multiple sequence alignment engines for TOLViewer.
//!
//! Three engines are provided, all native Rust reimplementations from the
//! published descriptions of the algorithms:
//!
//! * [`Engine::Clustal`] - ClustalW-style progressive alignment. The accuracy
//!   baseline: full pairwise distances, neighbor-joining guide tree,
//!   tree-derived sequence weights, position-specific gap penalties.
//! * [`Engine::Muscle`] - MUSCLE-style draft, tree re-estimation with subtree
//!   reuse, then iterative refinement by tree-dependent restricted
//!   partitioning.
//! * [`Engine::Mafft`] - MAFFT FFT-NS-2-style: homologous segments located by
//!   FFT correlation of volume/polarity signals, a banded group-to-group DP,
//!   and two progressive passes.
//!
//! Everything long-running takes a [`Progress`] and returns
//! [`tolviewer_core::Error::Cancelled`] when asked to stop.
//!
//! ```no_run
//! use tolviewer_align::{align, AlignParams, Engine, NoProgress};
//! # fn main() -> tolviewer_core::Result<()> {
//! # let input = tolviewer_core::Alignment::default();
//! let params = AlignParams::for_engine(Engine::Muscle);
//! let aligned = align(&input, &params, &NoProgress)?;
//! # Ok(()) }
//! ```
#![forbid(unsafe_code)]

pub mod distance;
pub mod fft;
pub mod matrix;
pub mod pairwise;
pub mod tree;

mod clustal;
mod mafft;
mod muscle;
mod profile;

use std::ops::Range;

use tolviewer_core::alphabet::is_gap;
use tolviewer_core::{Alignment, Alphabet, Error, Result, Sequence, GAP};

use crate::matrix::SubstMatrix;
use crate::profile::AlignCtx;

/// Progress reporting and cancellation for long operations.
///
/// Implementations must be cheap and thread-safe: `tick` is called often, and
/// although the alignment engines only call it from the thread that entered
/// the crate, [`distance::matrix`] and the engines share one implementation.
pub trait Progress: Sync {
    /// `fraction` is 0.0..=1.0. Return `false` to request cancellation; the
    /// algorithm then unwinds and returns [`Error::Cancelled`].
    fn tick(&self, fraction: f32, message: &str) -> bool;
}

/// A [`Progress`] that never cancels and discards messages.
pub struct NoProgress;

impl Progress for NoProgress {
    fn tick(&self, _fraction: f32, _message: &str) -> bool {
        true
    }
}

/// Which alignment algorithm to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Engine {
    /// Progressive alignment with a distance-based guide tree, in the style of
    /// ClustalW/Clustal Omega: pairwise distances -> NJ guide tree ->
    /// profile-profile progressive alignment with position-specific gap
    /// penalties.
    Clustal,
    /// Progressive draft then iterative refinement by tree-dependent restricted
    /// partitioning, in the style of MUSCLE: k-mer distance draft ->
    /// re-estimate tree from the draft -> re-align -> horizontal refinement.
    Muscle,
    /// FFT-accelerated group-to-group alignment in the style of MAFFT FFT-NS-2:
    /// residues mapped to volume/polarity vectors, homologous segments located
    /// by FFT correlation, progressive alignment, then a second pass on a tree
    /// re-estimated from the first alignment.
    Mafft,
}

impl Engine {
    /// Human-readable name for menus.
    pub fn name(self) -> &'static str {
        match self {
            Engine::Clustal => "Clustal",
            Engine::Muscle => "MUSCLE",
            Engine::Mafft => "MAFFT",
        }
    }

    /// All engines, for menu construction.
    pub fn all() -> &'static [Engine] {
        &[Engine::Clustal, Engine::Muscle, Engine::Mafft]
    }
}

/// Which substitution matrix to score with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixChoice {
    /// Pick by alphabet: IUB for nucleotide data, BLOSUM62 for protein.
    Auto,
    Blosum62,
    Blosum45,
    Blosum80,
    Pam250,
    Identity,
    /// ClustalW's IUB DNA matrix.
    Iub,
    /// ClustalW's transition-weighted DNA matrix.
    ClustalDna,
}

/// Guide tree construction method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeMethod {
    /// Neighbor joining (Saitou & Nei 1987). Handles unequal rates better.
    NeighborJoining,
    /// Average-linkage clustering. Cheap; what MUSCLE uses.
    Upgma,
}

/// Tuning for [`align`] and friends.
#[derive(Debug, Clone)]
pub struct AlignParams {
    pub engine: Engine,
    pub matrix: MatrixChoice,
    /// Affine gap penalties, positive numbers (they are subtracted). A gap of
    /// length `L` costs `gap_open + L * gap_extend`.
    ///
    /// The numbers are expressed in BLOSUM62 units and rescaled internally to
    /// whichever matrix is in use, so `gap_open = 10` means the same strength
    /// of penalty for protein and for DNA.
    pub gap_open: f32,
    pub gap_extend: f32,
    /// Terminal gaps cost this fraction of the normal penalty (0.0 = free ends).
    pub terminal_gap_factor: f32,
    /// Refinement iterations (Muscle/Mafft); 0 disables refinement.
    pub iterations: usize,
    /// Guide tree method for the first pass.
    pub tree: TreeMethod,
    /// Worker threads; 0 = all available cores.
    pub threads: usize,
}

impl Default for AlignParams {
    fn default() -> Self {
        AlignParams::for_engine(Engine::Clustal)
    }
}

impl AlignParams {
    /// Sensible defaults for one engine.
    pub fn for_engine(engine: Engine) -> Self {
        AlignParams {
            engine,
            matrix: MatrixChoice::Auto,
            // Chosen by sweeping gap_open in 3..20 and gap_extend in 0.1..2
            // over the simulated benchmark in `tests/accuracy.rs`. The surface
            // is flat between 3 and 6; 5.0 sits in the middle of the plateau
            // and is high enough not to shred noisy real-world data, which the
            // simulated benchmark does not model.
            gap_open: 5.0,
            gap_extend: 1.0,
            // Ragged ends are the norm in real data, so terminal gaps are
            // charged at half rate by default rather than in full.
            terminal_gap_factor: 0.5,
            iterations: match engine {
                Engine::Clustal => 0,
                Engine::Muscle | Engine::Mafft => 2,
            },
            tree: match engine {
                Engine::Clustal => TreeMethod::NeighborJoining,
                Engine::Muscle | Engine::Mafft => TreeMethod::Upgma,
            },
            threads: 0,
        }
    }
}

/// Run `f` on a private rayon pool of `threads` workers, or on the global pool
/// when `threads == 0`. The global pool is never reconfigured, so an
/// application that already set it up keeps its settings.
fn in_pool<T: Send>(threads: usize, f: impl FnOnce() -> T + Send) -> T {
    if threads == 0 {
        return f();
    }
    match rayon::ThreadPoolBuilder::new().num_threads(threads).build() {
        Ok(pool) => pool.install(f),
        // A pool we cannot build is not a reason to fail the alignment.
        Err(_) => f(),
    }
}

/// Guess the alphabet without needing a mutable alignment.
fn alphabet_of(aln: &Alignment) -> Alphabet {
    aln.alphabet_hint().unwrap_or_else(|| {
        Alphabet::guess(
            aln.sequences.iter().take(100).flat_map(|s| s.residues.iter().take(5000).copied()),
        )
    })
}

/// Quality scores at the non-gap positions of a row.
fn ungapped_quality(seq: &Sequence) -> Option<Vec<u8>> {
    let q = seq.quality.as_ref()?;
    Some(
        seq.residues
            .iter()
            .enumerate()
            .filter(|(_, &c)| !is_gap(c))
            .map(|(i, _)| q.get(i).copied().unwrap_or(0))
            .collect(),
    )
}

/// Rebuild an alignment from gapped rows, keeping the input's names, order,
/// descriptions, hidden flags and (re-gapped) quality scores.
fn rebuild(aln: &Alignment, rows: Vec<Vec<u8>>, alphabet: Alphabet) -> Alignment {
    let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let sequences = aln
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut out = s.clone();
            let row = rows.get(i).cloned().unwrap_or_default();
            let uq = ungapped_quality(s);
            out.residues = row;
            out.residues.resize(width, GAP);
            out.quality = uq.map(|q| {
                let mut k = 0usize;
                out.residues
                    .iter()
                    .map(|&c| {
                        if is_gap(c) {
                            0
                        } else {
                            let v = q.get(k).copied().unwrap_or(0);
                            k += 1;
                            v
                        }
                    })
                    .collect()
            });
            out
        })
        .collect();
    let mut out = Alignment::new(aln.name.clone(), sequences);
    out.set_alphabet(alphabet);
    out
}

/// Drop columns that are gaps in every row.
fn drop_all_gap_columns(rows: &mut [Vec<u8>]) {
    let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    if width == 0 || rows.is_empty() {
        return;
    }
    let keep: Vec<bool> =
        (0..width).map(|c| rows.iter().any(|r| r.get(c).is_some_and(|&ch| !is_gap(ch)))).collect();
    if keep.iter().all(|&k| k) {
        return;
    }
    for r in rows.iter_mut() {
        let mut out = Vec::with_capacity(width);
        for (c, &k) in keep.iter().enumerate() {
            if k {
                out.push(r.get(c).copied().unwrap_or(GAP));
            }
        }
        *r = out;
    }
}

/// Build the shared alignment context from user parameters.
fn context(params: &AlignParams, alphabet: Alphabet) -> AlignCtx {
    let mat = SubstMatrix::choose(params.matrix, alphabet);
    AlignCtx::new(
        mat,
        params.gap_open.max(0.0),
        params.gap_extend.max(0.0),
        params.terminal_gap_factor,
        alphabet,
    )
}

/// Align a set of ungapped sequences with the chosen engine.
fn run_engine(
    seqs: &[Vec<u8>],
    params: &AlignParams,
    alphabet: Alphabet,
    progress: &dyn Progress,
) -> Result<Vec<Vec<u8>>> {
    let ctx = context(params, alphabet);
    match params.engine {
        Engine::Clustal => clustal::align(seqs, params, alphabet, &ctx, progress),
        Engine::Muscle => muscle::align(seqs, params, alphabet, &ctx, progress),
        Engine::Mafft => mafft::align(seqs, params, alphabet, &ctx, progress),
    }
}

/// Align the (possibly already gapped) sequences.
///
/// Existing gaps are stripped first. Row order, names, descriptions and hidden
/// flags are preserved, and the result is always rectangular - including for
/// the degenerate cases (no sequences, one sequence, duplicate sequences,
/// zero-length sequences, rows that are entirely gaps).
pub fn align(aln: &Alignment, params: &AlignParams, progress: &dyn Progress) -> Result<Alignment> {
    let alphabet = alphabet_of(aln);
    let n = aln.len();
    if n == 0 {
        let mut out = Alignment::new(aln.name.clone(), Vec::new());
        out.set_alphabet(alphabet);
        return Ok(out);
    }
    let seqs: Vec<Vec<u8>> = aln.sequences.iter().map(|s| s.ungapped()).collect();
    if n == 1 {
        return Ok(rebuild(aln, seqs, alphabet));
    }
    // Nothing to align if every row is empty.
    if seqs.iter().all(|s| s.is_empty()) {
        return Ok(rebuild(aln, seqs, alphabet));
    }

    let mut rows = in_pool(params.threads, || run_engine(&seqs, params, alphabet, progress))?;
    drop_all_gap_columns(&mut rows);

    // The engines must never lose or reorder residues; if one ever did, say so
    // instead of handing back a corrupt alignment.
    for (i, row) in rows.iter().enumerate() {
        let ungapped: Vec<u8> = row.iter().copied().filter(|&c| !is_gap(c)).collect();
        if ungapped != seqs[i] {
            return Err(Error::algorithm(format!(
                "internal error: row {i} lost residues during alignment"
            )));
        }
    }
    if !progress.tick(1.0, "alignment complete") {
        return Err(Error::Cancelled);
    }
    Ok(rebuild(aln, rows, alphabet))
}

/// Re-align only columns `cols` of an existing alignment, keeping the rest
/// fixed.
///
/// Columns outside `cols` come out byte-identical to the input, and every row
/// still holds the same residues in the same order; only the gapping inside
/// the selected block changes. Used by "realign selection" in the GUI.
pub fn realign_region(
    aln: &Alignment,
    cols: Range<usize>,
    params: &AlignParams,
    progress: &dyn Progress,
) -> Result<Alignment> {
    let alphabet = alphabet_of(aln);
    let mut padded = aln.clone();
    padded.pad_to_width();
    let width = padded.width();
    let start = cols.start.min(width);
    let end = cols.end.clamp(start, width);
    if padded.is_empty() || start == end {
        return Ok(rebuild(
            &padded,
            padded.sequences.iter().map(|s| s.residues.clone()).collect(),
            alphabet,
        ));
    }

    let block: Vec<Vec<u8>> = padded
        .sequences
        .iter()
        .map(|s| s.residues[start..end].iter().copied().filter(|&c| !is_gap(c)).collect())
        .collect();

    let mut aligned = if block.iter().all(|b| b.is_empty()) || block.len() == 1 {
        block.clone()
    } else {
        in_pool(params.threads, || run_engine(&block, params, alphabet, progress))?
    };
    drop_all_gap_columns(&mut aligned);
    let bw = aligned.iter().map(|r| r.len()).max().unwrap_or(0);

    let rows: Vec<Vec<u8>> = padded
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut out = Vec::with_capacity(start + bw + (width - end));
            out.extend_from_slice(&s.residues[..start]);
            let mut mid = aligned.get(i).cloned().unwrap_or_default();
            mid.resize(bw, GAP);
            out.extend_from_slice(&mid);
            out.extend_from_slice(&s.residues[end..]);
            out
        })
        .collect();
    Ok(rebuild(&padded, rows, alphabet))
}

/// Align `query` sequences onto an existing alignment (profile-sequence).
///
/// The profile's own rows keep their relative alignment exactly: the only
/// change to them is that whole gap columns may be inserted where a query
/// sequence has residues the profile does not. Deleting the columns that are
/// gaps in every profile row therefore recovers the input profile unchanged.
/// Row order is profile rows first, then the queries in the order given.
pub fn add_to_alignment(
    profile: &Alignment,
    query: &[Sequence],
    params: &AlignParams,
    progress: &dyn Progress,
) -> Result<Alignment> {
    let alphabet = alphabet_of(profile);
    let mut base = profile.clone();
    base.pad_to_width();

    if query.is_empty() {
        return Ok(rebuild(
            &base,
            base.sequences.iter().map(|s| s.residues.clone()).collect(),
            alphabet,
        ));
    }
    let queries: Vec<Vec<u8>> = query.iter().map(|s| s.ungapped()).collect();

    let combined = Alignment::new(
        base.name.clone(),
        base.sequences.iter().cloned().chain(query.iter().cloned()).collect(),
    );

    if base.is_empty() {
        // Nothing to add to: this is a plain alignment of the queries.
        let mut rows = if queries.len() < 2 {
            queries.clone()
        } else {
            in_pool(params.threads, || run_engine(&queries, params, alphabet, progress))?
        };
        drop_all_gap_columns(&mut rows);
        return Ok(rebuild(&combined, rows, alphabet));
    }

    let ctx = context(params, alphabet);
    let base_rows: Vec<Vec<u8>> = base.sequences.iter().map(|s| s.residues.clone()).collect();
    let k = base_rows.len();

    // Align the queries among themselves first, so that a group of sequences
    // is added as one profile rather than one at a time.
    let mut query_rows = if queries.len() < 2 {
        queries.clone()
    } else {
        in_pool(params.threads, || run_engine(&queries, params, alphabet, progress))?
    };
    drop_all_gap_columns(&mut query_rows);

    let base_weights = profile::henikoff_weights(&base_rows);
    let query_weights = profile::henikoff_weights(&query_rows);
    let pa = profile::Profile::new(base_rows, (0..k).collect(), base_weights, &ctx);
    let pb =
        profile::Profile::new(query_rows, (k..k + queries.len()).collect(), query_weights, &ctx);
    let ops = in_pool(params.threads, || profile::align_profiles_auto(&pa, &pb, &ctx));
    let merged = profile::merge(&pa, &pb, &ops, &ctx);
    let rows = profile::to_rows(&merged, k + queries.len());
    if !progress.tick(1.0, "sequences added") {
        return Err(Error::Cancelled);
    }
    Ok(rebuild(&combined, rows, alphabet))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(rows: &[(&str, &[u8])]) -> Alignment {
        Alignment::new("test", rows.iter().map(|(id, r)| Sequence::new(*id, r.to_vec())).collect())
    }

    #[test]
    fn engine_names_and_list() {
        assert_eq!(Engine::all().len(), 3);
        assert_eq!(Engine::Muscle.name(), "MUSCLE");
    }

    #[test]
    fn defaults_differ_per_engine() {
        assert_eq!(AlignParams::for_engine(Engine::Clustal).iterations, 0);
        assert_eq!(AlignParams::for_engine(Engine::Muscle).iterations, 2);
        assert_eq!(AlignParams::default().engine, Engine::Clustal);
    }

    #[test]
    fn aligns_three_dna_sequences() {
        let a = aln(&[("a", b"ACGTACGTACGT"), ("b", b"ACGTTTACGTACGT"), ("c", b"ACGTACGTACGT")]);
        for &engine in Engine::all() {
            let p = AlignParams::for_engine(engine);
            let out = align(&a, &p, &NoProgress).expect("aligns");
            assert!(out.is_aligned(), "{} produced a ragged result", engine.name());
            assert_eq!(out.len(), 3);
            for (i, s) in out.sequences.iter().enumerate() {
                assert_eq!(s.id, a.sequences[i].id);
                assert_eq!(s.ungapped(), a.sequences[i].residues);
            }
        }
    }

    #[test]
    fn existing_gaps_are_stripped_first() {
        let a = aln(&[("a", b"AC--GTACGT"), ("b", b"ACGTACGT")]);
        let out = align(&a, &AlignParams::default(), &NoProgress).unwrap();
        assert!(out.is_aligned());
        assert_eq!(out.sequences[0].ungapped(), b"ACGTACGT".to_vec());
        assert_eq!(out.width(), 8);
    }

    #[test]
    fn quality_scores_follow_the_residues() {
        let mut a = aln(&[("a", b"ACGT"), ("b", b"ACGGT")]);
        a.sequences[0].quality = Some(vec![10, 20, 30, 40]);
        let out = align(&a, &AlignParams::default(), &NoProgress).unwrap();
        let q = out.sequences[0].quality.as_ref().expect("quality kept");
        assert_eq!(q.len(), out.sequences[0].residues.len());
        let kept: Vec<u8> = out.sequences[0]
            .residues
            .iter()
            .zip(q.iter())
            .filter(|(&c, _)| !is_gap(c))
            .map(|(_, &v)| v)
            .collect();
        assert_eq!(kept, vec![10, 20, 30, 40]);
    }
}
