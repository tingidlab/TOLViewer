//! ClustalW-style progressive alignment: the accuracy baseline.
//!
//! Thompson, Higgins & Gibson (1994, *Nucleic Acids Res* 22:4673-4680):
//! all-against-all pairwise distances, a neighbor-joining guide tree,
//! tree-derived sequence weights, then one pass of profile-profile alignment
//! up the tree with position-specific gap penalties.
//!
//! The only departure from the paper is the first step: computing `n(n-1)/2`
//! full pairwise alignments is quadratic in both sequence count and length, so
//! on large inputs the distances come from shared k-mer counts instead. The
//! switch is automatic and reported through `Progress`.

use tolviewer_core::{Alphabet, Result};

use crate::distance::{self, DistanceMethod};
use crate::profile::{self, AlignCtx};
use crate::{tree, AlignParams, Progress};

/// Above this many DP cells spent on the distance matrix, use k-mer distances.
const PAIRWISE_DISTANCE_BUDGET: u64 = 200_000_000;

/// Pick the distance method that will finish in reasonable time.
pub(crate) fn distance_method(seqs: &[Vec<u8>], alphabet: Alphabet) -> DistanceMethod {
    let n = seqs.len() as u64;
    let avg =
        if seqs.is_empty() { 0 } else { seqs.iter().map(|s| s.len()).sum::<usize>() / seqs.len() }
            as u64;
    let cells = n.saturating_mul(n.saturating_sub(1)) / 2 * avg.saturating_mul(avg);
    if cells > PAIRWISE_DISTANCE_BUDGET {
        DistanceMethod::Kmer { k: distance::default_k(alphabet) }
    } else {
        DistanceMethod::PairwiseAlignment
    }
}

/// Run the ClustalW-style engine.
pub(crate) fn align(
    seqs: &[Vec<u8>],
    params: &AlignParams,
    alphabet: Alphabet,
    ctx: &AlignCtx,
    progress: &dyn Progress,
) -> Result<Vec<Vec<u8>>> {
    let method = distance_method(seqs, alphabet);
    if matches!(method, DistanceMethod::Kmer { .. })
        && !progress.tick(0.0, "large input: using k-mer distances for the guide tree")
    {
        return Err(tolviewer_core::Error::Cancelled);
    }
    let d = distance::matrix(seqs, method, alphabet, progress)?;
    let t = tree::build(&d, params.tree, progress)?;
    let w = profile::tree_weights(&t, seqs.len());
    profile::progressive(seqs, &t, &w, ctx, progress, "Clustal: progressive alignment")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_inputs_use_full_pairwise_distances() {
        let seqs: Vec<Vec<u8>> = (0..10).map(|_| vec![b'A'; 300]).collect();
        assert_eq!(distance_method(&seqs, Alphabet::Dna), DistanceMethod::PairwiseAlignment);
    }

    #[test]
    fn large_inputs_fall_back_to_kmer_distances() {
        let seqs: Vec<Vec<u8>> = (0..500).map(|_| vec![b'A'; 2000]).collect();
        assert!(matches!(distance_method(&seqs, Alphabet::Dna), DistanceMethod::Kmer { .. }));
    }
}
