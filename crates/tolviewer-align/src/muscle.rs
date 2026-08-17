//! MUSCLE-style engine: draft, re-tree, then iterative refinement.
//!
//! Edgar (2004, *Nucleic Acids Res* 32:1792-1797) in three stages:
//!
//! * **Stage 1 (draft progressive)** - k-mer distances, UPGMA tree, progressive
//!   alignment. Cheap and rough.
//! * **Stage 2 (improved progressive)** - recompute the distances from the
//!   draft alignment with the Kimura correction, rebuild the tree, and align
//!   again, *reusing the alignment of every subtree whose topology did not
//!   change*. That reuse is what makes stage 2 nearly free when the draft tree
//!   was mostly right.
//! * **Stage 3 (refinement)** - tree-dependent restricted partitioning: delete
//!   one tree edge at a time, in order of decreasing distance from the root,
//!   re-align the two profiles the deletion produces, and keep the result only
//!   if the sum-of-pairs score improves.
//!
//! Every accept/reject decision in stages 2 and 3 uses the single
//! [`profile::sp_score`] function, so the objective cannot drift between the
//! step that proposes a change and the step that judges it.

use tolviewer_core::{Alphabet, Result, GAP};

use crate::distance::{self, DistanceMethod};
use crate::profile::{self, AlignCtx, Profile, SubtreeCache};
use crate::tree::{self, GuideTree};
use crate::{AlignParams, Progress, TreeMethod};

/// Run the MUSCLE-style engine.
pub(crate) fn align(
    seqs: &[Vec<u8>],
    params: &AlignParams,
    alphabet: Alphabet,
    ctx: &AlignCtx,
    progress: &dyn Progress,
) -> Result<Vec<Vec<u8>>> {
    let n = seqs.len();

    // ---- Stage 1: draft.
    let d1 = distance::matrix(
        seqs,
        DistanceMethod::Kmer { k: distance::default_k(alphabet) },
        alphabet,
        progress,
    )?;
    let t1 = tree::build(&d1, TreeMethod::Upgma, progress)?;
    let w1 = profile::tree_weights(&t1, n);
    let mut cache = SubtreeCache::enabled();
    let mut reused = 0usize;
    let draft = profile::progressive_cached(
        seqs,
        &t1,
        &w1,
        ctx,
        progress,
        "MUSCLE: draft alignment",
        &mut cache,
        &mut reused,
    )?;
    let mut rows = profile::to_rows(&draft, n);

    // Weights for the objective are fixed from here on, so that the score of
    // two candidate alignments is directly comparable.
    let weights = profile::henikoff_weights(&rows);
    let mut best_score = profile::sp_score(&rows, &weights, ctx);
    let mut best_tree = t1;

    // ---- Stage 2: re-estimate the tree from the draft and re-align.
    let d2 = distance::matrix(&rows, DistanceMethod::FromAlignment, alphabet, progress)?;
    let t2 = tree::build(&d2, TreeMethod::Upgma, progress)?;
    if t2 != best_tree {
        let w2 = profile::tree_weights(&t2, n);
        // `cache` still holds every subtree alignment from stage 1; subtrees
        // whose topology survived are taken straight from it.
        let mut reused2 = 0usize;
        let improved = profile::progressive_cached(
            seqs,
            &t2,
            &w2,
            ctx,
            progress,
            "MUSCLE: improved alignment",
            &mut cache,
            &mut reused2,
        )?;
        if !progress
            .tick(0.0, &format!("MUSCLE: reused {reused2} unchanged subtrees from the draft"))
        {
            return Err(tolviewer_core::Error::Cancelled);
        }
        let cand = profile::to_rows(&improved, n);
        let score = profile::sp_score(&cand, &weights, ctx);
        if score > best_score {
            rows = cand;
            best_score = score;
            best_tree = t2;
        }
    }

    // ---- Stage 3: tree-dependent restricted partitioning.
    if params.iterations > 0 && n > 2 {
        refine(&mut rows, &best_tree, &weights, ctx, params.iterations, &mut best_score, progress)?;
    }
    Ok(rows)
}

/// Every non-root node of the tree, as `(depth, leaves under it)`, deepest
/// first: deleting the edge above a node splits the alignment into that node's
/// leaves and everything else.
fn edges(t: &GuideTree) -> Vec<(usize, Vec<usize>)> {
    let mut out: Vec<(usize, Vec<usize>)> = Vec::new();
    fn walk(t: &GuideTree, depth: usize, out: &mut Vec<(usize, Vec<usize>)>) -> Vec<usize> {
        match t {
            GuideTree::Leaf(i) => {
                let set = vec![*i];
                if depth > 0 {
                    out.push((depth, set.clone()));
                }
                set
            }
            GuideTree::Node { left, right, .. } => {
                let mut set = walk(left, depth + 1, out);
                set.extend(walk(right, depth + 1, out));
                if depth > 0 {
                    out.push((depth, set.clone()));
                }
                set
            }
        }
    }
    walk(t, 0, &mut out);
    // Decreasing distance from the root.
    out.sort_by_key(|e| std::cmp::Reverse(e.0));
    out
}

/// Extract the rows for `ids`, dropping columns that are all gaps within the
/// selection, and wrap them in a profile.
fn sub_profile(rows: &[Vec<u8>], ids: &[usize], weights: &[f32], ctx: &AlignCtx) -> Profile {
    let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let keep: Vec<bool> = (0..width)
        .map(|c| {
            ids.iter().any(|&i| {
                rows.get(i)
                    .and_then(|r| r.get(c))
                    .is_some_and(|&ch| !tolviewer_core::alphabet::is_gap(ch))
            })
        })
        .collect();
    let sub: Vec<Vec<u8>> = ids
        .iter()
        .map(|&i| {
            let r = rows.get(i);
            (0..width)
                .filter(|&c| keep[c])
                .map(|c| r.and_then(|r| r.get(c)).copied().unwrap_or(GAP))
                .collect()
        })
        .collect();
    let w: Vec<f32> = ids.iter().map(|&i| weights.get(i).copied().unwrap_or(1.0)).collect();
    Profile::new(sub, ids.to_vec(), w, ctx)
}

/// Stage 3: for each tree edge, split, re-align and keep the better of the two.
fn refine(
    rows: &mut Vec<Vec<u8>>,
    tree: &GuideTree,
    weights: &[f32],
    ctx: &AlignCtx,
    iterations: usize,
    best_score: &mut f32,
    progress: &dyn Progress,
) -> Result<()> {
    let n = rows.len();
    let all: Vec<usize> = (0..n).collect();
    let ed = edges(tree);
    if ed.is_empty() {
        return Ok(());
    }
    let total = (ed.len() * iterations).max(1);
    let mut done = 0usize;

    for _round in 0..iterations {
        let mut improved_any = false;
        for (_, inside) in &ed {
            done += 1;
            let outside: Vec<usize> = all.iter().copied().filter(|i| !inside.contains(i)).collect();
            if inside.is_empty() || outside.is_empty() {
                continue;
            }
            let a = sub_profile(rows, inside, weights, ctx);
            let b = sub_profile(rows, &outside, weights, ctx);
            let ops = profile::align_profiles_auto(&a, &b, ctx);
            let merged = profile::merge(&a, &b, &ops, ctx);
            let cand = profile::to_rows(&merged, n);
            let score = profile::sp_score(&cand, weights, ctx);
            if score > *best_score + 1e-4 {
                *rows = cand;
                *best_score = score;
                improved_any = true;
            }
            if !progress.tick(done as f32 / total as f32, "MUSCLE: refining") {
                return Err(tolviewer_core::Error::Cancelled);
            }
        }
        if !improved_any {
            break; // converged
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree3() -> GuideTree {
        GuideTree::Node {
            left: Box::new(GuideTree::Node {
                left: Box::new(GuideTree::Leaf(0)),
                right: Box::new(GuideTree::Leaf(1)),
                left_len: 0.1,
                right_len: 0.1,
            }),
            right: Box::new(GuideTree::Leaf(2)),
            left_len: 0.2,
            right_len: 0.2,
        }
    }

    #[test]
    fn edges_are_deepest_first_and_exclude_the_root() {
        let e = edges(&tree3());
        // Three leaves plus one internal node below the root.
        assert_eq!(e.len(), 4);
        assert!(e.windows(2).all(|w| w[0].0 >= w[1].0));
        assert!(e.iter().any(|(_, s)| s == &vec![0, 1]));
        assert!(!e.iter().any(|(_, s)| s.len() == 3));
    }

    #[test]
    fn sub_profile_drops_columns_that_are_all_gaps_inside_the_selection() {
        let ctx =
            AlignCtx::new(crate::matrix::SubstMatrix::identity(), 4.0, 0.5, 1.0, Alphabet::Dna);
        let rows = vec![b"AC--GT".to_vec(), b"AC--GT".to_vec(), b"ACTTGT".to_vec()];
        let p = sub_profile(&rows, &[0, 1], &[1.0; 3], &ctx);
        assert_eq!(p.width, 4);
        assert_eq!(p.rows[0], b"ACGT".to_vec());
        let q = sub_profile(&rows, &[2], &[1.0; 3], &ctx);
        assert_eq!(q.rows[0], b"ACTTGT".to_vec());
    }
}
