//! Guide trees: neighbor joining, UPGMA and Newick output.

use tolviewer_core::{Error, Result};

use crate::distance::DistMatrix;
use crate::Progress;

/// Above this many taxa, neighbor joining's cubic behaviour stops being worth
/// the wait and the engines fall back to UPGMA.
pub const NJ_TAXA_LIMIT: usize = 2000;

/// Rooted binary guide tree over leaf indices.
#[derive(Debug, Clone, PartialEq)]
pub enum GuideTree {
    /// A single input sequence, identified by its index.
    Leaf(usize),
    /// An internal node with its two children and their branch lengths.
    Node { left: Box<GuideTree>, right: Box<GuideTree>, left_len: f32, right_len: f32 },
}

impl GuideTree {
    /// Leaf indices in left-to-right order.
    pub fn leaves(&self) -> Vec<usize> {
        let mut out = Vec::new();
        let mut stack = vec![self];
        // Push the right child first so the left subtree is popped, and hence
        // visited, first; the result is left-to-right leaf order.
        while let Some(node) = stack.pop() {
            match node {
                GuideTree::Leaf(i) => out.push(*i),
                GuideTree::Node { left, right, .. } => {
                    stack.push(right);
                    stack.push(left);
                }
            }
        }
        out
    }

    /// Number of leaves under this node.
    pub fn leaf_count(&self) -> usize {
        match self {
            GuideTree::Leaf(_) => 1,
            GuideTree::Node { left, right, .. } => left.leaf_count() + right.leaf_count(),
        }
    }

    /// Height of the tree in edges; useful for sizing progress reports.
    pub fn depth(&self) -> usize {
        match self {
            GuideTree::Leaf(_) => 0,
            GuideTree::Node { left, right, .. } => 1 + left.depth().max(right.depth()),
        }
    }

    /// Newick representation, terminated by `;`. Names are looked up by leaf
    /// index; characters Newick cannot carry are replaced by `_`.
    pub fn to_newick(&self, names: &[String]) -> String {
        let mut s = String::new();
        self.write_newick(names, &mut s);
        s.push(';');
        s
    }

    fn write_newick(&self, names: &[String], out: &mut String) {
        match self {
            GuideTree::Leaf(i) => {
                let name = names.get(*i).map(String::as_str).unwrap_or("");
                if name.is_empty() {
                    out.push_str(&format!("seq{i}"));
                } else {
                    out.extend(name.chars().map(|c| match c {
                        ' ' | '\t' | '(' | ')' | '[' | ']' | ':' | ';' | ',' | '\'' => '_',
                        other => other,
                    }));
                }
            }
            GuideTree::Node { left, right, left_len, right_len } => {
                out.push('(');
                left.write_newick(names, out);
                out.push_str(&format!(":{left_len:.5}"));
                out.push(',');
                right.write_newick(names, out);
                out.push_str(&format!(":{right_len:.5}"));
                out.push(')');
            }
        }
    }
}

/// Neighbor joining (Saitou & Nei 1987, *Mol Biol Evol* 4:406-425).
///
/// The result is rooted at the final join, which is what a progressive aligner
/// wants even though NJ itself produces an unrooted topology. Negative branch
/// lengths, which NJ can produce on non-additive data, are clamped to zero.
///
/// Cost is O(n^3): each of the `n - 2` iterations scans the O(n^2) Q matrix
/// once, with the divergence sums `r[i]` maintained incrementally rather than
/// recomputed inside the scan.
pub fn neighbor_joining(d: &DistMatrix) -> Result<GuideTree> {
    let n = d.len();
    if n == 0 {
        return Err(Error::algorithm("cannot build a guide tree from zero sequences"));
    }
    if n == 1 {
        return Ok(GuideTree::Leaf(0));
    }

    // Working copy as a dense square matrix; slot `i` is reused for the node
    // created by merging into it.
    let mut m: Vec<Vec<f32>> = (0..n).map(|i| (0..n).map(|j| d.get(i, j)).collect()).collect();
    let mut nodes: Vec<Option<GuideTree>> = (0..n).map(|i| Some(GuideTree::Leaf(i))).collect();
    let mut active: Vec<usize> = (0..n).collect();
    let mut r: Vec<f32> = (0..n).map(|i| active.iter().map(|&j| m[i][j]).sum::<f32>()).collect();

    while active.len() > 2 {
        let k = active.len();
        let denom = (k - 2) as f32;
        let mut best = f32::INFINITY;
        let mut best_pair = (active[0], active[1]);
        for (ai, &i) in active.iter().enumerate() {
            for &j in &active[ai + 1..] {
                let q = denom * m[i][j] - r[i] - r[j];
                if q < best {
                    best = q;
                    best_pair = (i, j);
                }
            }
        }
        let (i, j) = best_pair;
        let dij = m[i][j];
        let mut li = 0.5 * dij + (r[i] - r[j]) / (2.0 * denom);
        let mut lj = dij - li;
        if li < 0.0 {
            lj -= li;
            li = 0.0;
        }
        if lj < 0.0 {
            li -= lj;
            lj = 0.0;
        }
        let left = nodes[i].take().unwrap_or(GuideTree::Leaf(i));
        let right = nodes[j].take().unwrap_or(GuideTree::Leaf(j));
        nodes[i] = Some(GuideTree::Node {
            left: Box::new(left),
            right: Box::new(right),
            left_len: li.max(0.0),
            right_len: lj.max(0.0),
        });

        // New distances, and incremental maintenance of the divergence sums.
        active.retain(|&x| x != j);
        let mut new_ri = 0.0f32;
        for &x in &active {
            if x == i {
                continue;
            }
            let nd = 0.5 * (m[i][x] + m[j][x] - dij);
            r[x] += nd - m[i][x] - m[j][x];
            m[i][x] = nd;
            m[x][i] = nd;
            new_ri += nd;
        }
        r[i] = new_ri;
        m[i][i] = 0.0;
    }

    let (i, j) = (active[0], active[1]);
    let half = (m[i][j] / 2.0).max(0.0);
    let left = nodes[i].take().unwrap_or(GuideTree::Leaf(i));
    let right = nodes[j].take().unwrap_or(GuideTree::Leaf(j));
    Ok(GuideTree::Node {
        left: Box::new(left),
        right: Box::new(right),
        left_len: half,
        right_len: half,
    })
}

/// UPGMA (average linkage) clustering.
///
/// Branch lengths are the difference between the height of a node and the
/// height of its child, so the tree is ultrametric by construction.
pub fn upgma(d: &DistMatrix) -> Result<GuideTree> {
    let n = d.len();
    if n == 0 {
        return Err(Error::algorithm("cannot build a guide tree from zero sequences"));
    }
    if n == 1 {
        return Ok(GuideTree::Leaf(0));
    }
    let mut m: Vec<Vec<f32>> = (0..n).map(|i| (0..n).map(|j| d.get(i, j)).collect()).collect();
    let mut nodes: Vec<Option<GuideTree>> = (0..n).map(|i| Some(GuideTree::Leaf(i))).collect();
    let mut size: Vec<f32> = vec![1.0; n];
    let mut height: Vec<f32> = vec![0.0; n];
    let mut active: Vec<usize> = (0..n).collect();

    while active.len() > 1 {
        let mut best = f32::INFINITY;
        let mut best_pair = (active[0], active[1]);
        for (ai, &i) in active.iter().enumerate() {
            for &j in &active[ai + 1..] {
                if m[i][j] < best {
                    best = m[i][j];
                    best_pair = (i, j);
                }
            }
        }
        let (i, j) = best_pair;
        let h = best / 2.0;
        let left = nodes[i].take().unwrap_or(GuideTree::Leaf(i));
        let right = nodes[j].take().unwrap_or(GuideTree::Leaf(j));
        nodes[i] = Some(GuideTree::Node {
            left: Box::new(left),
            right: Box::new(right),
            left_len: (h - height[i]).max(0.0),
            right_len: (h - height[j]).max(0.0),
        });
        let (si, sj) = (size[i], size[j]);
        active.retain(|&x| x != j);
        for &x in &active {
            if x == i {
                continue;
            }
            let nd = (si * m[i][x] + sj * m[j][x]) / (si + sj);
            m[i][x] = nd;
            m[x][i] = nd;
        }
        size[i] = si + sj;
        height[i] = h;
    }
    nodes[active[0]].take().ok_or_else(|| Error::algorithm("UPGMA produced no root"))
}

/// Build a guide tree with the requested method, dropping to UPGMA when
/// neighbor joining would be too slow. Reports the substitution through
/// `progress`, so the GUI can tell the user why the tree is not what they asked
/// for.
pub(crate) fn build(
    d: &DistMatrix,
    method: crate::TreeMethod,
    progress: &dyn Progress,
) -> Result<GuideTree> {
    let want_nj = method == crate::TreeMethod::NeighborJoining;
    if want_nj && d.len() > NJ_TAXA_LIMIT {
        if !progress
            .tick(0.0, "too many sequences for neighbor joining; using UPGMA for the guide tree")
        {
            return Err(Error::Cancelled);
        }
        return upgma(d);
    }
    if want_nj {
        neighbor_joining(d)
    } else {
        upgma(d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(rows: &[&[f32]]) -> DistMatrix {
        DistMatrix::from_square(&rows.iter().map(|r| r.to_vec()).collect::<Vec<_>>())
    }

    fn clades(t: &GuideTree) -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        fn walk(t: &GuideTree, out: &mut Vec<Vec<usize>>) -> Vec<usize> {
            match t {
                GuideTree::Leaf(i) => vec![*i],
                GuideTree::Node { left, right, .. } => {
                    let mut l = walk(left, out);
                    let r = walk(right, out);
                    l.extend(r);
                    l.sort_unstable();
                    out.push(l.clone());
                    l
                }
            }
        }
        walk(t, &mut out);
        out.sort();
        out
    }

    /// The unrooted topology, as the set of non-trivial bipartitions. NJ is an
    /// unrooted method, so this - not the rooted clade set - is what a
    /// published NJ example pins down.
    fn splits(t: &GuideTree) -> Vec<Vec<usize>> {
        let mut all = t.leaves();
        all.sort_unstable();
        let n = all.len();
        let mut out: Vec<Vec<usize>> = Vec::new();
        for c in clades(t) {
            if c.len() < 2 || c.len() > n - 2 {
                continue;
            }
            let comp: Vec<usize> = all.iter().copied().filter(|x| !c.contains(x)).collect();
            // Record both halves so a caller can look up whichever side it
            // finds natural.
            out.push(c);
            out.push(comp);
        }
        out.sort();
        out.dedup();
        out
    }

    /// The worked example from Saitou & Nei's original description, as it
    /// appears in most textbooks: five taxa whose NJ tree is ((A,B),(D,E),C)
    /// with C attached to the internal edge.
    #[test]
    fn neighbor_joining_recovers_the_textbook_topology() {
        //      A    B    C    D    E
        let d = square(&[
            &[0.0, 5.0, 4.0, 7.0, 6.0],
            &[5.0, 0.0, 7.0, 10.0, 9.0],
            &[4.0, 7.0, 0.0, 7.0, 6.0],
            &[7.0, 10.0, 7.0, 0.0, 5.0],
            &[6.0, 9.0, 6.0, 5.0, 0.0],
        ]);
        let t = neighbor_joining(&d).expect("five taxa");
        let c = splits(&t);
        assert!(c.contains(&vec![0, 1]), "A and B must be neighbours: {c:?}");
        assert!(c.contains(&vec![3, 4]), "D and E must be neighbours: {c:?}");
        assert_eq!(t.leaf_count(), 5);
    }

    /// A second textbook matrix (Wikipedia's NJ example) with the known
    /// answer ((A,B),(C,(D,E))).
    #[test]
    fn neighbor_joining_second_textbook_case() {
        let d = square(&[
            &[0.0, 5.0, 9.0, 9.0, 8.0],
            &[5.0, 0.0, 10.0, 10.0, 9.0],
            &[9.0, 10.0, 0.0, 8.0, 7.0],
            &[9.0, 10.0, 8.0, 0.0, 3.0],
            &[8.0, 9.0, 7.0, 3.0, 0.0],
        ]);
        let t = neighbor_joining(&d).expect("five taxa");
        let c = splits(&t);
        assert!(c.contains(&vec![0, 1]), "{c:?}");
        assert!(c.contains(&vec![3, 4]), "{c:?}");
    }

    /// On an additive matrix NJ must reproduce the branch lengths exactly.
    #[test]
    fn neighbor_joining_recovers_additive_branch_lengths() {
        // Tree: ((A:1,B:2):1,(C:3,D:4):1) -- pairwise sums of path lengths.
        let d = square(&[
            &[0.0, 3.0, 6.0, 7.0],
            &[3.0, 0.0, 7.0, 8.0],
            &[6.0, 7.0, 0.0, 7.0],
            &[7.0, 8.0, 7.0, 0.0],
        ]);
        let t = neighbor_joining(&d).expect("four taxa");
        let c = splits(&t);
        // On four taxa there is exactly one non-trivial split: {A,B}|{C,D}.
        assert_eq!(c, vec![vec![0, 1], vec![2, 3]], "{c:?}");
        // The A/B split must carry lengths 1 and 2.
        fn find(t: &GuideTree) -> Option<(f32, f32)> {
            match t {
                GuideTree::Leaf(_) => None,
                GuideTree::Node { left, right, left_len, right_len } => {
                    if matches!(**left, GuideTree::Leaf(0)) && matches!(**right, GuideTree::Leaf(1))
                    {
                        return Some((*left_len, *right_len));
                    }
                    find(left).or_else(|| find(right))
                }
            }
        }
        let (a, b) = find(&t).expect("A/B cherry");
        assert!((a - 1.0).abs() < 1e-4 && (b - 2.0).abs() < 1e-4, "{a} {b}");
    }

    #[test]
    fn upgma_on_an_ultrametric_matrix() {
        // ((A,B):h=1, (C,D):h=1) joined at h=3.
        let d = square(&[
            &[0.0, 2.0, 6.0, 6.0],
            &[2.0, 0.0, 6.0, 6.0],
            &[6.0, 6.0, 0.0, 2.0],
            &[6.0, 6.0, 2.0, 0.0],
        ]);
        let t = upgma(&d).expect("four taxa");
        let c = clades(&t);
        assert!(c.contains(&vec![0, 1]), "{c:?}");
        assert!(c.contains(&vec![2, 3]), "{c:?}");
        // Root branch lengths: 3 - 1 = 2 on both sides.
        if let GuideTree::Node { left_len, right_len, .. } = &t {
            assert!((left_len - 2.0).abs() < 1e-5, "{left_len}");
            assert!((right_len - 2.0).abs() < 1e-5, "{right_len}");
        } else {
            panic!("root must be internal");
        }
    }

    #[test]
    fn degenerate_sizes() {
        assert!(neighbor_joining(&DistMatrix::zeros(0)).is_err());
        assert!(upgma(&DistMatrix::zeros(0)).is_err());
        assert_eq!(neighbor_joining(&DistMatrix::zeros(1)).unwrap(), GuideTree::Leaf(0));
        assert_eq!(upgma(&DistMatrix::zeros(1)).unwrap(), GuideTree::Leaf(0));
        let t = neighbor_joining(&DistMatrix::zeros(2)).unwrap();
        assert_eq!(t.leaf_count(), 2);
    }

    #[test]
    fn leaves_covers_every_index_once() {
        let d = square(&[
            &[0.0, 5.0, 9.0, 9.0, 8.0],
            &[5.0, 0.0, 10.0, 10.0, 9.0],
            &[9.0, 10.0, 0.0, 8.0, 7.0],
            &[9.0, 10.0, 8.0, 0.0, 3.0],
            &[8.0, 9.0, 7.0, 3.0, 0.0],
        ]);
        let t = neighbor_joining(&d).unwrap();
        let mut l = t.leaves();
        assert_eq!(l.len(), 5);
        l.sort_unstable();
        assert_eq!(l, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn newick_round_trips_names() {
        let t = GuideTree::Node {
            left: Box::new(GuideTree::Leaf(0)),
            right: Box::new(GuideTree::Leaf(1)),
            left_len: 0.5,
            right_len: 0.25,
        };
        let names = vec!["Homo sapiens".to_string(), "Pan(troglodytes)".to_string()];
        let nw = t.to_newick(&names);
        assert_eq!(nw, "(Homo_sapiens:0.50000,Pan_troglodytes_:0.25000);");
    }
}
