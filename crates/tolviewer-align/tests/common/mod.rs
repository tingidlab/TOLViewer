//! Shared helpers for the integration tests: a deterministic PRNG, a sequence
//! simulator with a known true alignment, and the two standard accuracy
//! measures.
//!
//! Nothing here depends on `rand`; the generator is a plain xorshift64* so
//! every run of the test suite sees exactly the same data.

#![allow(dead_code)]

use tolviewer_core::{Alignment, Sequence};

/// xorshift64* - small, fast, and good enough to simulate sequence evolution.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(if seed == 0 { 0x9E37_79B9_7F4A_7C15 } else { seed })
    }
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Uniform in `[0, n)`.
    pub fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() % n as u64) as usize
        }
    }
    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }
    pub fn pick(&mut self, from: &[u8]) -> u8 {
        from[self.below(from.len())]
    }
}

/// One residue plus the ordering key that identifies its homologous column.
#[derive(Debug, Clone, Copy)]
struct Site {
    key: f64,
    residue: u8,
}

/// Substitution and indel rates for a single branch of the simulated tree.
#[derive(Debug, Clone, Copy)]
pub struct Rates {
    /// Per-site probability of a substitution to a uniformly chosen residue.
    pub substitution: f64,
    /// Per-site probability of starting a deletion.
    pub deletion: f64,
    /// Per-site probability of starting an insertion.
    pub insertion: f64,
    /// Maximum length of an indel, in residues.
    pub indel_length: usize,
}

impl Rates {
    /// Closely related sequences: easy for any aligner.
    pub fn low() -> Self {
        Rates { substitution: 0.03, deletion: 0.004, insertion: 0.004, indel_length: 3 }
    }
    /// Moderately diverged: roughly 55-65 % pairwise identity after a few
    /// branches, which is where progressive aligners start to differ.
    pub fn moderate() -> Self {
        Rates { substitution: 0.10, deletion: 0.012, insertion: 0.012, indel_length: 5 }
    }
}

fn evolve(seq: &[Site], rates: &Rates, alphabet: &[u8], rng: &mut Rng) -> Vec<Site> {
    let mut out: Vec<Site> = Vec::with_capacity(seq.len() + 8);
    let mut i = 0usize;
    while i < seq.len() {
        if rng.chance(rates.deletion) {
            i += 1 + rng.below(rates.indel_length);
            continue;
        }
        if rng.chance(rates.insertion) {
            let len = 1 + rng.below(rates.indel_length);
            let lo = out.last().map(|s| s.key).unwrap_or(seq[i].key - 1.0);
            let hi = seq[i].key;
            // A random sub-interval, so two lineages that independently insert
            // at the same place get different columns rather than being scored
            // as homologous.
            let span = (hi - lo) * (0.1 + 0.8 * rng.unit());
            for k in 0..len {
                let t = (k + 1) as f64 / (len + 1) as f64;
                out.push(Site { key: lo + span * t, residue: rng.pick(alphabet) });
            }
        }
        let mut s = seq[i];
        if rng.chance(rates.substitution) {
            s.residue = rng.pick(alphabet);
        }
        out.push(s);
        i += 1;
    }
    out
}

/// Simulate `n` sequences by evolving a random ancestor down a random binary
/// tree of the given depth, and return the *true* alignment implied by the
/// simulation.
pub fn simulate(
    n: usize,
    length: usize,
    rates: Rates,
    alphabet: &[u8],
    depth: usize,
    seed: u64,
) -> Alignment {
    let mut rng = Rng::new(seed);
    let root: Vec<Site> =
        (0..length).map(|i| Site { key: i as f64, residue: rng.pick(alphabet) }).collect();

    let mut lineages: Vec<Vec<Site>> = vec![root];
    for _ in 0..depth {
        if lineages.len() >= n {
            break;
        }
        let mut grown: Vec<Vec<Site>> = Vec::new();
        for l in &lineages {
            grown.push(evolve(l, &rates, alphabet, &mut rng));
            grown.push(evolve(l, &rates, alphabet, &mut rng));
        }
        lineages = grown;
    }
    while lineages.len() < n {
        let pick = rng.below(lineages.len());
        let child = evolve(&lineages[pick], &rates, alphabet, &mut rng);
        lineages.push(child);
    }
    lineages.truncate(n);

    let mut keys: Vec<f64> = lineages.iter().flat_map(|l| l.iter().map(|s| s.key)).collect();
    keys.sort_by(|a, b| a.partial_cmp(b).expect("keys are finite"));
    keys.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

    let sequences = lineages
        .iter()
        .enumerate()
        .map(|(i, l)| {
            let mut row = Vec::with_capacity(keys.len());
            let mut p = 0usize;
            for &k in &keys {
                if p < l.len() && (l[p].key - k).abs() < 1e-12 {
                    row.push(l[p].residue);
                    p += 1;
                } else {
                    row.push(b'-');
                }
            }
            Sequence::new(format!("t{i}"), row)
        })
        .collect();
    Alignment::new("simulated", sequences)
}

/// For each row, the alignment column of each of its residues.
fn residue_columns(aln: &Alignment) -> Vec<Vec<usize>> {
    aln.sequences
        .iter()
        .map(|s| {
            s.residues.iter().enumerate().filter(|(_, &c)| c != b'-').map(|(i, _)| i).collect()
        })
        .collect()
}

/// For each row and column, the residue index there (or `None` for a gap).
fn index_map(aln: &Alignment) -> Vec<Vec<Option<usize>>> {
    aln.sequences
        .iter()
        .map(|s| {
            let mut k = 0usize;
            s.residues
                .iter()
                .map(|&c| {
                    if c == b'-' {
                        None
                    } else {
                        let v = Some(k);
                        k += 1;
                        v
                    }
                })
                .collect()
        })
        .collect()
}

/// Occupants (row, residue index) of every column of `aln`.
fn occupants(aln: &Alignment) -> Vec<Vec<(usize, usize)>> {
    let map = index_map(aln);
    let width = aln.width();
    let mut out = vec![Vec::new(); width];
    for (r, row) in map.iter().enumerate() {
        for (c, cell) in row.iter().enumerate() {
            if let Some(i) = cell {
                out[c].push((r, *i));
            }
        }
    }
    out
}

/// Sum-of-pairs accuracy: the fraction of residue pairs aligned in the
/// reference that are also aligned in `test`. This is the standard "SP score"
/// used to grade aligners against a benchmark.
pub fn sp_accuracy(reference: &Alignment, test: &Alignment) -> f32 {
    let tc = residue_columns(test);
    let mut total = 0usize;
    let mut hit = 0usize;
    for col in occupants(reference) {
        for a in 0..col.len() {
            for b in (a + 1)..col.len() {
                let (r1, i1) = col[a];
                let (r2, i2) = col[b];
                if r1 >= tc.len() || r2 >= tc.len() {
                    continue;
                }
                total += 1;
                if tc[r1].get(i1) == tc[r2].get(i2) {
                    hit += 1;
                }
            }
        }
    }
    if total == 0 {
        1.0
    } else {
        hit as f32 / total as f32
    }
}

/// Column accuracy ("TC score"): the fraction of reference columns reproduced
/// exactly - same residues, in one test column, with nothing extra in it.
/// Columns holding a single residue make no alignment claim and are skipped.
pub fn column_accuracy(reference: &Alignment, test: &Alignment) -> f32 {
    let tc = residue_columns(test);
    let mut test_occupancy = vec![0usize; test.width() + 1];
    for row in &tc {
        for &c in row {
            test_occupancy[c] += 1;
        }
    }
    let mut good = 0usize;
    let mut counted = 0usize;
    for col in occupants(reference) {
        if col.len() < 2 {
            continue;
        }
        counted += 1;
        let first = tc.get(col[0].0).and_then(|r| r.get(col[0].1)).copied();
        let same = col.iter().all(|&(r, i)| tc.get(r).and_then(|row| row.get(i)).copied() == first);
        if same {
            if let Some(c) = first {
                if test_occupancy.get(c).copied().unwrap_or(0) == col.len() {
                    good += 1;
                }
            }
        }
    }
    if counted == 0 {
        1.0
    } else {
        good as f32 / counted as f32
    }
}

/// Mean pairwise identity of an alignment, for reporting how hard a case is.
pub fn mean_identity(aln: &Alignment) -> f32 {
    let n = aln.len();
    let mut sum = 0.0f32;
    let mut k = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if let Some(v) = tolviewer_core::stats::pairwise_identity(
                &aln.sequences[i].residues,
                &aln.sequences[j].residues,
            ) {
                sum += v;
                k += 1;
            }
        }
    }
    if k == 0 {
        1.0
    } else {
        sum / k as f32
    }
}

/// Strip the gaps from an alignment, giving the aligner's input.
pub fn degapped(aln: &Alignment) -> Alignment {
    let mut out = aln.clone();
    out.degap();
    out
}

/// Build a small reference alignment written inline as `(id, row)` pairs.
pub fn reference(rows: &[(&str, &str)]) -> Alignment {
    Alignment::new(
        "reference",
        rows.iter().map(|(id, r)| Sequence::new(*id, r.as_bytes().to_vec())).collect(),
    )
}
