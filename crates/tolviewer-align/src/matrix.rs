//! Substitution matrices, stored as a flat ASCII lookup table.
//!
//! Every matrix is expanded once, lazily, into a 128x128 table of `f32` so that
//! [`SubstMatrix::score`] is a single indexed load with no branching, case
//! folding or symbol-to-index mapping at scoring time. The tables are the NCBI
//! values (`ftp.ncbi.nlm.nih.gov/blast/matrices/`) including the `B`, `Z` and
//! `X` rows, plus ClustalW's two nucleotide matrices and a plain identity
//! matrix.

use std::sync::OnceLock;

use tolviewer_core::Alphabet;

use crate::MatrixChoice;

/// Number of distinct byte values indexed by the table. Bytes >= 128 are
/// folded into the range by masking, and land on entries that hold the
/// matrix's "unknown residue" score.
const SIZE: usize = 128;
const CELLS: usize = SIZE * SIZE;

/// Column/row order of the embedded NCBI tables.
const AA_ORDER: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

/// Nucleotide symbols used when expanding the DNA matrices.
const NT_ORDER: &[u8; 17] = b"ACGTURYSWKMBDHVN-";

/// A substitution matrix as a flat 128x128 ASCII lookup.
///
/// Construct one through the associated functions ([`SubstMatrix::blosum62`]
/// and friends); they return references to lazily built statics, so cloning is
/// never needed and `score` stays a single load.
pub struct SubstMatrix {
    name: &'static str,
    table: Box<[f32; CELLS]>,
    mean_diagonal: f32,
}

impl std::fmt::Debug for SubstMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SubstMatrix").field("name", &self.name).finish()
    }
}

impl SubstMatrix {
    /// Score for aligning residue `a` against residue `b`.
    ///
    /// Case-insensitive: lowercase input is folded when the table is built, so
    /// no work happens here. Unknown bytes score as the matrix's "unknown"
    /// entry (the `X` row for protein matrices, 0 elsewhere).
    #[inline]
    pub fn score(&self, a: u8, b: u8) -> f32 {
        let i = (a as usize) & (SIZE - 1);
        let j = (b as usize) & (SIZE - 1);
        // `i` and `j` are both < 128 after masking, so the index is < CELLS and
        // the bounds check folds away against the fixed-size array.
        self.table[(i << 7) | j]
    }

    /// Display name, e.g. `"BLOSUM62"`.
    pub fn name(&self) -> &str {
        self.name
    }

    /// Mean of the self-comparison scores over the matrix's core alphabet.
    ///
    /// Used to put gap penalties expressed in "BLOSUM62 units" onto the scale
    /// of whichever matrix is actually in use; see [`crate::AlignParams`].
    pub fn mean_diagonal(&self) -> f32 {
        self.mean_diagonal
    }

    /// BLOSUM62 (Henikoff & Henikoff 1992), the default protein matrix.
    pub fn blosum62() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(|| protein("BLOSUM62", &BLOSUM62))
    }

    /// BLOSUM45, for distantly related proteins.
    pub fn blosum45() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(|| protein("BLOSUM45", &BLOSUM45))
    }

    /// BLOSUM80, for closely related proteins.
    pub fn blosum80() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(|| protein("BLOSUM80", &BLOSUM80))
    }

    /// PAM250 (Dayhoff), log-odds at 250 accepted point mutations per 100
    /// residues.
    pub fn pam250() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(|| protein("PAM250", &PAM250))
    }

    /// ClustalW's IUB DNA matrix: every match scores 1.9, every mismatch 0,
    /// and `N`/`X` match any IUB ambiguity symbol.
    pub fn iub() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(build_iub)
    }

    /// ClustalW's transition-weighted DNA matrix (`swgapdnamt` in ClustalW's
    /// `matrices.h`), scaled to units of 1/100: match 0.91-1.00, transition
    /// -0.31, transversion -1.14 to -1.25.
    pub fn clustal_dna() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(build_clustal_dna)
    }

    /// Identity matrix: 1.0 for a match (case-insensitive, `T` == `U`), 0.0
    /// otherwise.
    pub fn identity() -> &'static SubstMatrix {
        static M: OnceLock<SubstMatrix> = OnceLock::new();
        M.get_or_init(build_identity)
    }

    /// Resolve a [`MatrixChoice`]. `Auto` picks the IUB matrix for nucleotide
    /// data and BLOSUM62 for protein.
    pub fn choose(choice: MatrixChoice, alphabet: Alphabet) -> &'static SubstMatrix {
        match choice {
            MatrixChoice::Auto => {
                if alphabet.is_nucleotide() {
                    SubstMatrix::iub()
                } else {
                    SubstMatrix::blosum62()
                }
            }
            MatrixChoice::Blosum62 => SubstMatrix::blosum62(),
            MatrixChoice::Blosum45 => SubstMatrix::blosum45(),
            MatrixChoice::Blosum80 => SubstMatrix::blosum80(),
            MatrixChoice::Pam250 => SubstMatrix::pam250(),
            MatrixChoice::Identity => SubstMatrix::identity(),
            MatrixChoice::Iub => SubstMatrix::iub(),
            MatrixChoice::ClustalDna => SubstMatrix::clustal_dna(),
        }
    }
}

/// Empty table filled with `fill`.
fn blank(fill: f32) -> Box<[f32; CELLS]> {
    let v = vec![fill; CELLS].into_boxed_slice();
    // The vector is created with exactly CELLS elements, so the conversion
    // cannot fail.
    v.try_into().unwrap_or_else(|_| unreachable!("vec![_; CELLS] has CELLS elements"))
}

#[inline]
fn put(table: &mut [f32; CELLS], a: u8, b: u8, v: f32) {
    let write = |t: &mut [f32; CELLS], x: u8, y: u8| {
        let i = (x as usize) & (SIZE - 1);
        let j = (y as usize) & (SIZE - 1);
        t[(i << 7) | j] = v;
    };
    let (al, au) = (a.to_ascii_lowercase(), a.to_ascii_uppercase());
    let (bl, bu) = (b.to_ascii_lowercase(), b.to_ascii_uppercase());
    for &x in &[au, al] {
        for &y in &[bu, bl] {
            write(table, x, y);
        }
    }
}

/// Expand one of the 24x24 NCBI protein tables.
fn protein(name: &'static str, data: &[i8; 24 * 24]) -> SubstMatrix {
    // Unknown bytes behave like `X` against `X`.
    let x = AA_ORDER.iter().position(|&c| c == b'X').unwrap_or(22);
    let mut table = blank(data[x * 24 + x] as f32);
    for (i, &ai) in AA_ORDER.iter().enumerate() {
        for (j, &bj) in AA_ORDER.iter().enumerate() {
            put(&mut table, ai, bj, data[i * 24 + j] as f32);
        }
    }
    // Residues outside the NCBI 24: treat `J` (Leu/Ile) and `?` as `X`,
    // selenocysteine `U` as `C` and pyrrolysine `O` as `K`, which is what most
    // aligners do rather than refusing the input.
    for (extra, like) in [(b'J', b'X'), (b'?', b'X'), (b'U', b'C'), (b'O', b'K')] {
        for (i, &ai) in AA_ORDER.iter().enumerate() {
            let v = data[i * 24 + AA_ORDER.iter().position(|&c| c == like).unwrap_or(x)] as f32;
            put(&mut table, ai, extra, v);
            put(&mut table, extra, ai, v);
        }
        for (other, like2) in [(b'J', b'X'), (b'?', b'X'), (b'U', b'C'), (b'O', b'K')] {
            let vi = AA_ORDER.iter().position(|&c| c == like).unwrap_or(x);
            let vj = AA_ORDER.iter().position(|&c| c == like2).unwrap_or(x);
            put(&mut table, extra, other, data[vi * 24 + vj] as f32);
        }
    }
    // Gaps never reach `score` in the DP (they are handled by the gap model),
    // but keep them at zero so a stray gap cannot skew a score.
    for c in 0u8..128 {
        put(&mut table, b'-', c, 0.0);
        put(&mut table, c, b'-', 0.0);
    }
    let mean = mean_diagonal(&table, b"ACDEFGHIKLMNPQRSTVWY");
    SubstMatrix { name, table, mean_diagonal: mean }
}

fn mean_diagonal(table: &[f32; CELLS], symbols: &[u8]) -> f32 {
    if symbols.is_empty() {
        return 1.0;
    }
    let mut sum = 0.0;
    for &c in symbols {
        let i = (c as usize) & (SIZE - 1);
        sum += table[(i << 7) | i];
    }
    sum / symbols.len() as f32
}

/// IUPAC nucleotide code -> the set of bases it stands for, as a 4-bit mask
/// over A, C, G, T/U.
fn iupac_mask(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101, // A/G
        b'Y' => 0b1010, // C/T
        b'S' => 0b0110, // C/G
        b'W' => 0b1001, // A/T
        b'K' => 0b1100, // G/T
        b'M' => 0b0011, // A/C
        b'B' => 0b1110, // C/G/T
        b'D' => 0b1101, // A/G/T
        b'H' => 0b1011, // A/C/T
        b'V' => 0b0111, // A/C/G
        b'N' | b'X' | b'?' => 0b1111,
        _ => 0,
    }
}

/// ClustalW's IUB matrix. From the ClustalW documentation: "all matches score
/// 1.9; all mismatches for IUB symbols score 0", with `N`/`X` treated as a
/// match to any symbol.
fn build_iub() -> SubstMatrix {
    let mut table = blank(0.0);
    let mut symbols: Vec<u8> = NT_ORDER.to_vec();
    symbols.push(b'X');
    symbols.push(b'?');
    for &a in &symbols {
        for &b in &symbols {
            let (ma, mb) = (iupac_mask(a), iupac_mask(b));
            if ma == 0 || mb == 0 {
                continue;
            }
            let ambiguous = ma == 0b1111 || mb == 0b1111;
            let same = a.eq_ignore_ascii_case(&b)
                || (matches!(a.to_ascii_uppercase(), b'T' | b'U')
                    && matches!(b.to_ascii_uppercase(), b'T' | b'U'));
            let v = if same || ambiguous { 1.9 } else { 0.0 };
            put(&mut table, a, b, v);
        }
    }
    for c in 0u8..128 {
        put(&mut table, b'-', c, 0.0);
        put(&mut table, c, b'-', 0.0);
    }
    let mean = mean_diagonal(&table, b"ACGT");
    SubstMatrix { name: "IUB", table, mean_diagonal: mean }
}

/// ClustalW's `swgapdnamt` DNA matrix (ClustalW `matrices.h`), lower triangle
/// over A, C, G, T, divided by 100.
fn build_clustal_dna() -> SubstMatrix {
    // A    C     G     T
    const D: [[f32; 4]; 4] = [
        [0.91, -1.14, -0.31, -1.23],
        [-1.14, 1.00, -1.25, -0.31],
        [-0.31, -1.25, 1.00, -1.14],
        [-1.23, -0.31, -1.14, 0.91],
    ];
    let bases = b"ACGT";
    let mut table = blank(0.0);
    let symbols: Vec<u8> = NT_ORDER.iter().copied().chain(*b"X?").collect();
    for &a in &symbols {
        for &b in &symbols {
            let (ma, mb) = (iupac_mask(a), iupac_mask(b));
            if ma == 0 || mb == 0 {
                continue;
            }
            // Ambiguity codes score as the mean over the bases they cover,
            // which reduces to the plain entry for unambiguous symbols.
            let mut sum = 0.0;
            let mut n = 0.0;
            for (i, _) in bases.iter().enumerate() {
                if ma & (1 << i) == 0 {
                    continue;
                }
                for (j, _) in bases.iter().enumerate() {
                    if mb & (1 << j) == 0 {
                        continue;
                    }
                    sum += D[i][j];
                    n += 1.0;
                }
            }
            put(&mut table, a, b, if n > 0.0 { sum / n } else { 0.0 });
        }
    }
    for c in 0u8..128 {
        put(&mut table, b'-', c, 0.0);
        put(&mut table, c, b'-', 0.0);
    }
    let mean = mean_diagonal(&table, b"ACGT");
    SubstMatrix { name: "ClustalW DNA", table, mean_diagonal: mean }
}

fn build_identity() -> SubstMatrix {
    let mut table = blank(0.0);
    for c in b'A'..=b'Z' {
        put(&mut table, c, c, 1.0);
    }
    for c in b'0'..=b'9' {
        put(&mut table, c, c, 1.0);
    }
    // T and U are the same residue in an identity matrix.
    put(&mut table, b'T', b'U', 1.0);
    put(&mut table, b'U', b'T', 1.0);
    for c in 0u8..128 {
        put(&mut table, b'-', c, 0.0);
        put(&mut table, c, b'-', 0.0);
    }
    SubstMatrix { name: "identity", table, mean_diagonal: 1.0 }
}

// ---------------------------------------------------------------------------
// NCBI tables, row/column order A R N D C Q E G H I L K M F P S T W Y V B Z X *
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const BLOSUM62: [i8; 24 * 24] = [
     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,
    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,
    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1,
];

#[rustfmt::skip]
const BLOSUM45: [i8; 24 * 24] = [
     5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -1, -1,  0, -5,
    -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -1,  0, -1, -5,
    -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3,  4,  0, -1, -5,
    -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3,  5,  1, -1, -5,
    -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -3, -2, -5,
    -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3,  0,  4, -1, -5,
    -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3,  1,  4, -1, -5,
     0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -5,
    -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3,  0,  0, -1, -5,
    -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -3, -3, -1, -5,
    -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -3, -2, -1, -5,
    -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2,  0,  1, -1, -5,
    -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -2, -1, -1, -5,
    -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -3, -3, -1, -5,
    -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -2, -1, -1, -5,
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1,  0,  0,  0, -5,
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0,  0, -1,  0, -5,
    -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -4, -2, -2, -5,
    -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -2, -2, -1, -5,
     0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -3, -3, -1, -5,
    -1, -1,  4,  5, -2,  0,  1, -1,  0, -3, -3,  0, -2, -3, -2,  0,  0, -4, -2, -3,  4,  2, -1, -5,
    -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -2, -3,  2,  4, -1, -5,
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -2, -1, -1, -1, -1, -1, -5,
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1,
];

#[rustfmt::skip]
const BLOSUM80: [i8; 24 * 24] = [
     5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -6,
    -2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -1,  0, -1, -6,
    -2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4,  5,  0, -1, -6,
    -2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4,  5,  1, -1, -6,
    -1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -4, -1, -6,
    -1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3,  0,  4, -1, -6,
    -1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3,  1,  5, -1, -6,
     0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -3, -1, -6,
    -2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -1,  0, -1, -6,
    -2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -4, -4, -1, -6,
    -2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -4, -3, -1, -6,
    -1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -1,  1, -1, -6,
    -1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -3, -2, -1, -6,
    -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -4, -4, -1, -6,
    -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -2, -2, -1, -6,
     1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2,  0,  0, -1, -6,
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -1, -1, -1, -6,
    -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5, -3, -1, -6,
    -2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -3, -3, -1, -6,
     0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -4, -3, -1, -6,
    -2, -1,  5,  5, -4,  0,  1, -1, -1, -4, -4, -1, -3, -4, -2,  0, -1, -5, -3, -4,  5,  0, -1, -6,
    -1,  0,  0,  1, -4,  4,  5, -3,  0, -4, -3,  1, -2, -4, -2,  0, -1, -3, -3, -3,  0,  5, -1, -6,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6,
    -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1,
];

#[rustfmt::skip]
const PAM250: [i8; 24 * 24] = [
     2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0, -8,
    -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1, -8,
     0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0, -8,
     0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8,
    -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3, -8,
     0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1, -8,
     0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8,
     1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1, -8,
    -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1, -8,
    -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1, -8,
    -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1, -8,
    -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1, -8,
    -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1, -8,
    -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2, -8,
     1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1, -8,
     1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0, -8,
     1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0, -8,
    -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4, -8,
    -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2, -8,
     0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1, -8,
     0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1, -8,
     0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1, -8,
     0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1, -8,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1,
];

#[cfg(test)]
mod tests {
    use super::*;

    fn all() -> Vec<&'static SubstMatrix> {
        vec![
            SubstMatrix::blosum62(),
            SubstMatrix::blosum45(),
            SubstMatrix::blosum80(),
            SubstMatrix::pam250(),
            SubstMatrix::iub(),
            SubstMatrix::clustal_dna(),
            SubstMatrix::identity(),
        ]
    }

    #[test]
    fn blosum62_spot_checks() {
        let m = SubstMatrix::blosum62();
        assert_eq!(m.score(b'W', b'W'), 11.0);
        assert_eq!(m.score(b'C', b'C'), 9.0);
        assert_eq!(m.score(b'A', b'A'), 4.0);
        assert_eq!(m.score(b'H', b'H'), 8.0);
        assert_eq!(m.score(b'P', b'W'), -4.0);
        assert_eq!(m.score(b'F', b'Y'), 3.0);
        // B (Asx), Z (Glx) and X rows.
        assert_eq!(m.score(b'B', b'D'), 4.0);
        assert_eq!(m.score(b'B', b'N'), 3.0);
        assert_eq!(m.score(b'B', b'B'), 4.0);
        assert_eq!(m.score(b'Z', b'E'), 4.0);
        assert_eq!(m.score(b'Z', b'Q'), 3.0);
        assert_eq!(m.score(b'X', b'X'), -1.0);
        assert_eq!(m.score(b'X', b'A'), 0.0);
        assert_eq!(m.score(b'X', b'C'), -2.0);
        assert_eq!(m.score(b'*', b'*'), 1.0);
        assert_eq!(m.score(b'*', b'A'), -4.0);
    }

    #[test]
    fn blosum45_and_80_spot_checks() {
        let b45 = SubstMatrix::blosum45();
        assert_eq!(b45.score(b'C', b'C'), 12.0);
        assert_eq!(b45.score(b'W', b'W'), 15.0);
        assert_eq!(b45.score(b'R', b'K'), 3.0);
        assert_eq!(b45.score(b'B', b'D'), 5.0);
        let b80 = SubstMatrix::blosum80();
        assert_eq!(b80.score(b'C', b'C'), 9.0);
        assert_eq!(b80.score(b'W', b'W'), 11.0);
        assert_eq!(b80.score(b'A', b'A'), 5.0);
        assert_eq!(b80.score(b'Z', b'E'), 5.0);
    }

    #[test]
    fn pam250_spot_checks() {
        let m = SubstMatrix::pam250();
        assert_eq!(m.score(b'W', b'W'), 17.0);
        assert_eq!(m.score(b'C', b'C'), 12.0);
        assert_eq!(m.score(b'Y', b'F'), 7.0);
        assert_eq!(m.score(b'A', b'A'), 2.0);
        assert_eq!(m.score(b'B', b'D'), 3.0);
    }

    /// A transcription slip in one of the big tables almost always breaks
    /// symmetry, so check every matrix is symmetric over the printable range.
    #[test]
    fn matrices_are_symmetric() {
        for m in all() {
            for a in b'A'..=b'Z' {
                for b in b'A'..=b'Z' {
                    assert!(
                        (m.score(a, b) - m.score(b, a)).abs() < 1e-6,
                        "{} is asymmetric at {}/{}: {} vs {}",
                        m.name(),
                        a as char,
                        b as char,
                        m.score(a, b),
                        m.score(b, a)
                    );
                }
            }
        }
    }

    #[test]
    fn lowercase_folds_to_uppercase() {
        for m in all() {
            for a in b'A'..=b'Z' {
                for b in b'A'..=b'Z' {
                    let up = m.score(a, b);
                    assert_eq!(m.score(a.to_ascii_lowercase(), b), up);
                    assert_eq!(m.score(a, b.to_ascii_lowercase()), up);
                    assert_eq!(m.score(a.to_ascii_lowercase(), b.to_ascii_lowercase()), up);
                }
            }
        }
    }

    #[test]
    fn iub_matches_and_mismatches() {
        let m = SubstMatrix::iub();
        assert_eq!(m.score(b'A', b'A'), 1.9);
        assert_eq!(m.score(b'A', b'G'), 0.0);
        assert_eq!(m.score(b'T', b'U'), 1.9);
        assert_eq!(m.score(b'N', b'A'), 1.9);
        assert_eq!(m.score(b'a', b'a'), 1.9);
        assert_eq!(m.score(b'-', b'A'), 0.0);
    }

    #[test]
    fn clustal_dna_prefers_transitions_over_transversions() {
        let m = SubstMatrix::clustal_dna();
        assert!(m.score(b'A', b'A') > 0.0);
        assert!(m.score(b'A', b'G') > m.score(b'A', b'C'));
        assert!(m.score(b'C', b'T') > m.score(b'C', b'G'));
    }

    #[test]
    fn identity_is_one_or_zero() {
        let m = SubstMatrix::identity();
        assert_eq!(m.score(b'K', b'K'), 1.0);
        assert_eq!(m.score(b'K', b'R'), 0.0);
        assert_eq!(m.score(b'T', b'u'), 1.0);
    }

    #[test]
    fn choose_auto_by_alphabet() {
        assert_eq!(SubstMatrix::choose(MatrixChoice::Auto, Alphabet::Dna).name(), "IUB");
        assert_eq!(SubstMatrix::choose(MatrixChoice::Auto, Alphabet::Rna).name(), "IUB");
        assert_eq!(SubstMatrix::choose(MatrixChoice::Auto, Alphabet::Protein).name(), "BLOSUM62");
        assert_eq!(SubstMatrix::choose(MatrixChoice::Blosum80, Alphabet::Dna).name(), "BLOSUM80");
    }

    #[test]
    fn unknown_bytes_do_not_panic() {
        let m = SubstMatrix::blosum62();
        for a in 0u8..=255 {
            let _ = m.score(a, b'A');
            let _ = m.score(b'A', a);
            let _ = m.score(a, a);
        }
    }

    #[test]
    fn mean_diagonal_is_sane() {
        assert!((SubstMatrix::blosum62().mean_diagonal() - 5.8).abs() < 0.05);
        assert!((SubstMatrix::iub().mean_diagonal() - 1.9).abs() < 1e-6);
    }
}
