# TOLViewer architecture and inter-crate API contract

TOLViewer is a Rust workspace. Every crate is written against `tolviewer-core`,
which has no dependencies of its own. **This document is the contract: the
public signatures below are fixed, because other crates are written against them
in parallel. Do not change them without changing this file.**

```
tolviewer-core   data model: Sequence, Alignment, Alphabet, EditOp/UndoStack, Consensus/ColumnStats, Error
      |
      +-- tolviewer-io      readers/writers: FASTA, FASTQ, PHYLIP, NEXUS, Clustal, Stockholm, GenBank, MSF
      +-- tolviewer-align   MSA engines: Clustal-style, MUSCLE-style, MAFFT-style; matrices, trees, pairwise
      +-- tolviewer-clean   Gblocks column selection
      |
      +-- tolviewer-app     eframe/egui GUI, binary `tolviewer`
```

## Ground rules for every crate

* No `unsafe`. No panics on user input — return `tolviewer_core::Error`.
  `unwrap()`/`expect()` only where a comment justifies the invariant.
* `#![forbid(unsafe_code)]` at the top of every `lib.rs`.
* Residues are ASCII `u8`, **case preserved** (lowercase is used by some tools
  to mark masked/low-quality regions; never uppercase the user's data in place).
  Compare case-insensitively.
* The gap character is `tolviewer_core::GAP` (`b'-'`). `.` and `~` are
  normalised to it on input by `Sequence::new`.
* Long-running work must be interruptible and must report progress: take
  `&dyn Progress` (defined in `tolviewer-align`, re-exported where needed).
* Every public item gets a doc comment. Tests go in `#[cfg(test)] mod tests`
  in the same file, plus `tests/` for round-trip / integration tests.
* Prefer plain `std` + `rayon`. Do not add dependencies beyond those already in
  the crate's `Cargo.toml` without a strong reason; if you must, keep them
  well-maintained and permissively licensed, and note why in the PR/summary.

## `tolviewer-core` (already implemented — read the source)

Key items you will use:

```rust
pub enum Alphabet { Dna, Rna, Protein }
impl Alphabet {
    fn is_nucleotide(self) -> bool;
    fn name(self) -> &'static str;
    fn nexus_datatype(self) -> &'static str;
    fn symbols(self) -> &'static [u8];       // incl. IUPAC ambiguity
    fn core_symbols(self) -> &'static [u8];  // ACGT / ACGU / 20 aa
    fn is_valid(self, c: u8) -> bool;
    fn is_ambiguous(self, c: u8) -> bool;
    fn guess(residues: impl IntoIterator<Item = u8>) -> Alphabet;
    fn complement(self, c: u8) -> u8;
}
pub const GAP: u8;  pub const MISSING: u8;  pub fn is_gap(c: u8) -> bool;

pub struct Sequence {
    pub id: String, pub description: String,
    pub residues: Vec<u8>,           // gaps included
    pub quality: Option<Vec<u8>>,    // Phred, same length as residues, 0 at gaps
    pub hidden: bool,
}
impl Sequence {
    fn new(id: impl Into<String>, residues: impl Into<Vec<u8>>) -> Self;
    fn len(&self) -> usize;  fn ungapped_len(&self) -> usize;  fn ungapped(&self) -> Vec<u8>;
    fn header(&self) -> String;  fn set_header(&mut self, header: &str);
    fn residue_index_at(&self, column: usize) -> Option<usize>;
    fn reverse_complement(&mut self, alphabet: Alphabet);
    fn pad_to(&mut self, width: usize);
    fn mean_quality(&self) -> Option<f32>;
    fn ambiguity_fraction(&self, alphabet: Alphabet) -> f32;
}

pub struct Alignment { pub name: String, pub sequences: Vec<Sequence>, /* private alphabet cache */ }
impl Alignment {
    fn new(name: impl Into<String>, sequences: Vec<Sequence>) -> Self;
    fn len(&self) -> usize;              // rows
    fn width(&self) -> usize;            // columns (longest row)
    fn is_aligned(&self) -> bool;        // all rows equal length
    fn alphabet(&mut self) -> Alphabet;  // guesses + caches
    fn alphabet_hint(&self) -> Option<Alphabet>;
    fn set_alphabet(&mut self, a: Alphabet);
    fn pad_to_width(&mut self);
    fn require_aligned(&self) -> Result<()>;
    fn get(&self, row: usize, col: usize) -> Option<u8>;
    fn set(&mut self, row: usize, col: usize, residue: u8) -> Result<u8>;
    fn insert_columns(&mut self, at: usize, count: usize) -> Result<()>;
    fn delete_columns(&mut self, start: usize, end: usize) -> Result<Vec<Vec<u8>>>;
    fn restore_columns(&mut self, at: usize, columns: &[Vec<u8>]) -> Result<()>;
    fn keep_columns(&mut self, mask: &[bool]) -> Result<usize>;
    fn insert_gap(&mut self, row: usize, col: usize) -> Result<()>;
    fn delete_at(&mut self, row: usize, col: usize) -> Result<u8>;
    fn all_gap_columns(&self) -> Vec<usize>;
    fn remove_all_gap_columns(&mut self) -> usize;
    fn degap(&mut self);
    fn remove_sequence(&mut self, row: usize) -> Result<Sequence>;
    fn insert_sequence(&mut self, row: usize, seq: Sequence) -> Result<()>;
    fn move_sequence(&mut self, from: usize, to: usize) -> Result<()>;
    fn subset(&self, rows: &[usize], cols: Range<usize>) -> Alignment;
    fn find_by_id(&self, id: &str) -> Option<usize>;
    fn deduplicate_ids(&mut self) -> usize;
    fn column(&self, col: usize) -> impl Iterator<Item = u8> + '_;
}

pub enum Error { Io(..), Parse{format,line,message}, Format(String), NotAligned,
                 OutOfRange(String), Algorithm(String), Cancelled }
pub type Result<T> = std::result::Result<T, Error>;
// constructors: Error::parse(fmt, line, msg), Error::format(msg),
//               Error::algorithm(msg), Error::out_of_range(msg)

pub struct ColumnStats { pub occupancy: u32, pub rows: u32, pub majority: Option<u8>,
                         pub majority_count: u32, pub distinct: u8, pub ambiguous: u32 }
pub struct Consensus { pub residues: Vec<u8>, pub columns: Vec<ColumnStats> }
impl Consensus { fn compute(&Alignment, Alphabet, threshold: f32, min_occupancy: f32) -> Consensus; }
pub fn pairwise_identity(a: &[u8], b: &[u8]) -> Option<f32>;
```

## `tolviewer-io` — required public API

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Format { Fasta, Fastq, Phylip, PhylipRelaxed, Nexus, Clustal, Stockholm, Msf, Genbank }

impl Format {
    /// Human name for menus, e.g. "FASTA".
    pub fn name(self) -> &'static str;
    /// Lowercase extensions without the dot, first is the default for saving.
    pub fn extensions(self) -> &'static [&'static str];
    pub fn can_read(self) -> bool;
    pub fn can_write(self) -> bool;
    /// All formats, for menu construction.
    pub fn all() -> &'static [Format];
    /// Guess from a file extension.
    pub fn from_path(path: &Path) -> Option<Format>;
    /// Guess from the first bytes of the file. Prefer this over the extension
    /// when they disagree and the content is unambiguous.
    pub fn sniff(bytes: &[u8]) -> Option<Format>;
}

/// Read a file, detecting the format from content then extension.
pub fn read_file(path: &Path) -> Result<Alignment>;
/// Read a file with an explicitly chosen format.
pub fn read_file_as(path: &Path, format: Format) -> Result<Alignment>;
/// Parse in-memory bytes. `name` becomes `Alignment::name`.
pub fn parse(bytes: &[u8], format: Format, name: &str) -> Result<Alignment>;

#[derive(Debug, Clone)]
pub struct WriteOptions {
    /// Residues per line for sequential formats (FASTA); 0 = one line per seq.
    pub line_width: usize,
    /// Write interleaved blocks (PHYLIP/NEXUS/Clustal).
    pub interleaved: bool,
    /// Residues per interleaved block.
    pub block_width: usize,
    /// Uppercase all residues on output.
    pub uppercase: bool,
    /// Truncate/pad names to 10 chars for strict PHYLIP.
    pub strict_phylip_names: bool,
    /// Include rows flagged `hidden`.
    pub include_hidden: bool,
    /// Replace characters illegal in the target format (whitespace, quotes,
    /// parentheses, colons, semicolons) in names with `_`.
    pub sanitize_names: bool,
    /// Line ending to emit.
    pub line_ending: LineEnding,
}
impl Default for WriteOptions;  // line_width 60, interleaved false, block_width 60,
                                // uppercase false, strict_phylip_names false,
                                // include_hidden false, sanitize_names true, Lf
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineEnding { Lf, Crlf }

pub fn write_file(aln: &Alignment, path: &Path, format: Format, opts: &WriteOptions) -> Result<()>;
pub fn write_string(aln: &Alignment, format: Format, opts: &WriteOptions) -> Result<String>;
```

Writing must fail with `Error::Format` (not panic, not silently truncate) when
the data cannot be represented — e.g. PHYLIP/NEXUS given a ragged (unaligned)
set. `write_file` writes atomically: write a temp file in the destination
directory, then rename over the target.

## `tolviewer-align` — required public API

```rust
/// Progress + cancellation for long operations. Return `false` from `tick` to
/// request cancellation; algorithms must then return `Err(Error::Cancelled)`.
pub trait Progress: Sync {
    /// `fraction` is 0.0..=1.0. Called from worker threads; implementations
    /// must be cheap and thread-safe.
    fn tick(&self, fraction: f32, message: &str) -> bool;
}
/// A `Progress` that never cancels and discards messages.
pub struct NoProgress;
impl Progress for NoProgress { .. }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Engine {
    /// Progressive alignment with a distance-based guide tree, in the style of
    /// ClustalW/Clustal Omega: pairwise distances -> NJ guide tree ->
    /// profile-profile progressive alignment with position-specific gap penalties.
    Clustal,
    /// Progressive draft then iterative refinement by tree-dependent
    /// restricted partitioning, in the style of MUSCLE: k-mer distance draft ->
    /// re-estimate tree from the draft -> re-align -> horizontal refinement.
    Muscle,
    /// FFT-accelerated group-to-group alignment in the style of MAFFT FFT-NS-2:
    /// residues mapped to volume/polarity vectors, homologous segments located
    /// by FFT correlation, progressive alignment, then a second pass on a tree
    /// re-estimated from the first alignment.
    Mafft,
}
impl Engine { pub fn name(self) -> &'static str; pub fn all() -> &'static [Engine]; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixChoice { Auto, Blosum62, Blosum45, Blosum80, Pam250, Identity, Iub, ClustalDna }

#[derive(Debug, Clone)]
pub struct AlignParams {
    pub engine: Engine,
    pub matrix: MatrixChoice,
    /// Affine gap penalties, positive numbers (they are subtracted).
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
impl Default for AlignParams;
impl AlignParams { pub fn for_engine(engine: Engine) -> Self; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeMethod { NeighborJoining, Upgma }

/// Align the (possibly already gapped) sequences. Existing gaps are stripped
/// first. Row order and names are preserved; the result is always aligned.
pub fn align(aln: &Alignment, params: &AlignParams, progress: &dyn Progress) -> Result<Alignment>;

/// Re-align only columns `cols` of an existing alignment, keeping the rest
/// fixed. Used by "realign selection" in the GUI.
pub fn realign_region(aln: &Alignment, cols: Range<usize>, params: &AlignParams,
                      progress: &dyn Progress) -> Result<Alignment>;

/// Align `query` sequences onto an existing fixed alignment (profile-sequence).
pub fn add_to_alignment(profile: &Alignment, query: &[Sequence], params: &AlignParams,
                        progress: &dyn Progress) -> Result<Alignment>;

pub mod matrix {
    pub struct SubstMatrix { /* 256x256 lookup over ASCII */ }
    impl SubstMatrix {
        pub fn score(&self, a: u8, b: u8) -> f32;
        pub fn name(&self) -> &str;
        pub fn blosum62() -> &'static SubstMatrix;
        pub fn blosum45() -> &'static SubstMatrix;
        pub fn blosum80() -> &'static SubstMatrix;
        pub fn pam250() -> &'static SubstMatrix;
        /// ClustalW's IUB DNA matrix (match +1.9, mismatch 0).
        pub fn iub() -> &'static SubstMatrix;
        pub fn identity() -> &'static SubstMatrix;
        pub fn choose(choice: MatrixChoice, alphabet: Alphabet) -> &'static SubstMatrix;
    }
}

pub mod pairwise {
    /// Global alignment with affine gaps (Gotoh). Returns the two gapped rows
    /// and the score.
    pub fn global(a: &[u8], b: &[u8], m: &SubstMatrix, gap_open: f32, gap_extend: f32)
        -> (Vec<u8>, Vec<u8>, f32);
    /// Local alignment (Smith-Waterman) with affine gaps.
    pub fn local(a: &[u8], b: &[u8], m: &SubstMatrix, gap_open: f32, gap_extend: f32)
        -> (Vec<u8>, Vec<u8>, f32, Range<usize>, Range<usize>);
}

pub mod distance {
    /// Fraction of differing residues over aligned, non-gap positions.
    pub fn p_distance(a: &[u8], b: &[u8]) -> f32;
    /// Jukes-Cantor corrected distance; saturates at `max` for p >= 0.75.
    pub fn jukes_cantor(p: f32) -> f32;
    /// Kimura's protein correction, as used by ClustalW.
    pub fn kimura_protein(p: f32) -> f32;
    /// Fast alignment-free distance from shared k-mer counts (MUSCLE stage 1).
    pub fn kmer_distance(a: &[u8], b: &[u8], k: usize, alphabet: Alphabet) -> f32;
    /// Full pairwise distance matrix, lower triangle, parallel over pairs.
    pub fn matrix(seqs: &[Vec<u8>], method: DistanceMethod, alphabet: Alphabet,
                  progress: &dyn Progress) -> Result<DistMatrix>;
    pub enum DistanceMethod { Kmer { k: usize }, PairwiseAlignment, FromAlignment }
    pub struct DistMatrix { n: usize, data: Vec<f32> }
    impl DistMatrix { pub fn get(&self, i: usize, j: usize) -> f32; pub fn len(&self) -> usize; }
}

pub mod tree {
    /// Rooted guide tree over leaf indices.
    pub enum GuideTree { Leaf(usize), Node { left: Box<GuideTree>, right: Box<GuideTree>,
                                             left_len: f32, right_len: f32 } }
    impl GuideTree {
        pub fn leaves(&self) -> Vec<usize>;
        pub fn to_newick(&self, names: &[String]) -> String;
    }
    pub fn neighbor_joining(d: &DistMatrix) -> Result<GuideTree>;
    pub fn upgma(d: &DistMatrix) -> Result<GuideTree>;
}
```

Accuracy targets (checked by tests in `crates/tolviewer-align/tests/`): on
BAliBASE-style toy cases and on simulated sequences with known ancestry, each
engine must recover the reference alignment columns at high rates, and
`Engine::Muscle`/`Engine::Mafft` with `iterations >= 2` must be no worse than
`Engine::Clustal` on average. Include a test that 200 sequences x 1000 columns
aligns in reasonable time (mark `#[ignore]` if slow in debug).

## `tolviewer-clean` — required public API

Implements the Gblocks algorithm (Castresana 2000, *Mol Biol Evol* 17:540-552;
Talavera & Castresana 2007 for the relaxed settings), reimplemented from the
published description.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapPolicy { /// no gaps allowed in any kept column
                     None,
                     /// gaps allowed in up to half the sequences
                     Half,
                     /// gaps allowed in any number of sequences
                     All }

#[derive(Debug, Clone)]
pub struct GblocksParams {
    /// b1: minimum number of sequences for a conserved position. Gblocks
    /// requires a strict majority; the default is `n / 2 + 1` with integer
    /// division, which is what the original program uses ("50% + 1").
    pub min_seqs_conserved: usize,
    /// b2: minimum number of sequences for a flank position. Must be >= b1;
    /// default is ceil(n * 0.85).
    pub min_seqs_flank: usize,
    /// b3: maximum number of contiguous non-conserved positions. Default 8.
    pub max_contiguous_nonconserved: usize,
    /// b4: minimum length of a block. Default 10.
    pub min_block_length: usize,
    /// b5: allowed gap positions. Default `GapPolicy::None`.
    pub gaps: GapPolicy,
    /// Treat similar residues (positive substitution score) as conserved for
    /// the flank test, as Gblocks does for protein.
    pub use_similarity: bool,
}
impl GblocksParams {
    /// Gblocks' own defaults for `n` sequences.
    pub fn defaults(n_seqs: usize) -> Self;
    /// The relaxed settings of Talavera & Castresana (2007): b2 = b1,
    /// b3 = 8, b4 = 5, b5 = Half.
    pub fn relaxed(n_seqs: usize) -> Self;
    /// Reject impossible combinations (b1 <= n/2, b2 < b1, b4 < 2, ...).
    pub fn validate(&self, n_seqs: usize) -> Result<()>;
}

/// Per-column classification, for the GUI's cleaning track.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnFlag { Conserved, HighlyConserved, NonConserved, GapRich }

#[derive(Debug, Clone)]
pub struct GblocksResult {
    /// One entry per alignment column: true = keep.
    pub mask: Vec<bool>,
    /// Contiguous kept ranges.
    pub blocks: Vec<Range<usize>>,
    pub flags: Vec<ColumnFlag>,
    pub kept: usize,
    pub total: usize,
}
impl GblocksResult {
    pub fn kept_fraction(&self) -> f32;
    /// A new alignment with only the kept columns.
    pub fn apply(&self, aln: &Alignment) -> Result<Alignment>;
    /// The Gblocks-style mask line ("  ####  ...") for display/export.
    pub fn mask_line(&self) -> String;
}

pub fn gblocks(aln: &Alignment, params: &GblocksParams) -> Result<GblocksResult>;

/// Simple complementary filters the GUI also offers.
pub fn remove_gappy_columns(aln: &Alignment, max_gap_fraction: f32) -> Vec<bool>;
pub fn remove_gappy_sequences(aln: &Alignment, max_gap_fraction: f32) -> Vec<bool>;
/// Trim ragged 5'/3' ends to the first/last column with `min_occupancy` coverage.
pub fn trim_ends(aln: &Alignment, min_occupancy: f32) -> Range<usize>;
```

## `tolviewer-app`

eframe/egui. Owns the document model (open alignments, selection, undo stacks),
a custom-painted alignment canvas that only draws visible cells, and background
threads for align/clean that report through `Progress`.
