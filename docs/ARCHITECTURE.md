# TOLViewer architecture and inter-crate API contract

TOLViewer is a Rust workspace. Every crate is written against `tolviewer-core`,
which has no dependencies of its own. **This document is the contract: the
public signatures below are fixed, because other crates are written against them
in parallel. Do not change them without changing this file.**

```
tolviewer-core   data model: Sequence, Alignment, Alphabet, EditOp/UndoStack, Consensus/ColumnStats, Error
      |
      +-- tolviewer-io      readers/writers: FASTA, FASTQ, PHYLIP, NEXUS, Clustal, Stockholm, GenBank, MSF, AB1
      +-- tolviewer-align   MSA engines: Clustal-style, MUSCLE-style, MAFFT-style; matrices, trees, pairwise
      +-- tolviewer-clean   Gblocks column selection
      |     |
      |     +-- tolviewer-library   project library: folder tree over files, primers, concatenation
      |                             (depends on core + io)
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
pub enum Format { Fasta, Fastq, Phylip, PhylipRelaxed, Nexus, Clustal, Stockholm, Msf, Genbank,
                  Ab1 }

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

### AB1 traces

`Format::Ab1` is read-only and yields a one-row alignment of the base calls.
The chromatogram behind them is a separate entry point, because most callers
want the sequence and only the trace viewer wants the signal.

```rust
/// The format `read_file` would use, without parsing the whole file.
pub fn sniff_file(path: &Path) -> Result<Format>;

pub mod ab1 {
    /// A Sanger chromatogram. `channels` are all the same length; `calls`,
    /// `peaks` and `quality` are all the same length as each other.
    pub struct Trace {
        pub sample_name: String,
        /// The base each channel carries, in the file's channel order (`FWO_`).
        pub channel_bases: [u8; 4],
        pub channels: [Vec<u16>; 4],
        pub calls: Vec<u8>,
        /// The sample index in `channels` each call peaks at.
        pub peaks: Vec<u32>,
        pub quality: Option<Vec<u8>>,
        pub comment: String,
    }
    impl Trace {
        pub fn len(&self) -> usize;          // calls
        pub fn samples(&self) -> usize;      // trace length
        pub fn channel_for(&self, base: u8) -> Option<usize>;
        pub fn signal(&self, base: u8, sample: usize) -> u16;
        pub fn peak_signal(&self) -> u16;    // for scaling, never 0
        pub fn mean_peak_spacing(&self) -> f32;  // samples per call, for zoom decisions
        pub fn to_sequence(&self, id: &str) -> Sequence;
        pub fn reverse_complement(&mut self);
        pub fn set_call(&mut self, index: usize, base: u8) -> Result<()>;
        pub fn insert_call(&mut self, index: usize, base: u8) -> Result<()>;
        pub fn remove_call(&mut self, index: usize) -> Result<u8>;
        pub fn trim(&mut self, calls: Range<usize>) -> Result<()>;
        pub fn quality_trim_range(&self, window: usize, min_mean: f32) -> Range<usize>;
    }
    pub fn parse(bytes: &[u8], name: &str) -> Result<Trace>;
    pub fn read_file(path: &Path) -> Result<Trace>;
}
```

Files in the wild disagree about tag lengths; the reader trims to the shortest
consistent set rather than refusing a read. A missing `DATA` 9–12, `PBAS` or
`PLOC` is a real error, because there is nothing to show without them.

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

## `tolviewer-library` — required public API

A project library: a tree of folders over sequence files **that stay where they
are**. Adding a file reads it once and remembers where it is; the library never
moves, copies or rewrites the lab's data on its own. That rule is the reason
this crate exists and everything below serves it.

```rust
/// Stable handle into the tree. Ids are never reused, so a stale one from the
/// GUI resolves to `None` rather than to whatever took its place.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(u32);
impl NodeId { pub fn raw(self) -> u32; }

pub enum EntryKind { Sequences, Alignment, Trace }
impl EntryKind {
    pub fn name(self) -> &'static str;
    pub fn of(format: Format, alignment: &Alignment) -> EntryKind;
}

pub struct Entry {
    /// The lab's file. Never written without a confirmed `SaveChoice::Overwrite`.
    pub origin: PathBuf,
    pub format: Format,
    /// Set the first time edits are diverted to a copy; later saves go here.
    pub working: Option<PathBuf>,
    /// Which sequences of the file this entry is, by id. `None` = all of them.
    pub select: Option<Vec<String>>,
    /// Show reverse complemented. A flag, not a rewrite.
    pub reversed: bool,
    pub kind: EntryKind,
    pub note: String,
}
impl Entry {
    pub fn source_path(&self) -> &Path;      // the working copy once there is one
    pub fn load(&self) -> Result<Alignment>; // applies `select` and `reversed`
    pub fn load_trace(&self) -> Result<ab1::Trace>;
    pub fn save_target(&self) -> SaveTarget;
    pub fn suggested_copy(&self) -> PathBuf; // `foo.edited.fasta`, made unique
}

pub struct Folder { pub children: Vec<NodeId>, pub expanded: bool, pub note: String }
pub enum NodeKind { Folder(Folder), Entry(Box<Entry>) }
pub struct Node { pub name: String, pub parent: Option<NodeId>, pub kind: NodeKind }
impl Node {
    pub fn is_folder(&self) -> bool;
    pub fn entry(&self) -> Option<&Entry>;      pub fn entry_mut(&mut self) -> Option<&mut Entry>;
    pub fn folder(&self) -> Option<&Folder>;    pub fn folder_mut(&mut self) -> Option<&mut Folder>;
}

pub struct Library { pub name: String, pub primers: PrimerSet, pub path: Option<PathBuf>, /* arena */ }
impl Library {
    pub fn new(name: impl Into<String>) -> Self;
    // reading
    pub fn roots(&self) -> &[NodeId];
    pub fn get(&self, id: NodeId) -> Option<&Node>;
    pub fn get_mut(&mut self, id: NodeId) -> Option<&mut Node>;
    pub fn entry(&self, id: NodeId) -> Option<&Entry>;
    pub fn children(&self, parent: Option<NodeId>) -> &[NodeId];
    pub fn walk(&self) -> Vec<(NodeId, usize)>;           // depth first, with depth
    pub fn entries_under(&self, id: Option<NodeId>) -> Vec<NodeId>;
    pub fn path_of(&self, id: NodeId) -> String;          // "Project / 18S / TL-2213"
    pub fn is_ancestor(&self, ancestor: NodeId, id: NodeId) -> bool;
    pub fn broken_entries(&self) -> Vec<NodeId>;          // files that have gone missing
    pub fn revision(&self) -> u64;  pub fn is_dirty(&self) -> bool;
    pub fn mark_saved(&mut self);   pub fn touch(&mut self);
    // building
    pub fn add_folder(&mut self, parent: Option<NodeId>, name: impl Into<String>) -> Result<NodeId>;
    pub fn add_file(&mut self, parent: Option<NodeId>, path: &Path) -> Result<NodeId>;
    pub fn add_selection(&mut self, parent: Option<NodeId>, source: NodeId,
                         ids: Vec<String>, name: impl Into<String>) -> Result<NodeId>;
    pub fn remove(&mut self, id: NodeId) -> Vec<NodeId>;  // files are NOT deleted
    pub fn rename(&mut self, id: NodeId, name: impl Into<String>) -> Result<()>;
    pub fn move_to(&mut self, id: NodeId, parent: Option<NodeId>, index: usize) -> Result<()>;
    // using
    pub fn load(&self, id: NodeId) -> Result<Alignment>;
    /// Every entry under `ids`, names deduplicated, ready to align. Entries
    /// that fail to read are reported instead of losing the rest of the batch.
    pub fn gather(&self, ids: &[NodeId], name: &str) -> (Alignment, Vec<(NodeId, Error)>);
    pub fn set_reversed(&mut self, id: NodeId, reversed: bool) -> Result<()>;
    pub fn save_target(&self, id: NodeId) -> Result<SaveTarget>;
    pub fn save_entry(&mut self, id: NodeId, alignment: &Alignment,
                      choice: SaveChoice, options: &WriteOptions) -> Result<PathBuf>;
}
```

### The in-situ save policy

This is the contract that keeps a library safe to point at a shared drive.

```rust
pub enum CopyReason { PartOfAFile, ReadOnlyFormat }
impl CopyReason { pub fn explain(self) -> &'static str; }

pub enum SaveTarget {
    /// A copy the library already made. No confirmation needed.
    WorkingCopy(PathBuf),
    /// The lab's own file. The GUI must ask before this is written.
    Original(PathBuf),
    /// Only a copy is possible, at the suggested path.
    MustCopy(PathBuf, CopyReason),
}
impl SaveTarget {
    pub fn needs_confirmation(&self) -> bool;  // false only for WorkingCopy
    pub fn can_overwrite(&self) -> bool;       // false only for MustCopy
    pub fn path(&self) -> &Path;
}

pub enum SaveChoice { Overwrite, NewCopy(PathBuf) }
```

`SaveChoice::NewCopy` writes the copy **and remembers it**, so the question is
asked once per entry, not once per save, and editing a sequence ten times leaves
one extra file rather than ten. `SaveChoice::Overwrite` against a `MustCopy`
target is an error, not a silent copy.

Any successful save clears `select` and `reversed`: they describe how to turn
the file's contents into what the entry shows, and what was written *is* what
the entry shows. Leaving `reversed` set would reverse-complement the saved
sequence a second time the next time it was read.

### Primers

```rust
pub enum Strand { Forward, Reverse }
pub struct Primer { pub name: String, pub sequence: Vec<u8> }   // IUPAC, uppercase
impl Primer {
    pub fn new(name: impl Into<String>, sequence: &str) -> Result<Primer>;
    pub fn reverse_complement(&self) -> Vec<u8>;
    pub fn mismatch_budget(&self, fraction: f32) -> usize;      // never the whole primer
}
pub struct PrimerHit { pub primer: usize, pub name: String, pub strand: Strand,
                       pub range: Range<usize>, pub mismatches: usize }
impl PrimerHit { pub fn identity(&self) -> f32; }
pub struct PrimerSet { /* … */ }
impl PrimerSet {
    pub fn new(primers: Vec<Primer>) -> Self;
    pub fn primers(&self) -> &[Primer];
    pub fn push(&mut self, primer: Primer);
    pub fn remove(&mut self, index: usize) -> Option<Primer>;
    pub fn get(&self, index: usize) -> Option<&Primer>;
    pub fn find_by_name(&self, name: &str) -> Option<usize>;
    /// Every binding site on either strand, best first, overlaps collapsed.
    pub fn map(&self, seq: &[u8], max_mismatch_fraction: f32) -> Vec<PrimerHit>;
}

pub struct TrimOptions { pub max_mismatch_fraction: f32, pub search_window: usize,
                         pub keep_primers: bool }   // default 0.2, 120, false
pub struct TrimPlan { pub range: Range<usize>,
                      pub start_hit: Option<PrimerHit>, pub end_hit: Option<PrimerHit> }
impl TrimPlan { pub fn trims_anything(&self) -> bool; pub fn describe(&self, len: usize) -> String; }
pub fn plan_trim(set: &PrimerSet, seq: &[u8], opts: &TrimOptions) -> TrimPlan;
```

Matching compares IUPAC *code sets*, so a degenerate primer matches what it
stands for and an `N` in the read costs nothing — an uncalled base may well be
the base the primer wants, and the primer sits where the basecaller is least
sure. A gap never matches. The 3' hit must be a different binding site from the
5' one, so a read that is nothing but primer is left alone rather than emptied.

### Name matching and concatenation

```rust
pub struct MatchOptions { pub normalize: bool, pub strip_suffixes: bool,
                          pub extra_suffixes: Vec<String> }
impl MatchOptions { pub fn exact() -> Self; }
pub const DEFAULT_SUFFIXES: &[&str];
/// The sample a sequence id belongs to. Never empty for a non-empty name.
pub fn sample_key(name: &str, opts: &MatchOptions) -> String;
pub fn group<'a>(names: impl IntoIterator<Item = &'a str>, opts: &MatchOptions)
    -> Vec<(String, Vec<&'a str>)>;

pub struct ConcatOptions { pub matching: MatchOptions, pub include_partial: bool }
pub struct Partition { pub name: String, pub range: Range<usize> }
impl Partition { pub fn as_charset(&self) -> String; }     // 1-based, inclusive
pub struct MissingSample { pub sample: String, pub absent_from: Vec<String> }
pub struct ConcatResult {
    pub alignment: Alignment, pub partitions: Vec<Partition>,
    pub missing: Vec<MissingSample>, pub complete: usize, pub dropped: Vec<String>,
}
impl ConcatResult {
    pub fn nexus_charsets(&self) -> String;
    pub fn raxml_partitions(&self, model: &str) -> String;
}
pub fn concatenate(parts: &[&Alignment], opts: &ConcatOptions) -> Result<ConcatResult>;

pub struct SamplePreview { pub key: String, pub display: String,
                           pub found_in: Vec<(usize, String)> }
impl SamplePreview { pub fn is_complete(&self, total: usize) -> bool; }
/// What concatenating would do, without building the matrix. The GUI shows
/// this first, because a name-matching mistake is invisible in the result.
pub fn preview(parts: &[&Alignment], opts: &ConcatOptions) -> Vec<SamplePreview>;
```

Every input must already be aligned — concatenating ragged rows would silently
shift every locus after the first. Two rows of one locus that reduce to the same
sample is an `Error::Format` naming both, not a silent pick. Row order is the
first alignment's, then each later one's newcomers in its own order, so the same
files always give a byte-identical matrix.

### The library file

```rust
pub mod store {
    pub const EXTENSION: &str;                 // "tolvlib"
    pub fn write_string(library: &Library, at: &Path) -> String;
    pub fn parse(text: &str, at: &Path) -> Result<Library>;
    pub fn save(library: &mut Library, path: &Path) -> Result<()>;   // atomic
    pub fn load(path: &Path) -> Result<Library>;
}
```

Plain text, tab-indented to mirror the tree, versioned by its first line. Paths
under the library file's own directory are stored relative to it so a project
folder can be moved or shared; anything outside stays absolute. Unknown keys
inside a known block are **skipped, not rejected**, so a library written by a
later version still opens. Bump the version only when an old reader could
*misread* a new file — adding keys does not need it.

## `tolviewer-app`

eframe/egui. Owns the document model (open alignments, selection, undo stacks),
a custom-painted alignment canvas that only draws visible cells, and background
threads for align/clean that report through `Progress`.

It also owns the library panel and the chromatogram. Two rules there are worth
stating, because both are easy to break by trying to be clever:

* **A `Document` opened from a library entry saves through the library**, not
  through `tolviewer_io::write_file`. `Document::origin` is what marks it, and
  `save_through_library` is what intercepts it. "Save as" deliberately clears
  `origin`: the user picked a destination themselves, so later saves must not
  surprise them by going somewhere else.
* **The chromatogram does not track edits; it re-derives them.**
  `TraceView::relink` finds the offset at which the row's residues sit inside
  the file's calls, every time the document's revision changes. That is correct
  after an edit, after an undo and after a redo, because nothing is remembered
  that could go stale. Do not replace it with a second undo stack for traces.
  Trimming a trace document is an ordinary `EditOp` on the row — the signal is
  evidence and is never cut — and the link follows it.
