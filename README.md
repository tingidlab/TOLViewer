# TOLViewer

A cross-platform desktop viewer and editor for DNA, RNA and protein sequences
and alignments, with a project library for the pile of files a sequencing run
leaves behind. It does the everyday work you would otherwise open Geneious
for — keep a project's reads organised, vet the base calls against the
chromatogram, trim the primers off, align what you have, concatenate the loci,
and export to whatever your phylogenetics pipeline expects — in a single
self-contained binary with no runtime to install.

Everything runs natively. Alignment and cleaning are Rust reimplementations of
the published algorithms, so there are no external executables to download,
license or keep on `$PATH`.

## What it does

**Organise**
* A project library down the left-hand side: folders and subfolders you
  arrange — "Lace bug project" with "18S" and "28S" under it — so different
  genes and sequencing batches stay apart.
* **The files stay where they are.** Adding a file to a library reads it once
  and remembers where it is; nothing is copied into a database and nothing is
  moved. A library over a read-only archive or a shared drive works fine.
* Select across folders and align the lot in one gesture, or open several
  sequences side by side.
* Saving an edit that would replace one of the sequencing facility's own files
  **asks first**, and offers to keep the edit in a copy instead. Say yes once
  and every later edit goes to that copy without asking again — so an edited
  sequence leaves one extra file behind, not a pile of them.
* Save the library as a `.tolvlib` file. Paths inside the project folder are
  stored relative to it, so the whole project can be zipped, moved or shared
  and still open.

**Sanger traces**
* Reads Applied Biosystems `.ab1` files and draws the chromatogram under the
  sequence, so a doubtful call can be looked at rather than guessed about.
* Calls the instrument was unsure of are flagged, and so are calls a person has
  since changed — retyping one is an ordinary undoable edit, and the signal
  underneath is left alone as the evidence it is.
* Reverse complement a read to account for sequencing from the reverse primer.
  The trace is flipped with it, and the file on disk is not touched.

**Primers**
* Keep the project's primers in the library, degenerate IUPAC codes and all.
* Map them onto reads to see where they bind, on either strand, with a
  mismatch budget — the start of a Sanger read is where the basecaller is
  least sure, so an exact match is not the common case.
* Trim reads back to the amplicon. The trim is an ordinary edit you can look
  at and undo, and it goes through the same save question as anything else.

**Assemble**
* Extract individual sequences out of a multiple alignment as library entries
  of their own, without duplicating anything on disk.
* Concatenate per-locus alignments into a supermatrix, matching the same
  specimen across loci by name: `TL-2213_18S_F` and `TL_2213_28S` are
  recognised as one animal. Specimens missing a locus are gapped and
  **reported**, because a row that is mostly gaps is usually a name that failed
  to match rather than a gap in the sampling.
* Partition boundaries come out as NEXUS `charset` lines or RAxML partitions,
  ready for a partitioned analysis.

**View**
* A virtualised alignment canvas: only the visible cells are drawn, so a
  5,000 x 100,000 alignment scrolls as smoothly as a small one.
* Frozen name gutter and frozen position ruler, per-column quality track,
  consensus row, and a cleaning-mask track.
* Colour by residue, Clustal X palette, differences from consensus,
  per-column conservation, or Phred quality; dots for consensus matches.
* Per-sequence ungapped length, ambiguity fraction and mean Phred score.

**Edit**
* Change, insert or delete single residues; insert and delete whole alignment
  columns; blank a rectangular selection to gaps.
* Add, remove, reorder, rename, hide and reverse-complement sequences.
* Unlimited undo/redo with exact inverses — including the row padding that
  column edits imply, which is the part most editors get wrong.
* Select by cell, by column (drag in the ruler) or by row (drag in the names).

**Align**
* Three engines, all native:
  * **Clustal** — progressive alignment on a neighbour-joining guide tree with
    sequence weighting and position-specific gap penalties, after
    Thompson *et al.* (1994).
  * **MUSCLE** — k-mer draft, tree re-estimation, then iterative refinement by
    tree-dependent restricted partitioning, after Edgar (2004).
  * **MAFFT** — FFT-accelerated group-to-group alignment (FFT-NS-2), after
    Katoh *et al.* (2002).
* Align everything, realign just the selected columns, or add sequences to an
  existing alignment. Jobs run in the background with a progress bar and a
  working Cancel button.
* Long sequences are handled: the global aligner switches to Hirschberg's
  linear-space algorithm past ~4M DP cells, so whole organelle genomes align
  without exhausting memory.

Measured on the simulated benchmarks in `crates/tolviewer-align/tests/`
(release build, 8-core desktop). Accuracy is the sum-of-pairs score against
the known true alignment:

| | Clustal | MUSCLE | MAFFT |
| --- | --- | --- | --- |
| Low-divergence DNA (SP) | 0.94–0.97 | 0.94–0.97 | 0.94–0.97 |
| Low-divergence protein (SP) | 0.970 | 0.974 | 0.972 |
| Moderate divergence, ~69% id (SP) | 0.65 | 0.68 | 0.67 |
| 200 seqs × ~1000 cols | 2.3 s | 34 s | 4.5 s |
| 2 seqs × ~20 kb | 4.0 s | 4.0 s | 8.0 s |

MUSCLE's refinement is what makes it both the most accurate on divergent sets
and much the slowest; set refinement rounds to 0 to get its draft alignment at
Clustal-like speed. The align dialog warns before a run that will take
minutes.

**Clean**
* **Gblocks** (Castresana 2000) with all five parameters exposed, the original
  defaults, and the relaxed settings of Talavera & Castresana (2007). Results
  are previewed as a track over the alignment before you commit to them.
* Simpler filters too: drop all-gap columns, drop columns over a gap
  threshold, trim ragged ends.

**Read and write**

| Format | Read | Write |
| --- | :---: | :---: |
| FASTA | yes | yes |
| FASTQ (Phred+33/+64) | yes | yes |
| PHYLIP, strict and relaxed | yes | yes |
| NEXUS | yes | yes |
| Clustal | yes | yes |
| Stockholm | yes | yes |
| MSF / GCG | yes | — |
| GenBank | yes | — |
| AB1 (Sanger trace) | yes | — |

Files are detected by content first and extension second, and saving is atomic
(write to a temporary file, then rename), so an interrupted save cannot leave
you with a truncated alignment.

## Install

Download a build for your platform from the
[releases page](https://github.com/tingidlab/TOLViewer/releases), or build from
source:

```sh
git clone https://github.com/tingidlab/TOLViewer
cd TOLViewer
cargo build --release
./target/release/tolviewer            # or target\release\tolviewer.exe
```

You need a [Rust](https://rustup.rs) toolchain of 1.95 or newer, which is what
egui requires. On Linux you also need the windowing development packages:

```sh
# Debian / Ubuntu
sudo apt install libxkbcommon-dev libwayland-dev libx11-dev libxcursor-dev \
                 libxi-dev libxrandr-dev libgl1-mesa-dev
# Fedora
sudo dnf install libxkbcommon-devel wayland-devel libX11-devel libXcursor-devel \
                 libXi-devel libXrandr-devel mesa-libGL-devel
```

macOS and Windows need nothing beyond the toolchain.

## Use

```sh
tolviewer                       # start empty
tolviewer alignment.fasta       # open a file
tolviewer *.nex                 # open several, one per tab
tolviewer reads/*.ab1           # traces, with their chromatograms
tolviewer lace-bugs.tolvlib     # open a project library
```

You can also drop files onto the window. The library you had open last is
reopened on the next run.

### Keyboard

| | |
| --- | --- |
| arrows | move the caret |
| shift + arrows | extend the selection |
| Home / End | start / end of the row |
| PageUp / PageDown | scroll a screen |
| a letter, or `-` | overwrite the residue and step right |
| Insert | insert a gap (a whole column when columns are selected) |
| Delete | remove columns or rows; blank a cell selection to gaps |
| Backspace | as Delete, stepping left first |
| Ctrl/Cmd + Z / Shift+Z | undo / redo |
| Ctrl/Cmd + A | select all |
| Ctrl/Cmd + C | copy the selection as FASTA |
| Ctrl/Cmd + V | paste FASTA as new sequences |
| Ctrl/Cmd + O / S / Shift+S | open / save / save as |
| Ctrl/Cmd + G | go to column |
| Ctrl/Cmd + `+` / `-` | zoom |

## How it is put together

```
tolviewer-core    data model — sequences, alignments, editing, undo, statistics
tolviewer-io      readers and writers for the formats above
tolviewer-align   Clustal-, MUSCLE- and MAFFT-style alignment engines
tolviewer-clean   Gblocks and the simpler column/row filters
tolviewer-library the project library, primers, concatenation, the save policy
tolviewer-app     the egui/eframe desktop application
```

`tolviewer-core` has no dependencies at all, and the four other engine crates
depend only on each other (plus `rayon` for parallelism), so they are usable as
libraries independently of the GUI — and they build on Rust 1.85, well below
the 1.95 the GUI needs, which CI checks on every push.
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) is the API contract between them.

## A note on the reimplementations

The alignment and cleaning engines are written from the published descriptions
of the algorithms. They follow the same heuristics and will usually agree with
the original programs, but they are not forks and **they will not produce
byte-identical output**. If you are reproducing a published analysis that used
a specific version of ClustalW, MUSCLE, MAFFT or Gblocks, run that program.
If you are aligning sequences to look at them, curate them, and build a tree,
these are the same algorithms and are meant to be trusted for that.

Citations for the originals:

* Thompson JD, Higgins DG, Gibson TJ (1994) CLUSTAL W. *Nucleic Acids Res* 22:4673-4680.
* Edgar RC (2004) MUSCLE. *Nucleic Acids Res* 32:1792-1797.
* Katoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT. *Nucleic Acids Res* 30:3059-3066.
* Castresana J (2000) Selection of conserved blocks. *Mol Biol Evol* 17:540-552.
* Talavera G, Castresana J (2007) *Syst Biol* 56:564-577.

## Contributing

```sh
cargo test --workspace
cargo test -p tolviewer-align --release   # the accuracy tests need optimisation
cargo clippy --workspace --all-targets
cargo fmt --all
```

CI runs all of the above on Linux, macOS and Windows, plus a build at the
minimum supported Rust version.

## License

MIT. See [LICENSE](LICENSE).
