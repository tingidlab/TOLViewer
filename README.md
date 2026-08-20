# TOLViewer

A cross-platform desktop viewer and editor for DNA, RNA and protein sequences
and alignments. It does the everyday work you would otherwise open Geneious
for — look at an alignment, check whether it is any good, fix the places where
it is not, align it, trim it, and export it to whatever your phylogenetics
pipeline expects — in a single self-contained binary with no runtime to
install.

Everything runs natively. Alignment and cleaning are Rust reimplementations of
the published algorithms, so there are no external executables to download,
license or keep on `$PATH`.

## What it does

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
```

You can also drop files onto the window.

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
tolviewer-app     the egui/eframe desktop application
```

`tolviewer-core` has no dependencies at all, and the three engine crates depend
only on it (plus `rayon` for parallelism), so they are usable as libraries
independently of the GUI — and they build on Rust 1.85, well below the 1.95
the GUI needs, which CI checks on every push.
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
