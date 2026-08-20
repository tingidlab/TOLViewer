# TOLViewer

A cross-platform desktop viewer/editor for DNA, RNA and protein alignments,
with a project library over the files a sequencing run produces. Rust
workspace, egui GUI, six crates, no runtime dependencies and no external
executables.

**Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) before changing anything.**
It holds the crate layout, the ground rules (no `unsafe`, no panics on user
input, case-preserving residues, interruptible long work) and the inter-crate
API contract. Those public signatures are fixed; changing one means changing
that document in the same commit.

## Commands

```sh
cargo test --workspace                      # 466 tests, ~1 min
cargo test -p tolviewer-align --release -- --ignored   # alignment accuracy + timing
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all
cargo run -- testdata/examples/tingidae_COI.fasta      # DNA, unaligned
cargo run -- testdata/examples/wingless_protein.fasta  # protein
cargo run -- testdata/ab1/tingidae_COI_F.ab1           # trace + chromatogram
cargo run -- <project>.tolvlib                         # opens as a library
```

CI runs all of the above on Linux, macOS and Windows. Everything must be green
before pushing; the workflows are the source of truth for what "done" means.

## Two-tier MSRV

This is CI-enforced and easy to break by accident:

* **`tolviewer-app` needs Rust 1.95**, because egui 0.36 does.
* **`tolviewer-core`, `-io`, `-align`, `-clean`, `-library` build on Rust
  1.85.** They depend on nothing but each other and `rayon`, and the README
  offers them as libraries usable without the GUI.

So a newer-than-1.85 language feature (`is_multiple_of`, let-else in a new
position, a fresh `std` method) is fine in the app crate and a CI failure in the
five library crates. Verify with `cargo +1.85.0 build -p tolviewer-core -p
tolviewer-io -p tolviewer-align -p tolviewer-clean -p tolviewer-library`.

`tolviewer-core` has **zero** dependencies. Keep it that way.

## Things that will bite you

**egui 0.36 is the Ui-based app model, not the older Context-based one.**
`App::ui(&mut self, ui, frame)`, and panels are the unified
`egui::Panel::top(egui::Id::new("menu")).show(ui, ..)` — not `TopBottomPanel` /
`SidePanel`. Most egui examples you'll find online are for 0.31 and won't
compile. `ctx.available_rect()` and `InputState::screen_rect()` are gone.

**Headless render tests must call `output.textures_delta.clear()`.** egui
panics on drop otherwise. See the `paint()` helper in
`crates/tolviewer-app/tests/canvas_render.rs:35`, and the same pattern in
`tests/trace_render.rs`.

**Undo inverses must be exact, including side effects.** Column operations pad
ragged rows via `pad_to_width()`, and that padding is part of what undo has to
reverse — `UndoStack::apply` records the row shape first and wraps the inverse
in `InverseOp::Compound` when it changed
(`crates/tolviewer-core/src/edit.rs:173`). Any new `EditOp` that can alter row
lengths needs the same treatment, and a round-trip test asserting byte identity
after undo.

**PHYLIP, NEXUS, Clustal and Stockholm cannot represent ragged rows.** Writing
an unaligned set to them fails with `Error::Format`; the app pads and retries,
then tells the user it did. Don't turn that into a bare error.

**The alignment accuracy tests are meaningless in debug builds** — they're
`#[ignore]`d and need `--release`, or they take minutes.

**The library never writes to the lab's files without being asked, and that is
the whole point of it.** `Entry::save_target` returns `Original` for a file the
user imported and `MustCopy` for one that cannot be written back at all; the
GUI turns both into a question. `SaveChoice::NewCopy` then *remembers* the copy
so the question is asked once per entry, not once per save. If you add a code
path that writes an entry, route it through `Library::save_entry` — anything
that reaches `tolviewer_io::write_file` directly has escaped the policy.
`crates/tolviewer-library/tests/workflow.rs` asserts byte-for-byte that a full
session leaves every original untouched; keep that passing.

**Any successful save clears `Entry::select` and `Entry::reversed`.** Those two
say how to turn the file's contents into what the entry shows, and after a save
the file *is* what the entry shows. Leaving `reversed` set reverse-complements
the saved sequence a second time on the next read — which is silent, plausible
and wrong. This bit once already; `overwriting_a_reversed_entry_does_not_reverse_it_again`
pins it.

**The chromatogram re-derives its link to the sequence, it does not track it.**
`TraceView::relink` finds the offset at which the row's residues sit inside the
file's calls, whenever the document revision changes. It looks wasteful and
isn't: it is correct after an edit, an undo *and* a redo, which no amount of
bookkeeping alongside the undo stack would be. The app never trims a
`TraceView` — trimming a trace document is an `EditOp` on the row, and the
signal stays as evidence of what the instrument saw. (`Trace::trim` does exist
in `tolviewer-io`, for cropping a standalone trace; the GUI does not use it, and
wiring it in would break the relinking above.)

**AB1 is the only binary format, and `Format::sniff` settles it first**, before
anything is decoded as text. Files in the wild disagree about tag lengths, so
the reader trims to the shortest consistent set rather than refusing a read.
The test files in `testdata/ab1/` are generated by `testdata/ab1/make_ab1.py`;
regenerate rather than hand-editing them.

**Name matching for concatenation is a heuristic and is allowed to be wrong.**
`naming::sample_key` strips locus and direction suffixes, which is why
`TL-2213_18S_F` and `TL_2213_28S` become one specimen. `concat::preview` exists
so the GUI can show what it decided *before* building a matrix, because a
mismatch is invisible in the result and obvious in the table. Never let
`sample_key` return an empty string for a non-empty name — every oddly
punctuated name would collide on it.

## Rendering

The canvas is virtualised: only visible cells are painted, equal-colour runs are
coalesced into single rects, and glyphs are dropped below 6pt cell width. This
is a correctness property, not an optimisation —
`a_huge_alignment_costs_no_more_to_paint_than_a_screenful` pins it
quantitatively (4000x40000 must cost no more than 4x a 40x400 alignment). If you
touch `canvas.rs`, keep that test passing.

The chromatogram has the same rule and for the same reason: below
`MIN_LABEL_WIDTH` (6pt per call, matching the canvas) the base letters are
dropped, so painting a 1000-base read costs no more than painting a window of
it. `painting_a_window_costs_no_more_than_the_window_holds` pins that.

Rendering can't be checked visually here (Wayland; no working screenshot path),
so it's verified through headless `Context::run_ui` tests instead —
`canvas_render.rs` for the alignment, `trace_render.rs` for the chromatogram.
That is how rendering bugs get caught, and all four of the bugs found while
building 0.2.0 were found that way. Add to those tests rather than eyeballing
it.

## The engine reimplementations

`tolviewer-align` and `tolviewer-clean` are written from the published papers
(Thompson 1994, Edgar 2004, Katoh 2002, Castresana 2000), not ported from the
GPL C sources. They follow the same heuristics but **are not byte-identical** to
the original programs, and the README says so. Don't claim equivalence, and
don't add code that shells out to a `clustalw`/`muscle`/`mafft` binary — the
whole point is a single self-contained executable.

Accuracy is measured as sum-of-pairs against known true alignments in
`crates/tolviewer-align/tests/accuracy.rs`. Thresholds there are measured values
rounded down; if a change moves them, that's a real regression to explain, not a
number to adjust.

## Releases

Tagging `v*` triggers `.github/workflows/release.yml`, which builds
x86_64-linux-gnu, both macOS architectures and x86_64-windows-msvc, then
publishes a GitHub release with the four archives attached. Check CI is green on
the commit first — the release workflow does not re-run the tests.

`generate_release_notes: true` only produces anything useful once there is a
previous tag to diff against. For v0.1.0 it emitted a bare "Full Changelog"
link, and the notes were written by hand afterwards with `gh release edit
<tag> --notes-file`. Expect to do that again for anything worth announcing.

Both 0.1.0 gaps are now dealt with, one fully and one waiting on certificates:

* **Signing is wired up but needs secrets.** `release.yml` has a macOS
  codesign + notarytool step and a Windows signtool step. Both are skipped
  cleanly when the secrets are absent, so the workflow still produces a release
  either way — but an unsigned build is refused by Gatekeeper on first launch
  and warned about by SmartScreen, and **a user who hits that concludes the
  download is broken rather than unsigned**. So: until `MACOS_CERTIFICATE` and
  `WINDOWS_CERTIFICATE` (and the other secrets listed in the workflow) are
  configured, the release notes must keep saying the binaries are unsigned and
  how to get past it. Check the workflow log for whether the signing steps
  actually ran before writing notes that claim otherwise.
* **Stripping is settled: `strip = "debuginfo"`.** Measured on Linux with the
  release profile: 24.6 MB unstripped, 21.7 MB as configured, 18.5 MB with
  `strip = true`. The old note in this file claimed the binary was "almost
  entirely debug symbols" and that stripping would halve it; that was wrong.
  It is mostly `.text` — 13.6 MB of egui and the alignment engines. The 3 MB
  between the current setting and `true` is the symbol table, which is what
  makes a backtrace in a bug report name functions instead of listing
  addresses. Don't trade that for 3 MB, and don't expect a big win from
  stripping further.

## The library

`tolviewer-library` is the newest crate and the one whose *design* matters more
than its code. A library is an index over files that stay where the sequencing
facility put them: `add_file` reads a file once to see what it holds and then
remembers only where it is. Nothing is copied, moved or rewritten by filing it.

Two things are stored as flags rather than as edits, deliberately:

* `Entry::reversed` — flipping a read to account for a reverse primer is cheap,
  lossless and constantly undone. Rewriting the file for it would be slower and
  destructive.
* `Entry::select` — an entry extracted from an alignment names the rows it
  stands for and reads them back out of the same file. That is why extraction
  costs nothing on disk, and why saving an edit to one insists on a copy.

The on-disk format (`store.rs`) skips unknown keys inside known blocks, so a
library written by a later version still opens. Bump `VERSION` only when an old
reader could *misread* a new file; adding keys does not need it.

## Conventions

* Tests live in `#[cfg(test)] mod tests` in the same file, plus `tests/` for
  round-trip and cross-crate integration work.
* Prefer `std` + `rayon`. New dependencies need a real justification.
* Saving is atomic (temp file + rename). Preserve that.
* Commit messages explain *why*, in prose, wrapped at 72 columns.
