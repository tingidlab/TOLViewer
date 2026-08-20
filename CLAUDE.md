# TOLViewer

A cross-platform desktop viewer/editor for DNA, RNA and protein alignments.
Rust workspace, egui GUI, five crates, no runtime dependencies and no external
executables.

**Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) before changing anything.**
It holds the crate layout, the ground rules (no `unsafe`, no panics on user
input, case-preserving residues, interruptible long work) and the inter-crate
API contract. Those public signatures are fixed; changing one means changing
that document in the same commit.

## Commands

```sh
cargo test --workspace                      # 330 tests, ~1 min
cargo test -p tolviewer-align --release -- --ignored   # alignment accuracy + timing
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all
cargo run -- testdata/examples/tingidae_COI.fasta      # DNA, unaligned
cargo run -- testdata/examples/wingless_protein.fasta  # protein
```

CI runs all of the above on Linux, macOS and Windows. Everything must be green
before pushing; the workflows are the source of truth for what "done" means.

## Two-tier MSRV

This is CI-enforced and easy to break by accident:

* **`tolviewer-app` needs Rust 1.95**, because egui 0.36 does.
* **`tolviewer-core`, `-io`, `-align`, `-clean` build on Rust 1.85.** They
  depend on nothing but each other and `rayon`, and the README offers them as
  libraries usable without the GUI.

So a newer-than-1.85 language feature (`is_multiple_of`, let-else in a new
position, a fresh `std` method) is fine in the app crate and a CI failure in the
four library crates. Verify with `cargo +1.85.0 build -p tolviewer-core -p
tolviewer-io -p tolviewer-align -p tolviewer-clean`.

`tolviewer-core` has **zero** dependencies. Keep it that way.

## Things that will bite you

**egui 0.36 is the Ui-based app model, not the older Context-based one.**
`App::ui(&mut self, ui, frame)`, and panels are the unified
`egui::Panel::top(egui::Id::new("menu")).show(ui, ..)` — not `TopBottomPanel` /
`SidePanel`. Most egui examples you'll find online are for 0.31 and won't
compile. `ctx.available_rect()` and `InputState::screen_rect()` are gone.

**Headless canvas tests must call `output.textures_delta.clear()`.** egui
panics on drop otherwise. See the `paint()` helper in
`crates/tolviewer-app/tests/canvas_render.rs:35`.

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

## Rendering

The canvas is virtualised: only visible cells are painted, equal-colour runs are
coalesced into single rects, and glyphs are dropped below 6pt cell width. This
is a correctness property, not an optimisation —
`a_huge_alignment_costs_no_more_to_paint_than_a_screenful` pins it
quantitatively (4000x40000 must cost no more than 4x a 40x400 alignment). If you
touch `canvas.rs`, keep that test passing.

Rendering can't be checked visually here (Wayland; no working screenshot path),
so it's verified through headless `Context::run_ui` tests instead. That is how
rendering bugs get caught — add to those tests rather than eyeballing it.

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

Two known gaps, both deliberate for 0.1.0 and both worth fixing before this is
handed to people who did not build it:

* **Nothing is signed or notarised.** macOS Gatekeeper refuses the binary on
  first launch unless the user does right-click → Open or strips the quarantine
  attribute, and SmartScreen warns on Windows. The v0.1.0 release notes say so
  explicitly; keep saying so until it is fixed, because a user who hits this
  silently concludes the download is broken.
* **The Linux binary is not stripped.** It ships at 7.8 MB against Windows'
  5.6 MB, almost entirely debug symbols. `strip = true` in `[profile.release]`
  roughly halves it, at the cost of useful backtraces in bug reports — which is
  the actual trade-off to think about, not a free win.

## Conventions

* Tests live in `#[cfg(test)] mod tests` in the same file, plus `tests/` for
  round-trip and cross-crate integration work.
* Prefer `std` + `rayon`. New dependencies need a real justification.
* Saving is atomic (temp file + rename). Preserve that.
* Commit messages explain *why*, in prose, wrapped at 72 columns.
