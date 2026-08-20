//! Turning sequence names into the sample they came from.
//!
//! Names arrive from the sequencing facility in whatever convention the lab
//! uses, and the same specimen appears under slightly different names in each
//! locus's file:
//!
//! ```text
//! TL-2213_18S_F.ab1     TL-2213_28S_F      TL2213 COI-r2
//! ```
//!
//! Concatenating those alignments means deciding they are all specimen
//! TL-2213. [`sample_key`] does that by normalising punctuation and case, then
//! stripping the trailing tokens that name a locus, a primer direction or a
//! read number rather than a specimen.
//!
//! It is a heuristic, so the GUI shows what it decided before anything is
//! concatenated, and [`MatchOptions`] lets a lab whose convention it gets wrong
//! turn the stripping off or add its own tokens.

/// How aggressively to reduce a name to its sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatchOptions {
    /// Fold case, and treat `-`, `.`, space and `_` as the same separator.
    pub normalize: bool,
    /// Strip trailing locus, direction and read-number tokens.
    pub strip_suffixes: bool,
    /// Extra trailing tokens to strip, in addition to [`DEFAULT_SUFFIXES`].
    /// Compared case-insensitively.
    pub extra_suffixes: Vec<String>,
}

impl Default for MatchOptions {
    fn default() -> Self {
        MatchOptions { normalize: true, strip_suffixes: true, extra_suffixes: Vec::new() }
    }
}

impl MatchOptions {
    /// Match on the name exactly as written. Use this when the lab's names are
    /// already the sample ids and any cleverness would only cause collisions.
    pub fn exact() -> Self {
        MatchOptions { normalize: false, strip_suffixes: false, extra_suffixes: Vec::new() }
    }
}

/// Trailing tokens that describe the read rather than the specimen.
///
/// Read direction, the common ribosomal and barcoding loci, and the plate/read
/// numbering that follows them. A token is only stripped when something is left
/// in front of it, so a sequence genuinely called "COI" keeps its name.
pub const DEFAULT_SUFFIXES: &[&str] = &[
    // Direction.
    "f",
    "r",
    "fw",
    "rv",
    "fwd",
    "rev",
    "forward",
    "reverse",
    "5",
    "3",
    // Loci a phylogenetics lab sequences by the plateful.
    "18s",
    "28s",
    "16s",
    "12s",
    "5.8s",
    "its",
    "its1",
    "its2",
    "coi",
    "co1",
    "cox1",
    "coii",
    "cytb",
    "ef1a",
    "h3",
    "wg",
    "wingless",
    "28sd2",
    "28sd3",
    "rag1",
    "nd1",
    "nd2",
    "atp6",
    // Housekeeping.
    "seq",
    "sequence",
    "consensus",
    "contig",
    "trimmed",
    "edited",
    "copy",
];

/// File extensions to drop before anything else, since names are often just
/// the trace file's name.
const EXTENSIONS: &[&str] = &[
    "ab1", "abi", "fsa", "fasta", "fa", "fas", "fna", "faa", "seq", "fastq", "fq", "phy", "nex",
    "aln", "sto", "gb", "txt",
];

/// The sample key for a name: what two names must share to be the same
/// specimen.
///
/// `name` is a sequence *id* — the first field of a FASTA header, or the name
/// column of a PHYLIP or NEXUS matrix. Strict PHYLIP names may contain spaces,
/// so whitespace is treated as one more separator rather than as the end of the
/// name.
///
/// Never returns an empty key for a non-empty name: a name that normalises or
/// strips away to nothing falls back to itself, so two unrelated oddly-punctuated
/// names cannot collide on `""` and be declared the same specimen.
pub fn sample_key(name: &str, opts: &MatchOptions) -> String {
    let trimmed = name.trim();
    if !opts.normalize {
        let key = if opts.strip_suffixes { strip(trimmed, opts) } else { trimmed.to_string() };
        return non_empty(key, trimmed);
    }
    let normalized = normalize(trimmed);
    let key = if opts.strip_suffixes { strip(&normalized, opts) } else { normalized };
    non_empty(key, trimmed)
}

/// `key` unless it came out empty, in which case the name stands for itself.
fn non_empty(key: String, original: &str) -> String {
    if key.is_empty() {
        original.to_string()
    } else {
        key
    }
}

/// Lowercase, and reduce every run of separator characters to a single `_`.
fn normalize(name: &str) -> String {
    let lower = name.to_ascii_lowercase();
    let mut out = String::with_capacity(lower.len());
    let mut pending = false;
    for c in lower.chars() {
        if matches!(c, '_' | '-' | ' ' | '.' | '|' | ':' | '/' | '\\' | '\t') {
            // A separator only counts once it has something to separate.
            pending = !out.is_empty();
        } else {
            if pending {
                out.push('_');
                pending = false;
            }
            out.push(c);
        }
    }
    out
}

/// Drop trailing tokens that name the read rather than the specimen.
fn strip(name: &str, opts: &MatchOptions) -> String {
    let mut parts: Vec<&str> = name.split('_').filter(|p| !p.is_empty()).collect();
    if parts.is_empty() {
        return name.trim().to_string();
    }
    // The last token is often a file extension left on the name.
    if parts.len() > 1 {
        if let Some(last) = parts.last() {
            if EXTENSIONS.contains(&last.to_ascii_lowercase().as_str()) {
                parts.pop();
            }
        }
    }
    while parts.len() > 1 {
        let last = parts[parts.len() - 1].to_ascii_lowercase();
        let is_suffix = DEFAULT_SUFFIXES.contains(&last.as_str())
            || opts.extra_suffixes.iter().any(|s| s.eq_ignore_ascii_case(&last))
            // A bare number is a read, plate or replicate index.
            || (last.chars().all(|c| c.is_ascii_digit()) && last.len() <= 3)
            // "f2", "r1", "rev3": a direction with a replicate number.
            || is_direction_with_number(&last);
        if !is_suffix {
            break;
        }
        parts.pop();
    }
    let joined = parts.join("_");
    if joined.is_empty() {
        name.trim().to_string()
    } else {
        joined
    }
}

/// `f2`, `r10`, `fwd3` — a direction marker with a replicate number stuck to it.
fn is_direction_with_number(token: &str) -> bool {
    let split = token.find(|c: char| c.is_ascii_digit());
    match split {
        Some(0) | None => false,
        Some(i) => {
            let (word, number) = token.split_at(i);
            number.chars().all(|c| c.is_ascii_digit())
                && matches!(word, "f" | "r" | "fw" | "rv" | "fwd" | "rev")
        }
    }
}

/// Group names by sample key, keeping the order they were first seen in.
///
/// The returned display name for each group is the first name that mapped to
/// it, which is what the concatenated alignment is labelled with.
pub fn group<'a>(
    names: impl IntoIterator<Item = &'a str>,
    opts: &MatchOptions,
) -> Vec<(String, Vec<&'a str>)> {
    let mut groups: Vec<(String, Vec<&'a str>)> = Vec::new();
    for name in names {
        let key = sample_key(name, opts);
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, members)) => members.push(name),
            None => groups.push((key, vec![name])),
        }
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(name: &str) -> String {
        sample_key(name, &MatchOptions::default())
    }

    #[test]
    fn one_specimen_across_three_loci_gets_one_key() {
        let names = ["TL-2213_18S_F.ab1", "TL_2213_28S", "TL 2213 COI-r2", "tl.2213"];
        let keys: Vec<String> = names.iter().map(|n| key(n)).collect();
        assert!(keys.windows(2).all(|w| w[0] == w[1]), "{keys:?}");
        assert_eq!(keys[0], "tl_2213");
    }

    #[test]
    fn different_specimens_stay_apart() {
        assert_ne!(key("TL-2213_18S_F"), key("TL-2214_18S_F"));
        assert_ne!(key("Corythucha_ciliata"), key("Corythucha_arcuata"));
    }

    #[test]
    fn whitespace_is_a_separator_not_the_end_of_the_name() {
        // Strict PHYLIP allows spaces in a name, so "TL 2213" is one specimen,
        // written the way "TL-2213" and "TL_2213" are elsewhere.
        assert_eq!(key("TL 2213 18S F"), "tl_2213");
        assert_eq!(key("TL 2213"), key("TL-2213"));
    }

    #[test]
    fn a_name_that_is_all_suffix_keeps_itself() {
        // Stripping everything would collapse these into one sample.
        assert_ne!(key("COI"), key("18S"));
        assert_eq!(key("COI"), "coi");
        // Punctuation-only names normalise to nothing, so they stand for
        // themselves rather than all collapsing onto one key.
        assert_eq!(key("_"), "_");
        assert_ne!(key("_"), key("--"));
        assert_eq!(key(""), "");
    }

    #[test]
    fn direction_markers_with_replicate_numbers_are_stripped() {
        assert_eq!(key("TL2213_f2"), "tl2213");
        assert_eq!(key("TL2213_rev10"), "tl2213");
        // Not a direction: keep it, or "sampleA" and "sampleB" would merge.
        assert_eq!(key("TL2213_x2"), "tl2213_x2");
        // A leading digit is part of the token, not a replicate number.
        assert!(!is_direction_with_number("2f"));
    }

    #[test]
    fn long_numbers_are_specimen_ids_not_read_numbers() {
        // Three digits or fewer reads as a replicate; more is an accession.
        assert_eq!(key("sample_12"), "sample");
        assert_eq!(key("sample_12345"), "sample_12345");
    }

    #[test]
    fn exact_matching_leaves_names_alone() {
        let o = MatchOptions::exact();
        assert_eq!(sample_key("TL-2213_18S_F", &o), "TL-2213_18S_F");
        assert_ne!(sample_key("TL-2213_18S_F", &o), sample_key("TL_2213_18S_F", &o));
    }

    #[test]
    fn normalising_without_stripping_unifies_punctuation_only() {
        let o = MatchOptions { normalize: true, strip_suffixes: false, ..Default::default() };
        assert_eq!(sample_key("TL-2213_18S_F", &o), "tl_2213_18s_f");
        assert_eq!(sample_key("TL.2213 18S/F", &o), "tl_2213_18s_f");
    }

    #[test]
    fn a_lab_can_add_its_own_suffixes() {
        let o = MatchOptions {
            extra_suffixes: vec!["plateA".to_string(), "run7".to_string()],
            ..Default::default()
        };
        assert_eq!(sample_key("TL2213_plateA_run7", &o), "tl2213");
        // Without them, they are part of the sample.
        assert_eq!(key("TL2213_plateA_run7"), "tl2213_platea_run7");
    }

    #[test]
    fn grouping_reports_members_in_order() {
        let names = ["a_18S_F", "b_18S_F", "a_28S_R", "a_COI"];
        let groups = group(names, &MatchOptions::default());
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].0, "a");
        assert_eq!(groups[0].1, vec!["a_18S_F", "a_28S_R", "a_COI"]);
        assert_eq!(groups[1].1, vec!["b_18S_F"]);
    }

    #[test]
    fn leading_and_repeated_separators_do_not_create_empty_tokens() {
        assert_eq!(key("__TL--2213__18S__"), "tl_2213");
        assert_eq!(normalize("--a--b--"), "a_b");
    }
}
