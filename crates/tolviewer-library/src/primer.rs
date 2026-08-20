//! Mapping PCR primers onto reads, and trimming reads back to the amplicon.
//!
//! A Sanger read starts inside the forward primer and, if it ran far enough,
//! ends inside the reverse primer's binding site. Both are laboratory artefacts
//! rather than sequence from the specimen, and both are the least reliable part
//! of the read, so they come off before anything is aligned.
//!
//! Primers are written in IUPAC, and so are the reads, so matching compares
//! *code sets*: `R` in the primer matches `A`, `G` or `R` in the read, and an
//! `N` in the read matches anything. Mismatches are allowed up to a budget,
//! because a primer sits at the noisy start of the trace where the basecaller
//! is least sure.
//!
//! Search is an exhaustive scan of every offset, which is O(read x primer).
//! With reads of ~1 kb and primers of ~25 bases that is tens of thousands of
//! comparisons per primer — far too little to be worth an index.

use std::ops::Range;

use tolviewer_core::{is_gap, Alphabet, Error, Result};

/// Which strand of the read a primer matched on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// The primer sequence itself was found: a forward primer at the 5' end.
    Forward,
    /// The primer's reverse complement was found, which is how a reverse
    /// primer appears at the 3' end of a forward read.
    Reverse,
}

impl Strand {
    pub fn name(self) -> &'static str {
        match self {
            Strand::Forward => "forward",
            Strand::Reverse => "reverse",
        }
    }
}

/// One PCR primer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Primer {
    pub name: String,
    /// IUPAC nucleotide codes, uppercase. Gaps and whitespace are dropped on
    /// construction.
    pub sequence: Vec<u8>,
}

impl Primer {
    /// Build a primer, rejecting anything that is not nucleotide IUPAC.
    ///
    /// Whitespace and gap characters are stripped, so a sequence pasted out of
    /// a supplier's order form works as typed.
    pub fn new(name: impl Into<String>, sequence: &str) -> Result<Primer> {
        let mut residues = Vec::with_capacity(sequence.len());
        for c in sequence.bytes() {
            if c.is_ascii_whitespace() || is_gap(c) {
                continue;
            }
            let up = c.to_ascii_uppercase();
            if !Alphabet::Dna.is_valid(up) {
                return Err(Error::format(format!(
                    "'{}' is not a nucleotide code, so it cannot be part of primer {}",
                    c as char,
                    name.into()
                )));
            }
            residues.push(up);
        }
        let name = name.into();
        if residues.is_empty() {
            return Err(Error::format(format!("primer {name} has no sequence")));
        }
        Ok(Primer { name, sequence: residues })
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// The primer as it would read on the opposite strand.
    pub fn reverse_complement(&self) -> Vec<u8> {
        self.sequence.iter().rev().map(|&c| Alphabet::Dna.complement(c)).collect()
    }

    /// How many mismatches this primer is allowed against a read, given a
    /// tolerated fraction. Always at least 0 and never the whole primer, so a
    /// generous fraction cannot make everything match.
    pub fn mismatch_budget(&self, fraction: f32) -> usize {
        let budget = (self.len() as f32 * fraction.clamp(0.0, 1.0)).floor() as usize;
        budget.min(self.len().saturating_sub(1))
    }
}

/// A place a primer was found.
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerHit {
    /// Index into the [`PrimerSet`] the hit came from.
    pub primer: usize,
    pub name: String,
    pub strand: Strand,
    /// Where the primer sits in the read, in the read's own coordinates.
    pub range: Range<usize>,
    pub mismatches: usize,
}

impl PrimerHit {
    /// Fraction of the primer that matched, for ranking and for display.
    pub fn identity(&self) -> f32 {
        let len = self.range.len();
        if len == 0 {
            return 0.0;
        }
        (len - self.mismatches.min(len)) as f32 / len as f32
    }
}

/// The primers a project uses.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PrimerSet {
    primers: Vec<Primer>,
}

impl PrimerSet {
    pub fn new(primers: Vec<Primer>) -> Self {
        PrimerSet { primers }
    }

    pub fn primers(&self) -> &[Primer] {
        &self.primers
    }

    pub fn len(&self) -> usize {
        self.primers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.primers.is_empty()
    }

    pub fn push(&mut self, primer: Primer) {
        self.primers.push(primer);
    }

    pub fn remove(&mut self, index: usize) -> Option<Primer> {
        (index < self.primers.len()).then(|| self.primers.remove(index))
    }

    pub fn get(&self, index: usize) -> Option<&Primer> {
        self.primers.get(index)
    }

    pub fn find_by_name(&self, name: &str) -> Option<usize> {
        self.primers.iter().position(|p| p.name.eq_ignore_ascii_case(name))
    }

    /// Every place any primer matches `seq` within its mismatch budget, on
    /// either strand, best first.
    ///
    /// Overlapping hits of the same primer on the same strand are collapsed to
    /// the best one, so a primer with an internal repeat reports one binding
    /// site rather than a smear.
    pub fn map(&self, seq: &[u8], max_mismatch_fraction: f32) -> Vec<PrimerHit> {
        let mut hits = Vec::new();
        for (i, primer) in self.primers.iter().enumerate() {
            let budget = primer.mismatch_budget(max_mismatch_fraction);
            for (strand, pattern) in [
                (Strand::Forward, primer.sequence.clone()),
                (Strand::Reverse, primer.reverse_complement()),
            ] {
                for (start, mismatches) in scan(&pattern, seq, budget) {
                    hits.push(PrimerHit {
                        primer: i,
                        name: primer.name.clone(),
                        strand,
                        range: start..start + pattern.len(),
                        mismatches,
                    });
                }
            }
        }
        // Fewest mismatches first, then leftmost, so the ordering is total and
        // does not depend on the order the primers were entered.
        hits.sort_by(|a, b| {
            a.mismatches
                .cmp(&b.mismatches)
                .then(a.range.start.cmp(&b.range.start))
                .then(a.primer.cmp(&b.primer))
        });
        dedup_overlaps(hits)
    }
}

/// Drop any hit that overlaps a better hit of the same primer on the same
/// strand. `hits` must already be sorted best-first.
fn dedup_overlaps(hits: Vec<PrimerHit>) -> Vec<PrimerHit> {
    let mut kept: Vec<PrimerHit> = Vec::with_capacity(hits.len());
    for hit in hits {
        let shadowed = kept.iter().any(|k| {
            k.primer == hit.primer
                && k.strand == hit.strand
                && k.range.start < hit.range.end
                && hit.range.start < k.range.end
        });
        if !shadowed {
            kept.push(hit);
        }
    }
    kept
}

/// Every offset where `pattern` matches `seq` with at most `budget`
/// mismatches, as `(start, mismatches)`.
fn scan(pattern: &[u8], seq: &[u8], budget: usize) -> Vec<(usize, usize)> {
    if pattern.is_empty() || seq.len() < pattern.len() {
        return Vec::new();
    }
    let mut out = Vec::new();
    for start in 0..=seq.len() - pattern.len() {
        let mut mismatches = 0;
        let mut ok = true;
        for (i, &p) in pattern.iter().enumerate() {
            if !compatible(p, seq[start + i]) {
                mismatches += 1;
                if mismatches > budget {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            out.push((start, mismatches));
        }
    }
    out
}

/// Do two IUPAC codes have any base in common?
///
/// A gap in the read never matches: a primer cannot bind a column that is not
/// there, and treating gaps as wildcards would let a primer "match" an
/// all-gap stretch of an aligned row.
fn compatible(a: u8, b: u8) -> bool {
    if is_gap(b) || is_gap(a) {
        return false;
    }
    code(a) & code(b) != 0
}

/// The set of bases an IUPAC code stands for, as a four-bit mask (A C G T).
fn code(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101, // A G
        b'Y' => 0b1010, // C T
        b'S' => 0b0110, // C G
        b'W' => 0b1001, // A T
        b'K' => 0b1100, // G T
        b'M' => 0b0011, // A C
        b'B' => 0b1110, // C G T
        b'D' => 0b1101, // A G T
        b'H' => 0b1011, // A C T
        b'V' => 0b0111, // A C G
        b'N' | b'X' | b'?' => 0b1111,
        // Anything else is not a base and matches nothing, so an unexpected
        // character counts as a mismatch instead of silently matching.
        _ => 0,
    }
}

/// How to decide where a read's amplicon starts and ends.
#[derive(Debug, Clone, PartialEq)]
pub struct TrimOptions {
    /// Fraction of a primer allowed to mismatch. 0.2 tolerates five errors in a
    /// 25-mer, which is about what the start of a Sanger trace produces.
    pub max_mismatch_fraction: f32,
    /// Only look for the 5' primer in this many leading bases, and the 3'
    /// primer in this many trailing bases. A primer-length match in the middle
    /// of a read is a repeat, not a binding site. 0 searches the whole read.
    pub search_window: usize,
    /// Keep the primer sequence in the read instead of cutting it off.
    pub keep_primers: bool,
}

impl Default for TrimOptions {
    fn default() -> Self {
        TrimOptions { max_mismatch_fraction: 0.2, search_window: 120, keep_primers: false }
    }
}

/// What trimming a read would do to it.
#[derive(Debug, Clone, PartialEq)]
pub struct TrimPlan {
    /// The bases to keep.
    pub range: Range<usize>,
    /// The hit that set the 5' cut, if one was found.
    pub start_hit: Option<PrimerHit>,
    /// The hit that set the 3' cut, if one was found.
    pub end_hit: Option<PrimerHit>,
}

impl TrimPlan {
    /// Would this actually change the read?
    pub fn trims_anything(&self) -> bool {
        self.start_hit.is_some() || self.end_hit.is_some()
    }

    /// A one-line account of what was found, for the log and the tree's
    /// hover text.
    pub fn describe(&self, original_len: usize) -> String {
        match (&self.start_hit, &self.end_hit) {
            (None, None) => "no primer found".to_string(),
            _ => {
                let mut parts = Vec::new();
                if let Some(h) = &self.start_hit {
                    parts.push(format!("5' {} ({} mismatch)", h.name, h.mismatches));
                }
                if let Some(h) = &self.end_hit {
                    parts.push(format!("3' {} ({} mismatch)", h.name, h.mismatches));
                }
                format!("{}: {} -> {} bases", parts.join(", "), original_len, self.range.len())
            }
        }
    }
}

/// Work out which part of `seq` is amplicon rather than primer.
///
/// The 5' cut comes from the best hit whose *start* is inside the leading
/// window; the 3' cut from the best hit whose *end* is inside the trailing one.
/// A read short enough for the two windows to overlap is searched as a whole,
/// and the two cuts are never allowed to cross.
pub fn plan_trim(set: &PrimerSet, seq: &[u8], opts: &TrimOptions) -> TrimPlan {
    let len = seq.len();
    let hits = set.map(seq, opts.max_mismatch_fraction);
    let window = if opts.search_window == 0 { len } else { opts.search_window.min(len) };

    // `hits` is best-first, so the first match in each window is the one to use.
    let start_hit = hits.iter().find(|h| h.range.start < window).cloned();
    // The 3' primer has to be a *different* binding site from the 5' one, or a
    // read consisting of a single primer would be trimmed from both ends at
    // once and vanish.
    let end_hit = hits
        .iter()
        .find(|h| {
            h.range.end + window >= len
                && start_hit.as_ref().is_none_or(|s| h.range.start >= s.range.end)
        })
        .cloned();

    let mut start = match (&start_hit, opts.keep_primers) {
        (Some(h), false) => h.range.end,
        (Some(h), true) => h.range.start,
        (None, _) => 0,
    };
    let mut end = match (&end_hit, opts.keep_primers) {
        (Some(h), false) => h.range.start,
        (Some(h), true) => h.range.end,
        (None, _) => len,
    };
    // Two primers that overlap each other would otherwise ask for a negative
    // length. Leave the read alone rather than emptying it.
    if start >= end {
        start = 0;
        end = len;
        return TrimPlan { range: start..end, start_hit: None, end_hit: None };
    }
    TrimPlan { range: start..end, start_hit, end_hit }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Real barcoding primers: Folmer's COI pair, one of which is degenerate.
    const LCO: &str = "GGTCAACAAATCATAAAGATATTGG";
    const HCO: &str = "TAAACTTCAGGGTGACCAAAAAATCA";

    fn set() -> PrimerSet {
        PrimerSet::new(vec![
            Primer::new("LCO1490", LCO).unwrap(),
            Primer::new("HCO2198", HCO).unwrap(),
        ])
    }

    fn rc(s: &str) -> String {
        s.bytes().rev().map(|c| Alphabet::Dna.complement(c) as char).collect()
    }

    /// A read: junk, forward primer, insert, reverse primer's complement, junk.
    fn amplicon(insert: &str) -> Vec<u8> {
        format!("NNCTG{LCO}{insert}{}GGNTA", rc(HCO)).into_bytes()
    }

    #[test]
    fn a_primer_rejects_non_nucleotide_sequence() {
        assert!(Primer::new("p", "GGTCAACAAATCATAAAGATATTGG").is_ok());
        assert!(Primer::new("p", "GGT CAA-CAA").is_ok(), "spacing and gaps are cosmetic");
        assert!(Primer::new("p", "GGTZAA").is_err());
        assert!(Primer::new("p", "").is_err());
        assert!(Primer::new("p", "   ").is_err());
    }

    #[test]
    fn iupac_codes_match_the_bases_they_stand_for() {
        assert!(compatible(b'R', b'A') && compatible(b'R', b'G'));
        assert!(!compatible(b'R', b'C') && !compatible(b'R', b'T'));
        assert!(compatible(b'N', b'A') && compatible(b'A', b'N'));
        assert!(compatible(b'a', b'A'), "case is not information");
        assert!(compatible(b'T', b'U'), "an RNA read still binds a DNA primer");
        assert!(!compatible(b'A', b'-'), "a primer cannot bind a gap");
        assert!(!compatible(b'A', b'@'), "junk is a mismatch, not a wildcard");
    }

    #[test]
    fn a_degenerate_primer_finds_both_variants() {
        // Y is C or T: both reads must match with no mismatches at all.
        let primer = Primer::new("deg", "ACGTYACGT").unwrap();
        let set = PrimerSet::new(vec![primer]);
        for base in ["C", "T"] {
            let read = format!("GGGG ACGT{base}ACGT GGGG").replace(' ', "").into_bytes();
            let hits = set.map(&read, 0.0);
            assert_eq!(hits.len(), 1, "{base}");
            assert_eq!(hits[0].mismatches, 0, "{base}");
            assert_eq!(hits[0].range, 4..13);
        }
        // A is not in Y, so it costs one mismatch on the forward strand. (The
        // primer's own reverse complement, ACGTRACGT, does match this read
        // exactly — which is why the strand is part of a hit.)
        let read = b"GGGGACGTAACGTGGGG";
        assert!(!set.map(read, 0.0).iter().any(|h| h.strand == Strand::Forward));
        let forward = set.map(read, 0.2).into_iter().find(|h| h.strand == Strand::Forward).unwrap();
        assert_eq!(forward.mismatches, 1);
    }

    #[test]
    fn primers_are_found_on_both_strands() {
        let read = amplicon("ACGTACGTACGTACGTACGTACGTACGTACGT");
        let hits = set().map(&read, 0.1);
        let lco = hits.iter().find(|h| h.name == "LCO1490").expect("LCO not found");
        let hco = hits.iter().find(|h| h.name == "HCO2198").expect("HCO not found");
        assert_eq!(lco.strand, Strand::Forward);
        assert_eq!(lco.range, 5..5 + LCO.len());
        assert_eq!(hco.strand, Strand::Reverse, "the reverse primer appears complemented");
        assert_eq!(hco.range.end, read.len() - 5);
        assert_eq!(lco.mismatches, 0);
        assert_eq!(hco.mismatches, 0);
        assert_eq!(lco.identity(), 1.0);
    }

    #[test]
    fn the_mismatch_budget_scales_with_the_primer_and_never_matches_everything() {
        let p = Primer::new("p", LCO).unwrap();
        assert_eq!(p.mismatch_budget(0.0), 0);
        assert_eq!(p.mismatch_budget(0.2), 5);
        assert_eq!(p.mismatch_budget(1.0), LCO.len() - 1, "a full budget is still not free");
        assert_eq!(p.mismatch_budget(-1.0), 0);
    }

    #[test]
    fn a_read_with_a_damaged_primer_still_maps() {
        let mut read = amplicon("ACGTACGTACGTACGTACGTACGTACGTACGT");
        // The first bases of a Sanger read are the worst. LCO starts at 5 and
        // reads GGTCAACAA...; miscall three of its bases outright.
        read[5] = b'A';
        read[7] = b'A';
        read[9] = b'T';
        let hits = set().map(&read, 0.2);
        let lco = hits.iter().find(|h| h.name == "LCO1490").expect("LCO lost to three errors");
        assert_eq!(lco.mismatches, 3);
        assert_eq!(lco.range.start, 5);
        // Three errors in a 25-mer is more than a 5% budget allows.
        assert!(!set().map(&read, 0.05).iter().any(|h| h.name == "LCO1490"));
    }

    #[test]
    fn an_uncalled_base_is_not_counted_against_a_primer() {
        // N is "unknown", not "wrong": it may well be the base the primer wants,
        // and the start of a trace where the primer sits is full of them.
        let mut read = amplicon("ACGTACGTACGTACGTACGTACGTACGTACGT");
        read[5] = b'N';
        read[7] = b'N';
        let hits = set().map(&read, 0.0);
        let lco = hits.iter().find(|h| h.name == "LCO1490").expect("Ns should cost nothing");
        assert_eq!(lco.mismatches, 0);
    }

    #[test]
    fn trimming_removes_the_primers_and_the_junk_outside_them() {
        let insert = "ACGTACGTACGTTTGGAACCTTGGAACCTTAA";
        let read = amplicon(insert);
        let plan = plan_trim(&set(), &read, &TrimOptions::default());
        assert_eq!(&read[plan.range.clone()], insert.as_bytes());
        assert!(plan.trims_anything());
        assert_eq!(plan.start_hit.as_ref().unwrap().name, "LCO1490");
        assert_eq!(plan.end_hit.as_ref().unwrap().name, "HCO2198");
        assert!(plan.describe(read.len()).contains("LCO1490"), "{}", plan.describe(read.len()));
    }

    #[test]
    fn keeping_primers_cuts_only_the_junk_around_them() {
        let insert = "ACGTACGTACGTTTGGAACCTTGGAACCTTAA";
        let read = amplicon(insert);
        let opts = TrimOptions { keep_primers: true, ..Default::default() };
        let plan = plan_trim(&set(), &read, &opts);
        assert_eq!(&read[plan.range.clone()], format!("{LCO}{insert}{}", rc(HCO)).as_bytes());
    }

    #[test]
    fn a_read_with_no_primer_is_left_whole() {
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let plan = plan_trim(&set(), &read, &TrimOptions::default());
        assert_eq!(plan.range, 0..read.len());
        assert!(!plan.trims_anything());
        assert_eq!(plan.describe(read.len()), "no primer found");
    }

    #[test]
    fn a_primer_found_only_at_one_end_trims_only_that_end() {
        let insert = "ACGTACGTACGTTTGGAACCTTGGAACCTTAA";
        let read = format!("NNCTG{LCO}{insert}").into_bytes();
        let plan = plan_trim(&set(), &read, &TrimOptions::default());
        assert_eq!(plan.range, 5 + LCO.len()..read.len());
        assert!(plan.start_hit.is_some() && plan.end_hit.is_none());
    }

    #[test]
    fn the_search_window_keeps_an_internal_repeat_from_trimming_the_read_away() {
        // The forward primer's sequence also occurs deep inside the insert.
        let insert = format!("ACGTACGTACGT{LCO}TTGGAACCTTAA");
        let read = format!("NNCTG{LCO}{insert}").into_bytes();
        let narrow =
            plan_trim(&set(), &read, &TrimOptions { search_window: 40, ..Default::default() });
        assert_eq!(narrow.range.start, 5 + LCO.len(), "the 5' copy is the binding site");
        // Searching the whole read, the internal copy is a candidate for the 3'
        // cut, which is exactly why the window exists.
        let wide =
            plan_trim(&set(), &read, &TrimOptions { search_window: 0, ..Default::default() });
        assert!(wide.range.len() <= narrow.range.len());
    }

    #[test]
    fn overlapping_hits_of_one_primer_collapse_to_the_best() {
        // A homopolymer primer matches at many offsets inside a homopolymer run.
        let set = PrimerSet::new(vec![Primer::new("polyA", "AAAAAAAA").unwrap()]);
        let read = b"CCCCAAAAAAAAAAAACCCC";
        let hits = set.map(read, 0.25);
        assert_eq!(hits.len(), 1, "expected one binding site, got {hits:?}");
        assert_eq!(hits[0].mismatches, 0);
    }

    #[test]
    fn a_primer_longer_than_the_read_matches_nothing() {
        assert!(set().map(b"ACGT", 0.5).is_empty());
        assert!(scan(b"ACGT", b"", 0).is_empty());
        assert!(scan(b"", b"ACGT", 0).is_empty());
    }

    #[test]
    fn a_read_that_is_nothing_but_primer_is_left_alone_rather_than_emptied() {
        // The one hit is both the leading and the trailing candidate, so the
        // two cuts would cross and leave nothing at all.
        let read = rc(HCO).into_bytes();
        let plan =
            plan_trim(&set(), &read, &TrimOptions { search_window: 0, ..Default::default() });
        assert_eq!(plan.range, 0..read.len());
        assert!(!plan.trims_anything());
    }

    #[test]
    fn one_binding_site_is_never_used_as_both_cuts() {
        // A read that is a primer plus a couple of bases: the primer comes off
        // the front, and nothing pretends to be a 3' primer as well.
        let read = format!("{LCO}AC").into_bytes();
        let plan =
            plan_trim(&set(), &read, &TrimOptions { search_window: 0, ..Default::default() });
        assert_eq!(plan.range, LCO.len()..read.len());
        assert!(plan.end_hit.is_none());
    }

    #[test]
    fn a_primer_set_can_be_edited_by_name() {
        let mut s = set();
        assert_eq!(s.find_by_name("lco1490"), Some(0));
        assert_eq!(s.find_by_name("nope"), None);
        assert_eq!(s.remove(0).unwrap().name, "LCO1490");
        assert_eq!(s.len(), 1);
        assert!(s.remove(9).is_none());
    }
}
