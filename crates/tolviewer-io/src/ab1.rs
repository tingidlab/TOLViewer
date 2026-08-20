//! Applied Biosystems trace files (`.ab1`), the raw output of Sanger
//! sequencing.
//!
//! An ABIF file is a tagged container: a fixed header points at a directory of
//! 28-byte entries, each naming a four-character tag, an element type and
//! either an offset to its data or — for anything four bytes or smaller — the
//! data inlined in the offset field itself. Everything TOLViewer needs lives in
//! a handful of those tags:
//!
//! | tag | what it holds |
//! |-----|---------------|
//! | `DATA` 9–12 | the four analysed trace channels |
//! | `FWO_` 1 | which base each of those channels belongs to |
//! | `PBAS` 1/2 | the base calls (2 is the edited set, when the instrument software made one) |
//! | `PLOC` 1/2 | the trace sample each call peaks at |
//! | `PCON` 1/2 | per-call quality |
//! | `SMPL` 1 | the sample name the operator typed |
//!
//! Everything else is skipped. Files from other vendors' basecallers vary in
//! which optional tags they write, so a missing tag is only an error when the
//! trace is unusable without it.

use std::ops::Range;
use std::path::Path;

use tolviewer_core::{Alphabet, Error, Result, Sequence};

/// The magic at the start of every ABIF file.
const MAGIC: &[u8; 4] = b"ABIF";
/// Size of one directory entry, and of the root entry embedded in the header.
const ENTRY_SIZE: usize = 28;
/// Offset of the root ("tdir") entry inside the header.
const ROOT_OFFSET: usize = 4 + 2;

/// A Sanger chromatogram: four signal channels plus the calls made from them.
///
/// The four channels are stored in the file's own order, which varies by
/// instrument; [`Trace::channel_bases`] says which base each one belongs to and
/// [`Trace::channel_for`] does the lookup. `channels` are all the same length,
/// and `calls`, `peaks` and `quality` are all the same length as each other.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trace {
    /// Sample name from the `SMPL` tag, or the file stem when it is absent.
    pub sample_name: String,
    /// The base each channel carries, in channel order (the `FWO_` tag).
    pub channel_bases: [u8; 4],
    /// The four analysed traces, all of the same length.
    pub channels: [Vec<u16>; 4],
    /// One base call per peak, uppercase ASCII, `N` where the caller gave up.
    pub calls: Vec<u8>,
    /// The sample index in `channels` each call peaks at.
    pub peaks: Vec<u32>,
    /// Per-call confidence (Phred-like), when the file carries `PCON`.
    pub quality: Option<Vec<u8>>,
    /// Instrument and run description, for the info panel.
    pub comment: String,
}

impl Default for Trace {
    fn default() -> Self {
        Trace {
            sample_name: String::new(),
            channel_bases: *b"GATC",
            channels: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            calls: Vec::new(),
            peaks: Vec::new(),
            quality: None,
            comment: String::new(),
        }
    }
}

impl Trace {
    /// Number of base calls.
    pub fn len(&self) -> usize {
        self.calls.len()
    }

    pub fn is_empty(&self) -> bool {
        self.calls.is_empty()
    }

    /// Length of each trace channel in samples.
    pub fn samples(&self) -> usize {
        self.channels[0].len()
    }

    /// Which channel carries `base`, comparing case-insensitively and treating
    /// U as T.
    pub fn channel_for(&self, base: u8) -> Option<usize> {
        let want = match base.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        };
        self.channel_bases.iter().position(|&b| b.to_ascii_uppercase() == want)
    }

    /// Signal for `base` at `sample`, or 0 outside the trace.
    pub fn signal(&self, base: u8, sample: usize) -> u16 {
        match self.channel_for(base) {
            Some(c) => self.channels[c].get(sample).copied().unwrap_or(0),
            None => 0,
        }
    }

    /// The tallest signal anywhere in the trace, for scaling the display. At
    /// least 1, so callers can divide by it.
    pub fn peak_signal(&self) -> u16 {
        self.channels.iter().flat_map(|c| c.iter().copied()).max().unwrap_or(0).max(1)
    }

    /// The calls as a [`Sequence`], carrying quality across when it is known.
    pub fn to_sequence(&self, id: &str) -> Sequence {
        let mut seq = Sequence::new(id, self.calls.clone());
        seq.quality = self.quality.clone();
        if !self.comment.is_empty() {
            seq.description = self.comment.clone();
        }
        seq
    }

    /// Reverse complement the whole chromatogram: the channels are reversed and
    /// relabelled, and the peak positions are mirrored so calls still line up
    /// with the signal underneath them.
    ///
    /// This is the operation you want after sequencing from the reverse primer.
    pub fn reverse_complement(&mut self) {
        for channel in &mut self.channels {
            channel.reverse();
        }
        for base in &mut self.channel_bases {
            *base = Alphabet::Dna.complement(*base);
        }
        let samples = self.samples();
        self.calls.reverse();
        for call in &mut self.calls {
            *call = Alphabet::Dna.complement(*call);
        }
        self.peaks.reverse();
        for peak in &mut self.peaks {
            // Mirror about the trace, saturating rather than wrapping on the
            // malformed files that place a peak past the end of the signal.
            *peak = (samples.saturating_sub(1) as u32).saturating_sub(*peak);
        }
        if let Some(q) = &mut self.quality {
            q.reverse();
        }
    }

    /// Replace the base at `index`. The trace is untouched, so the display
    /// still shows what the instrument saw under the new call.
    pub fn set_call(&mut self, index: usize, base: u8) -> Result<()> {
        let len = self.len();
        let call = self.calls.get_mut(index).ok_or_else(|| {
            Error::out_of_range(format!("no base {} in a trace of {len} calls", index + 1))
        })?;
        *call = base.to_ascii_uppercase();
        // The instrument's confidence described the old call, not this one.
        if let Some(q) = &mut self.quality {
            q[index] = 0;
        }
        Ok(())
    }

    /// Insert a call at `index`, positioned midway between its neighbours.
    /// Used when the operator can see a peak the basecaller missed.
    pub fn insert_call(&mut self, index: usize, base: u8) -> Result<()> {
        if index > self.len() {
            return Err(Error::out_of_range(format!(
                "cannot insert at base {} of {}",
                index + 1,
                self.len()
            )));
        }
        let before = index.checked_sub(1).and_then(|i| self.peaks.get(i).copied());
        let after = self.peaks.get(index).copied();
        let peak = match (before, after) {
            (Some(a), Some(b)) => a + (b - a) / 2,
            (Some(a), None) => a.saturating_add(self.mean_spacing()),
            (None, Some(b)) => b.saturating_sub(self.mean_spacing()),
            (None, None) => 0,
        };
        self.calls.insert(index, base.to_ascii_uppercase());
        self.peaks.insert(index, peak);
        if let Some(q) = &mut self.quality {
            q.insert(index, 0);
        }
        Ok(())
    }

    /// Delete the call at `index`. The trace itself is kept, so the peak is
    /// still visible — it just no longer has a base attached.
    pub fn remove_call(&mut self, index: usize) -> Result<u8> {
        if index >= self.len() {
            return Err(Error::out_of_range(format!(
                "no base {} in a trace of {} calls",
                index + 1,
                self.len()
            )));
        }
        self.peaks.remove(index);
        if let Some(q) = &mut self.quality {
            q.remove(index);
        }
        Ok(self.calls.remove(index))
    }

    /// Typical distance between peaks, in samples.
    ///
    /// Used to place an inserted call at the end of a read, and by a viewer to
    /// work out how much room a call has on screen. Falls back to 12 samples,
    /// a common value, on traces too short to measure.
    pub fn mean_peak_spacing(&self) -> f32 {
        self.mean_spacing() as f32
    }

    /// As [`Trace::mean_peak_spacing`], in whole samples.
    fn mean_spacing(&self) -> u32 {
        match (self.peaks.first(), self.peaks.last()) {
            (Some(&first), Some(&last)) if self.peaks.len() > 1 => {
                (last.saturating_sub(first)) / (self.peaks.len() as u32 - 1)
            }
            _ => 12,
        }
    }

    /// Keep only the calls in `calls`, cropping the signal to the span around
    /// them so the chromatogram still lines up.
    ///
    /// The kept signal reaches halfway to each discarded neighbour, so a
    /// trimmed trace still shows the shoulders of its end peaks.
    pub fn trim(&mut self, calls: Range<usize>) -> Result<()> {
        if calls.start > calls.end || calls.end > self.len() {
            return Err(Error::out_of_range(format!(
                "cannot keep bases {}..{} of a {}-base trace",
                calls.start,
                calls.end,
                self.len()
            )));
        }
        if calls.is_empty() {
            *self = Trace {
                sample_name: std::mem::take(&mut self.sample_name),
                channel_bases: self.channel_bases,
                comment: std::mem::take(&mut self.comment),
                quality: self.quality.as_ref().map(|_| Vec::new()),
                ..Trace::default()
            };
            return Ok(());
        }
        let half = self.mean_spacing() / 2;
        let first = self.peaks[calls.start].saturating_sub(half) as usize;
        let last =
            (self.peaks[calls.end - 1].saturating_add(half) as usize + 1).min(self.samples());
        let first = first.min(last);

        for channel in &mut self.channels {
            *channel = channel[first..last].to_vec();
        }
        self.calls = self.calls[calls.clone()].to_vec();
        self.peaks = self.peaks[calls.clone()]
            .iter()
            .map(|&p| (p as usize).saturating_sub(first) as u32)
            .collect();
        if let Some(q) = &mut self.quality {
            *q = q[calls].to_vec();
        }
        Ok(())
    }

    /// The range of calls left after trimming low-quality ends, using the
    /// sliding-window rule Phred and most assemblers use: walk in from each end
    /// until a window of `window` calls averages at least `min_mean`.
    ///
    /// The window then gets tightened, because on its own it cuts too little: a
    /// window straddling the boundary passes as soon as enough of it is over
    /// good sequence, leaving a tail of junk inside the kept range. So each
    /// edge is advanced past any individual call that is itself below
    /// `min_mean`. The result is the read a person would have trimmed by eye.
    ///
    /// Returns an empty range starting at 0 when no window anywhere is good
    /// enough, and the whole read when the file carries no quality at all.
    pub fn quality_trim_range(&self, window: usize, min_mean: f32) -> Range<usize> {
        let Some(q) = &self.quality else { return 0..self.len() };
        let window = window.max(1);
        if q.len() < window {
            return 0..0;
        }
        let good = |start: usize| -> bool {
            let sum: u32 = q[start..start + window].iter().map(|&s| s as u32).sum();
            sum as f32 / window as f32 >= min_mean
        };
        let last_start = q.len() - window;
        let Some(first) = (0..=last_start).find(|&i| good(i)) else { return 0..0 };
        // Walk back from the far end for the 3' cut, so a good patch after a
        // bad tail does not extend the read into the noise.
        let last = (first..=last_start).rev().find(|&i| good(i)).map_or(q.len(), |i| i + window);

        let mut start = first;
        while start < last && f32::from(q[start]) < min_mean {
            start += 1;
        }
        let mut end = last;
        while end > start && f32::from(q[end - 1]) < min_mean {
            end -= 1;
        }
        start..end
    }
}

/// Parse the bytes of an `.ab1` file. `name` is used as the sample name when
/// the file does not carry one.
pub fn parse(bytes: &[u8], name: &str) -> Result<Trace> {
    let dir = Directory::read(bytes)?;

    let mut trace = Trace {
        sample_name: dir
            .string(b"SMPL", 1)
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| name.to_string()),
        ..Trace::default()
    };

    // FWO_ is four ASCII bases naming DATA9..DATA12 in order. It is stored
    // inline as a char[4], so it comes back as a string.
    if let Some(order) = dir.string(b"FWO_", 1) {
        let order = order.trim();
        let bytes = order.as_bytes();
        if bytes.len() == 4 && bytes.iter().all(|b| b"ACGTacgt".contains(b)) {
            for (slot, &b) in trace.channel_bases.iter_mut().zip(bytes) {
                *slot = b.to_ascii_uppercase();
            }
        }
    }

    for (i, tag) in (9..=12).enumerate() {
        trace.channels[i] = dir
            .shorts(b"DATA", tag)
            .map(|v| v.into_iter().map(|s| s.max(0) as u16).collect())
            .unwrap_or_default();
    }
    if trace.channels.iter().all(|c| c.is_empty()) {
        return Err(Error::parse(
            "AB1",
            None,
            "the file has no analysed trace channels (DATA 9-12); it may be a \
             raw capture that was never basecalled",
        ));
    }
    // A truncated channel would put every peak lookup out of step; pad the
    // short ones rather than rejecting a file that is otherwise readable.
    let samples = trace.channels.iter().map(|c| c.len()).max().unwrap_or(0);
    for channel in &mut trace.channels {
        channel.resize(samples, 0);
    }

    // Tag 2 is the edited call set, written when the instrument software or an
    // operator revised tag 1. Prefer it, exactly as the vendor's viewer does.
    let calls = dir
        .bytes(b"PBAS", 2)
        .or_else(|| dir.bytes(b"PBAS", 1))
        .ok_or_else(|| Error::parse("AB1", None, "the file has no base calls (PBAS)"))?;
    trace.calls = calls.into_iter().map(normalise_call).collect();

    let peaks = dir
        .shorts(b"PLOC", 2)
        .or_else(|| dir.shorts(b"PLOC", 1))
        .ok_or_else(|| Error::parse("AB1", None, "the file has no peak positions (PLOC)"))?;
    trace.peaks = peaks.into_iter().map(|p| p.max(0) as u32).collect();

    trace.quality = dir.bytes(b"PCON", 2).or_else(|| dir.bytes(b"PCON", 1));

    let mut parts: Vec<String> = Vec::new();
    if let Some(model) = dir.string(b"MODL", 1) {
        parts.push(model.trim().to_string());
    }
    if let Some(date) = dir.date(b"RUND", 1) {
        parts.push(date);
    }
    trace.comment = parts.into_iter().filter(|p| !p.is_empty()).collect::<Vec<_>>().join(" ");

    reconcile(&mut trace)?;
    Ok(trace)
}

/// Read an `.ab1` file from disk.
pub fn read_file(path: &Path) -> Result<Trace> {
    let bytes = std::fs::read(path)?;
    let name = path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    parse(&bytes, &name)
}

/// Do the tag lengths agree? Files in the wild disagree by a base or two when
/// an operator edited the calls without the peak list being rewritten. Trim to
/// the shortest rather than refusing to open the read.
fn reconcile(trace: &mut Trace) -> Result<()> {
    let n = trace.calls.len().min(trace.peaks.len());
    if n == 0 {
        return Err(Error::parse(
            "AB1",
            None,
            "the file's base calls and peak positions do not overlap",
        ));
    }
    trace.calls.truncate(n);
    trace.peaks.truncate(n);
    if let Some(q) = &mut trace.quality {
        if q.len() < n {
            q.resize(n, 0);
        } else {
            q.truncate(n);
        }
    }
    Ok(())
}

/// Basecallers write `N` for no call, but also occasionally a space or a NUL.
fn normalise_call(c: u8) -> u8 {
    if c.is_ascii_alphabetic() {
        c.to_ascii_uppercase()
    } else {
        b'N'
    }
}

/// The parsed directory: every entry's tag, number and raw bytes.
struct Directory<'a> {
    bytes: &'a [u8],
    entries: Vec<Entry>,
}

#[derive(Clone, Copy)]
struct Entry {
    name: [u8; 4],
    number: i32,
    element_type: u16,
    element_size: u16,
    element_count: u32,
    data_size: u32,
    data_offset: u32,
    /// Where the four inline bytes live when `data_size <= 4`.
    inline_at: usize,
}

impl<'a> Directory<'a> {
    fn read(bytes: &'a [u8]) -> Result<Directory<'a>> {
        if bytes.len() < ROOT_OFFSET + ENTRY_SIZE {
            return Err(Error::parse("AB1", None, "the file is too short to be an ABIF container"));
        }
        if &bytes[..4] != MAGIC {
            return Err(Error::parse(
                "AB1",
                None,
                "the file does not start with the ABIF signature; it is not an .ab1 trace",
            ));
        }
        let root = Entry::read(bytes, ROOT_OFFSET)?;
        if &root.name != b"tdir" {
            return Err(Error::parse("AB1", None, "the ABIF header has no directory entry"));
        }
        let count = root.element_count as usize;
        // A corrupt count could otherwise ask for gigabytes of entries.
        let max = bytes.len() / ENTRY_SIZE + 1;
        if count > max {
            return Err(Error::parse(
                "AB1",
                None,
                format!("the directory claims {count} entries, more than the file can hold"),
            ));
        }
        let start = root.data_offset as usize;
        let mut entries = Vec::with_capacity(count);
        for i in 0..count {
            let at = start.checked_add(i * ENTRY_SIZE).ok_or_else(|| {
                Error::parse("AB1", None, "the directory offset overflows the file")
            })?;
            // Stop at the first entry that runs off the end: the rest of the
            // file is still usable, and truncated downloads are common.
            match Entry::read(bytes, at) {
                Ok(entry) => entries.push(entry),
                Err(_) => break,
            }
        }
        if entries.is_empty() {
            return Err(Error::parse("AB1", None, "the ABIF directory is empty or unreadable"));
        }
        Ok(Directory { bytes, entries })
    }

    fn find(&self, name: &[u8; 4], number: i32) -> Option<&Entry> {
        self.entries.iter().find(|e| &e.name == name && e.number == number)
    }

    /// The raw bytes of an entry's data, inline or out of line.
    fn data(&self, name: &[u8; 4], number: i32) -> Option<&'a [u8]> {
        let entry = self.find(name, number)?;
        let size = entry.data_size as usize;
        if size <= 4 {
            return self.bytes.get(entry.inline_at..entry.inline_at + size);
        }
        let start = entry.data_offset as usize;
        self.bytes.get(start..start.checked_add(size)?)
    }

    /// A byte array (element type 1).
    fn bytes(&self, name: &[u8; 4], number: i32) -> Option<Vec<u8>> {
        Some(self.data(name, number)?.to_vec())
    }

    /// A big-endian 16-bit array (element types 3 and 4).
    ///
    /// Returns `None` when the entry says its elements are some other width,
    /// so a mislabelled tag reads as absent rather than as noise.
    fn shorts(&self, name: &[u8; 4], number: i32) -> Option<Vec<i16>> {
        let entry = self.find(name, number)?;
        if entry.element_size != 2 {
            return None;
        }
        let data = self.data(name, number)?;
        Some(data.chunks_exact(2).map(|c| i16::from_be_bytes([c[0], c[1]])).collect())
    }

    /// A string: Pascal (type 18), C (type 19) or a plain char array.
    fn string(&self, name: &[u8; 4], number: i32) -> Option<String> {
        let entry = self.find(name, number)?;
        let data = self.data(name, number)?;
        let slice = match entry.element_type {
            // Pascal string: a length byte, then that many characters.
            18 => {
                let len = (*data.first()? as usize).min(data.len().saturating_sub(1));
                &data[1..1 + len]
            }
            // C string: NUL-terminated.
            19 => {
                let end = data.iter().position(|&b| b == 0).unwrap_or(data.len());
                &data[..end]
            }
            _ => data,
        };
        Some(String::from_utf8_lossy(slice).trim_end_matches('\0').to_string())
    }

    /// A date (element type 10): big-endian year, then month and day bytes.
    fn date(&self, name: &[u8; 4], number: i32) -> Option<String> {
        let data = self.data(name, number)?;
        if data.len() < 4 {
            return None;
        }
        let year = i16::from_be_bytes([data[0], data[1]]);
        Some(format!("{year:04}-{:02}-{:02}", data[2], data[3]))
    }
}

impl Entry {
    fn read(bytes: &[u8], at: usize) -> Result<Entry> {
        let raw = bytes
            .get(at..at + ENTRY_SIZE)
            .ok_or_else(|| Error::parse("AB1", None, "a directory entry runs past the end"))?;
        let be32 = |o: usize| i32::from_be_bytes([raw[o], raw[o + 1], raw[o + 2], raw[o + 3]]);
        let be16 = |o: usize| u16::from_be_bytes([raw[o], raw[o + 1]]);
        let element_count = be32(12);
        let data_size = be32(16);
        let data_offset = be32(20);
        if element_count < 0 || data_size < 0 || data_offset < 0 {
            return Err(Error::parse("AB1", None, "a directory entry has a negative length"));
        }
        Ok(Entry {
            name: [raw[0], raw[1], raw[2], raw[3]],
            number: be32(4),
            element_type: be16(8),
            element_size: be16(10),
            element_count: element_count as u32,
            data_size: data_size as u32,
            data_offset: data_offset as u32,
            inline_at: at + 20,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One directory entry to build: tag, number, element type, element size,
    /// payload.
    type Item<'a> = (&'a [u8; 4], i32, u16, u16, Vec<u8>);

    /// Assemble an ABIF container from `Item`s, so tests can build exactly the
    /// file they mean — including the broken ones.
    fn container(items: &[Item<'_>]) -> Vec<u8> {
        const HEADER: usize = 128;
        let dir_len = items.len() * ENTRY_SIZE;
        let data_start = HEADER + dir_len;
        let mut dir = Vec::new();
        let mut blob: Vec<u8> = Vec::new();
        for (name, number, etype, esize, payload) in items {
            let size = payload.len() as u32;
            let offset = if size <= 4 {
                let mut inline = [0u8; 4];
                inline[..payload.len()].copy_from_slice(payload);
                u32::from_be_bytes(inline)
            } else {
                let at = (data_start + blob.len()) as u32;
                blob.extend_from_slice(payload);
                at
            };
            let count = if *esize == 0 { size } else { size / *esize as u32 };
            dir.extend_from_slice(&entry(name, *number, *etype, *esize, count, size, offset));
        }
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&101i16.to_be_bytes());
        out.extend_from_slice(&entry(
            b"tdir",
            1,
            1023,
            ENTRY_SIZE as u16,
            items.len() as u32,
            dir_len as u32,
            HEADER as u32,
        ));
        out.resize(HEADER, 0);
        out.extend_from_slice(&dir);
        out.extend_from_slice(&blob);
        out
    }

    fn entry(
        name: &[u8; 4],
        number: i32,
        etype: u16,
        esize: u16,
        count: u32,
        size: u32,
        offset: u32,
    ) -> [u8; ENTRY_SIZE] {
        let mut e = [0u8; ENTRY_SIZE];
        e[0..4].copy_from_slice(name);
        e[4..8].copy_from_slice(&number.to_be_bytes());
        e[8..10].copy_from_slice(&etype.to_be_bytes());
        e[10..12].copy_from_slice(&esize.to_be_bytes());
        e[12..16].copy_from_slice(&count.to_be_bytes());
        e[16..20].copy_from_slice(&size.to_be_bytes());
        e[20..24].copy_from_slice(&offset.to_be_bytes());
        e
    }

    fn shorts(values: &[i16]) -> Vec<u8> {
        values.iter().flat_map(|v| v.to_be_bytes()).collect()
    }

    /// A four-base read: one clean Gaussian-ish peak per call, 8 samples apart.
    fn tiny() -> Vec<u8> {
        let calls = b"ACGT";
        let peaks: Vec<i16> = (0..4).map(|i| 4 + i * 8).collect();
        let channel = |base: u8| -> Vec<u8> {
            let mut samples = vec![0i16; 36];
            for (i, &c) in calls.iter().enumerate() {
                if c == base {
                    let centre = 4 + i * 8;
                    for (d, amp) in [(-1i32, 400i16), (0, 1000), (1, 400)] {
                        samples[(centre as i32 + d) as usize] = amp;
                    }
                }
            }
            shorts(&samples)
        };
        container(&[
            (b"DATA", 9, 4, 2, channel(b'G')),
            (b"DATA", 10, 4, 2, channel(b'A')),
            (b"DATA", 11, 4, 2, channel(b'T')),
            (b"DATA", 12, 4, 2, channel(b'C')),
            (b"FWO_", 1, 2, 1, b"GATC".to_vec()),
            (b"PBAS", 1, 2, 1, calls.to_vec()),
            (b"PLOC", 1, 4, 2, shorts(&peaks)),
            (b"PCON", 1, 1, 1, vec![40, 41, 42, 43]),
            (b"SMPL", 1, 18, 1, {
                let n = b"specimen 7";
                let mut v = vec![n.len() as u8];
                v.extend_from_slice(n);
                v
            }),
        ])
    }

    #[test]
    fn reads_calls_peaks_quality_and_name() {
        let t = parse(&tiny(), "fallback").unwrap();
        assert_eq!(t.calls, b"ACGT");
        assert_eq!(t.peaks, vec![4, 12, 20, 28]);
        assert_eq!(t.quality.as_deref(), Some(&[40u8, 41, 42, 43][..]));
        assert_eq!(t.sample_name, "specimen 7");
        assert_eq!(t.channel_bases, *b"GATC");
        assert_eq!(t.samples(), 36);
    }

    #[test]
    fn signal_is_looked_up_by_base_not_channel_order() {
        let t = parse(&tiny(), "x").unwrap();
        // 'A' is the second channel here, and its call peaks at sample 4.
        assert_eq!(t.channel_for(b'A'), Some(1));
        assert_eq!(t.signal(b'A', 4), 1000);
        assert_eq!(t.signal(b'C', 4), 0);
        // Lowercase and RNA spellings resolve to the same channel.
        assert_eq!(t.channel_for(b't'), t.channel_for(b'T'));
        assert_eq!(t.channel_for(b'U'), t.channel_for(b'T'));
        assert_eq!(t.channel_for(b'X'), None);
    }

    #[test]
    fn the_edited_call_set_wins_over_the_original() {
        let mut items: Vec<Item<'_>> = Vec::new();
        let base = tiny();
        let original = parse(&base, "x").unwrap();
        // Rebuild with both PBAS 1 and PBAS 2 present.
        let channel =
            |i: usize| shorts(&original.channels[i].iter().map(|&v| v as i16).collect::<Vec<_>>());
        for (n, i) in [(9, 0), (10, 1), (11, 2), (12, 3)] {
            items.push((b"DATA", n, 4, 2, channel(i)));
        }
        items.push((b"FWO_", 1, 2, 1, b"GATC".to_vec()));
        items.push((b"PBAS", 1, 2, 1, b"ACGT".to_vec()));
        items.push((b"PBAS", 2, 2, 1, b"ACNT".to_vec()));
        items.push((b"PLOC", 1, 4, 2, shorts(&[4, 12, 20, 28])));
        let t = parse(&container(&items), "x").unwrap();
        assert_eq!(t.calls, b"ACNT", "PBAS 2 is the operator's revision and must win");
    }

    #[test]
    fn a_non_abif_file_is_rejected_by_name() {
        let e = parse(b">seq\nACGT\n", "x").unwrap_err();
        match e {
            Error::Parse { message, .. } => assert!(message.contains("ABIF"), "{message}"),
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn truncation_and_nonsense_never_panic() {
        let full = tiny();
        for cut in 0..full.len() {
            let _ = parse(&full[..cut], "x");
        }
        let mut corrupt = full.clone();
        // A directory claiming a preposterous number of entries.
        corrupt[4 + 2 + 12..4 + 2 + 16].copy_from_slice(&0x7fff_ffffi32.to_be_bytes());
        assert!(parse(&corrupt, "x").is_err());
        // Offsets pointing off the end.
        let mut corrupt = full.clone();
        corrupt[128 + 20..128 + 24].copy_from_slice(&0xffff_fffeu32.to_be_bytes());
        let _ = parse(&corrupt, "x");
    }

    #[test]
    fn a_file_without_traces_says_so() {
        let items: Vec<Item<'_>> = vec![
            (b"PBAS", 1, 2, 1, b"ACGT".to_vec()),
            (b"PLOC", 1, 4, 2, shorts(&[4, 12, 20, 28])),
        ];
        let e = parse(&container(&items), "x").unwrap_err();
        match e {
            Error::Parse { message, .. } => assert!(message.contains("DATA 9-12"), "{message}"),
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn calls_longer_than_peaks_are_trimmed_not_rejected() {
        let items: Vec<Item<'_>> = vec![
            (b"DATA", 9, 4, 2, shorts(&[0; 36])),
            (b"DATA", 10, 4, 2, shorts(&[0; 36])),
            (b"DATA", 11, 4, 2, shorts(&[0; 36])),
            (b"DATA", 12, 4, 2, shorts(&[0; 36])),
            (b"PBAS", 1, 2, 1, b"ACGTACGT".to_vec()),
            (b"PLOC", 1, 4, 2, shorts(&[4, 12, 20])),
            (b"PCON", 1, 1, 1, vec![40; 8]),
        ];
        let t = parse(&container(&items), "x").unwrap();
        assert_eq!(t.calls, b"ACG");
        assert_eq!(t.peaks.len(), 3);
        assert_eq!(t.quality.as_ref().unwrap().len(), 3);
    }

    #[test]
    fn reverse_complement_keeps_calls_over_their_peaks() {
        let mut t = parse(&tiny(), "x").unwrap();
        let samples = t.samples();
        // Record what each call sits on top of, then check it still does.
        let before: Vec<(u8, u16)> =
            t.calls.iter().zip(&t.peaks).map(|(&c, &p)| (c, t.signal(c, p as usize))).collect();
        t.reverse_complement();
        assert_eq!(t.calls, b"ACGT", "ACGT is its own reverse complement");
        assert_eq!(t.samples(), samples);
        let after: Vec<(u8, u16)> =
            t.calls.iter().zip(&t.peaks).map(|(&c, &p)| (c, t.signal(c, p as usize))).collect();
        assert_eq!(before, after);
        assert_eq!(t.channel_bases, *b"CTAG", "the channels are relabelled, not reordered");
    }

    #[test]
    fn reverse_complement_twice_is_the_identity() {
        let original = parse(&tiny(), "x").unwrap();
        let mut t = original.clone();
        t.reverse_complement();
        t.reverse_complement();
        assert_eq!(t, original);
    }

    #[test]
    fn editing_a_call_clears_the_instruments_confidence_in_it() {
        let mut t = parse(&tiny(), "x").unwrap();
        t.set_call(2, b'a').unwrap();
        assert_eq!(t.calls, b"ACAT", "the call is stored uppercase");
        assert_eq!(t.quality.as_ref().unwrap()[2], 0);
        assert!(t.set_call(9, b'A').is_err());
    }

    #[test]
    fn an_inserted_call_lands_between_its_neighbours() {
        let mut t = parse(&tiny(), "x").unwrap();
        t.insert_call(2, b'n').unwrap();
        assert_eq!(t.calls, b"ACNGT");
        assert!(t.peaks[1] < t.peaks[2] && t.peaks[2] < t.peaks[3]);
        assert_eq!(t.quality.as_ref().unwrap().len(), 5);
        // At either end there is only one neighbour to work from.
        t.insert_call(0, b'A').unwrap();
        t.insert_call(6, b'A').unwrap();
        assert_eq!(t.len(), 7);
        assert!(t.insert_call(99, b'A').is_err());
    }

    #[test]
    fn removing_a_call_leaves_the_trace_alone() {
        let mut t = parse(&tiny(), "x").unwrap();
        let samples = t.samples();
        assert_eq!(t.remove_call(1).unwrap(), b'C');
        assert_eq!(t.calls, b"AGT");
        assert_eq!(t.samples(), samples, "the signal is evidence; only the call goes");
        assert!(t.remove_call(3).is_err());
    }

    #[test]
    fn trimming_crops_the_signal_with_the_calls() {
        let mut t = parse(&tiny(), "x").unwrap();
        let (call, peak) = (t.calls[1], t.peaks[1]);
        let signal_under_it = t.signal(call, peak as usize);
        t.trim(1..3).unwrap();
        assert_eq!(t.calls, b"CG");
        assert!(t.samples() < 36);
        assert_eq!(
            t.signal(t.calls[0], t.peaks[0] as usize),
            signal_under_it,
            "a trimmed call must still sit on its own peak"
        );
        assert_eq!(t.quality.as_ref().unwrap().len(), 2);
    }

    #[test]
    fn trimming_everything_away_is_legal_and_empty() {
        let mut t = parse(&tiny(), "x").unwrap();
        t.trim(2..2).unwrap();
        assert!(t.is_empty());
        assert_eq!(t.samples(), 0);
        assert_eq!(t.sample_name, "specimen 7", "the identity of the read survives");
        assert!(t.quality.is_some());
        let mut t = parse(&tiny(), "x").unwrap();
        assert!(t.trim(1..99).is_err());
        // Built by hand: a reversed range literal is a clippy error.
        assert!(t.trim(Range { start: 3, end: 1 }).is_err());
    }

    #[test]
    fn quality_trimming_finds_the_good_middle() {
        let mut t = parse(&tiny(), "x").unwrap();
        t.calls = vec![b'A'; 40];
        t.peaks = (0..40).map(|i| i * 8).collect();
        let mut q = vec![5u8; 40];
        q[10..30].fill(45);
        t.quality = Some(q);
        assert_eq!(t.quality_trim_range(5, 20.0), 10..30, "the edges must land on the good run");
    }

    #[test]
    fn quality_trimming_of_a_hopeless_read_keeps_nothing() {
        let mut t = parse(&tiny(), "x").unwrap();
        t.quality = Some(vec![2, 2, 2, 2]);
        assert!(t.quality_trim_range(3, 30.0).is_empty());
        // Shorter than the window: nothing can pass.
        assert!(t.quality_trim_range(50, 1.0).is_empty());
        // No quality at all means no basis for trimming, so keep everything.
        t.quality = None;
        assert_eq!(t.quality_trim_range(5, 30.0), 0..4);
    }

    #[test]
    fn the_checked_in_traces_read_back() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/ab1");
        let f = read_file(&dir.join("tingidae_COI_F.ab1")).unwrap();
        assert_eq!(f.sample_name, "TL-2213_COI_F");
        assert!(f.len() > 250, "{} calls", f.len());
        assert_eq!(f.len(), f.peaks.len());
        assert!(f.comment.contains("3730xl") && f.comment.contains("2026-03-14"), "{}", f.comment);

        // Every call should be the tallest channel at its own peak, except the
        // Ns, which are deliberately flat.
        let mut agree = 0;
        for (&call, &peak) in f.calls.iter().zip(&f.peaks) {
            if call == b'N' {
                continue;
            }
            let mine = f.signal(call, peak as usize);
            if b"ACGT".iter().all(|&b| b == call || f.signal(b, peak as usize) <= mine) {
                agree += 1;
            }
        }
        assert_eq!(agree, f.calls.iter().filter(|&&c| c != b'N').count());

        // The reverse read is the same fragment from the other primer.
        let mut r = read_file(&dir.join("tingidae_COI_R.ab1")).unwrap();
        r.reverse_complement();
        assert_eq!(r.calls, f.calls);
    }

    #[test]
    fn quality_trimming_the_real_read_cuts_both_messy_ends() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/ab1");
        let f = read_file(&dir.join("tingidae_COI_F.ab1")).unwrap();
        // The generator makes the first 18 and last 22 calls poor.
        let range = f.quality_trim_range(20, 20.0);
        assert_eq!(range, 18..f.len() - 22, "the messy ends must go, and nothing more");
    }
}
