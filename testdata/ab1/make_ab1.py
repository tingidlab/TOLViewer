"""Generate small but structurally faithful .ab1 test files.

Peaks are Gaussians on a fixed spacing with a little cross-channel bleed, which
is enough to exercise the parser, the display scaling and the trimming rules.
"""
import math, struct, sys

SPACING = 12
ORDER = b"GATC"   # FWO_: DATA9..DATA12 carry G, A, T, C in that order

def trace_for(seq, spacing=SPACING, noise=0.04):
    n = len(seq) * spacing + spacing
    chans = {b: [0.0] * n for b in b"GATC"}
    for i, base in enumerate(seq):
        centre = spacing // 2 + i * spacing
        strong = base if base in b"GATC" else None
        for b in b"GATC":
            amp = 0.0
            if strong is None:
                amp = 0.25          # an N: all four about equal
            elif b == strong:
                amp = 1.0
            for s in range(max(0, centre - spacing), min(n, centre + spacing)):
                d = (s - centre) / (spacing * 0.30)
                chans[b][s] += amp * math.exp(-0.5 * d * d)
    out = {}
    for b in b"GATC":
        vals = []
        for s, v in enumerate(chans[b]):
            # A deterministic ripple stands in for baseline noise.
            v += noise * (0.5 + 0.5 * math.sin(s * 0.7 + b))
            vals.append(max(0, min(32767, int(v * 1800))))
        out[b] = vals
    peaks = [spacing // 2 + i * spacing for i in range(len(seq))]
    return out, peaks

def entry(name, number, etype, esize, count, size, offset):
    return struct.pack(">4slhhllll", name, number, etype, esize, count, size, offset, 0)

def build(seq, quality, sample_name, model=b"3730xl", date=(2026, 3, 14)):
    chans, peaks = trace_for(seq)
    items = []   # (name, number, etype, esize, count, payload)
    for i, b in enumerate(ORDER):
        vals = chans[b]
        items.append((b"DATA", 9 + i, 4, 2, len(vals), struct.pack(">%dh" % len(vals), *vals)))
    items.append((b"FWO_", 1, 2, 1, 4, ORDER))
    items.append((b"PBAS", 1, 2, 1, len(seq), seq))
    items.append((b"PLOC", 1, 4, 2, len(peaks), struct.pack(">%dh" % len(peaks), *peaks)))
    items.append((b"PCON", 1, 1, 1, len(quality), bytes(quality)))
    items.append((b"SMPL", 1, 18, 1, len(sample_name) + 1,
                  bytes([len(sample_name)]) + sample_name))
    items.append((b"MODL", 1, 19, 1, len(model) + 1, model + b"\0"))
    items.append((b"RUND", 1, 10, 4, 1, struct.pack(">hBB", date[0], date[1], date[2])))
    items.sort(key=lambda it: (it[0], it[1]))

    header_len = 128
    dir_off = header_len
    dir_len = 28 * len(items)
    data_off = dir_off + dir_len
    dir_bytes, blob = b"", b""
    for name, number, etype, esize, count, payload in items:
        size = len(payload)
        if size <= 4:
            off = int.from_bytes(payload.ljust(4, b"\0"), "big")
            dir_bytes += entry(name, number, etype, esize, count, size, off)
        else:
            dir_bytes += entry(name, number, etype, esize, count, size, data_off + len(blob))
            blob += payload
            if len(blob) % 2:
                blob += b"\0"
    head = b"ABIF" + struct.pack(">h", 101)
    head += entry(b"tdir", 1, 1023, 28, len(items), dir_len, dir_off)
    head = head.ljust(header_len, b"\0")
    return head + dir_bytes + blob

def qual(seq, lead_bad=18, tail_bad=22):
    q = []
    for i in range(len(seq)):
        if i < lead_bad or i >= len(seq) - tail_bad:
            q.append(6 + (i * 7) % 9)
        elif seq[i:i+1] == b"N":
            q.append(3)
        else:
            q.append(48 + (i * 3) % 12)
    return q

if __name__ == "__main__":
    # A COI fragment with the messy leader and trailer a real read has, and one
    # ambiguous call in the good stretch.
    good = (b"TTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCATTAAGAATTTTAATTCGAGCAGAA"
            b"TTAGGACAACCAGGATCATTAATTGGAGATGATCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTA"
            b"TTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGATTAGTTCCTTTAATATT")
    lead = b"NNCTGNATCGNNTTANCGGTNA"
    tail = b"GGNTANCNGTTANCGNNATCGGNTANCGNTTANCGGATNC"
    seq = lead + good[:60] + b"N" + good[61:] + tail
    open("testdata/ab1/tingidae_COI_F.ab1", "wb").write(build(seq, qual(seq), b"TL-2213_COI_F"))

    # The same locus sequenced from the reverse primer: the reverse complement,
    # so a round trip through reverse_complement() has something to check.
    comp = {ord('A'): 'T', ord('C'): 'G', ord('G'): 'C', ord('T'): 'A', ord('N'): 'N'}
    rc = bytes(ord(comp[c]) for c in reversed(seq))
    open("testdata/ab1/tingidae_COI_R.ab1", "wb").write(build(rc, qual(rc), b"TL-2213_COI_R"))
    print("wrote", len(seq), "base reads")
