"""E15 data acquisition: fetch only the neocortex MRF entries from BrainSpan per-donor zip archives.

Source: http://download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10/mrf/<DONOR>.zip
(GRCh37, GENCODE v10; MRF = RSEQtools Mapped Read Format: alignment blocks only, no sequence).

Method:
- Parse each remote zip's central directory through HTTP Range requests (zipfile on a
  range-reading file object, so ZIP64 is handled).
- For every entry whose region is in NEOCORTEX, download the exact byte span of the local header plus
  the compressed data. Each span is saved as <entry>.zipentry.
- Record name, CRC-32, compressed/uncompressed sizes and offsets in manifest.tsv. The CRC is checked
  when the entry is decompressed (E15 junction extraction), and the byte length is checked here.

No analysis is done here. Selection is by file name only (region code), so no outcome is involved.
Resumable: entries already present with the correct length are skipped.

Usage: python -I E15_fetch_brainspan_mrf.py <out_dir>
"""
import io
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

BASE = "http://download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10/mrf/"
# BrainSpan neocortical areas: 11 adult-defined areas plus the early-prenatal composite or precursor
# areas used before the 11 areas are dissectable. Codes as they appear in the MRF entry names
# (checked on HSB112, 8 pcw: M1C-S1C, Ocx, PCx; TCx is included by analogy).
NEOCORTEX = {"A1C", "DFC", "IPC", "ITC", "M1C", "MFC", "OFC", "S1C", "STC", "V1C", "VFC",
             "M1C-S1C", "Ocx", "PCx", "TCx"}
DONORS = ("HSB103 HSB105 HSB107 HSB112 HSB113 HSB114 HSB118 HSB119 HSB121 HSB122 HSB123 HSB124 HSB126 "
          "HSB127 HSB130 HSB132 HSB135 HSB136 HSB139 HSB141 HSB143 HSB145 HSB148 HSB149 HSB150 HSB153 "
          "HSB154 HSB155 HSB159 HSB168 HSB169 HSB171 HSB172 HSB173 HSB174 HSB175 HSB178 HSB194 HSB92 "
          "HSB96 HSB97 HSB98").split()


def get_range(url, a, b, tries=30):
    """Fetch bytes a..b inclusive; on a dropped connection keep what arrived and request only the rest."""
    buf = bytearray()
    for t in range(tries):
        start = a + len(buf)
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{b}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                if r.status != 206:
                    raise IOError(f"expected 206, got {r.status}")
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    buf += chunk
        except Exception as e:  # network instability seen earlier in this project
            print(f"  retry {t + 1} {start}-{b} ({len(buf)} B kept): {e}", flush=True)
        if len(buf) == b - a + 1:
            return bytes(buf)
        if len(buf) > b - a + 1:
            raise IOError(f"range {a}-{b}: received more bytes than requested")
        time.sleep(min(60, 5 * (t + 1)))
    raise IOError(f"range {a}-{b} failed for {url}")


class RemoteFile(io.RawIOBase):
    def __init__(self, url):
        self.url, self.pos = url, 0
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=120) as r:
            self.size = int(r.headers["Content-Length"])

    def seekable(self):
        return True

    def readable(self):
        return True

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def tell(self):
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        b[:n] = get_range(self.url, self.pos, self.pos + n - 1)
        self.pos += n
        return n


def region_of(name):
    # e.g. HSB_103_A1C_R_RNASeq.mrf -> A1C
    parts = name.split("_")
    return parts[2] if len(parts) >= 5 else None


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    man = out / "manifest.tsv"
    if not man.exists():
        man.write_text("donor\tentry\tregion\tcrc32\tcompress_size\tfile_size\tspan_start\tspan_len\n")
    done = {l.split("\t")[1] for l in man.read_text().splitlines()[1:]}
    for d in DONORS:
        url = BASE + d + ".zip"
        rf = RemoteFile(url)
        zf = zipfile.ZipFile(io.BufferedReader(rf, buffer_size=1 << 16))
        infos = sorted(zf.infolist(), key=lambda i: i.header_offset)
        ends = [i.header_offset for i in infos[1:]] + [rf.size]
        for info, end in zip(infos, ends):
            reg = region_of(info.filename)
            if reg not in NEOCORTEX or info.filename in done:
                continue
            # local header: 30 bytes + name + extra (read to get the exact data start)
            lh = get_range(url, info.header_offset, info.header_offset + 29)
            n, e = int.from_bytes(lh[26:28], "little"), int.from_bytes(lh[28:30], "little")
            span = 30 + n + e + info.compress_size
            assert info.header_offset + span <= end, (info.filename, span, end)
            dest = out / (info.filename + ".zipentry")
            tmp = dest.with_suffix(".part")
            CH = 64 << 20
            with tmp.open("wb") as fh:
                for a in range(info.header_offset, info.header_offset + span, CH):
                    b = min(a + CH, info.header_offset + span) - 1
                    fh.write(get_range(url, a, b))
            assert tmp.stat().st_size == span
            tmp.replace(dest)
            with man.open("a") as fh:
                fh.write(f"{d}\t{info.filename}\t{reg}\t{info.CRC:08x}\t{info.compress_size}\t{info.file_size}"
                         f"\t{info.header_offset}\t{span}\n")
            print(f"{d} {info.filename} {span / 1e6:.0f} MB", flush=True)
        print(f"{d} done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
