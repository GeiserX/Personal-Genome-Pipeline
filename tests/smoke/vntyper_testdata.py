#!/usr/bin/env python3
"""Fetch members of VNtyper's public test data without the whole archive.

VNtyper publishes its test BAMs as one 1.2 GB zip on Zenodo (record 19181821,
"VNtyper 2.0 test data", the archive_file of tests/test_data_config.json in
VNtyper v2.0.47). A zip keeps its file list at the end, so this script reads
that list with an HTTP range request, then fetches and inflates only the
members it is asked for, about 10 MB instead of 1.2 GB (the .bai files are
not fetched: samtools index writes them). Each member is checked
against the CRC-32 the archive stores and against the sha256 pinned below.

Usage: vntyper_testdata.py OUTDIR NAME...   (NAME: a key of MEMBERS)
Standard library only: it runs in VNTYPER_IMAGE's Python and on the host.
"""
import hashlib
import os
import struct
import sys
import time
import urllib.request
import zlib

URL = "https://zenodo.org/records/19181821/files/data.zip?download=1"
SIZE = 1174975271  # bytes; Zenodo lists md5 fdd5958d58c733cdca452f3287898964 for the whole file

# name -> (path in the archive, sha256 of the inflated file)
MEMBERS = {
    # Positive: a MUC1 VNTR frameshift VNtyper calls High_Precision* (its
    # integration test example_dfc3_GRCh38_fast runs the same reads on NCBI
    # contig names). Reads realigned by VNtyper to GRCh38 chr1, UCSC names.
    "positive.bam": ("data/remapped/bwa/hg38/example_dfc3_hg38_bwa.bam",
                     "5934330a8ac2144ac0cc6ae8c84d7d9262e90b28127417e4068f7d662b7c64bb"),
    # Negative: no variant called (its integration test
    # example_40cf_hg38_subset_fast_gdp_guard), mean VNTR depth about 148x.
    "negative.bam": ("data/example_40cf_hg38_subset.bam",
                     "44aff8e3b2096b7ea3529989660052c61b29c8a26afd7f48de4988f7047d0652"),
}


def get(start, end):
    """Bytes start..end (inclusive) of the archive, with retries."""
    for attempt in range(4):
        try:
            req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                if r.status != 206:
                    raise OSError(f"HTTP {r.status}, want 206 (a range)")
                data = r.read()
            if len(data) != end - start + 1:
                raise OSError(f"got {len(data)} bytes, want {end - start + 1}")
            return data
        except OSError as e:
            if attempt == 3:
                raise
            print(f"range {start}-{end}: {e}; again in {15 * (attempt + 1)} s", file=sys.stderr)
            time.sleep(15 * (attempt + 1))
    raise AssertionError("unreachable")


def directory():
    """{path: (method, crc, compressed size, size, local header offset)}"""
    tail = get(SIZE - 65536, SIZE - 1)
    i = tail.rfind(b"PK\x05\x06")
    if i < 0:
        raise SystemExit("no end-of-central-directory record in the last 64 kB")
    _, _, _, _, n, cd_size, cd_off, _ = struct.unpack("<IHHHHIIH", tail[i:i + 22])
    cd = get(cd_off, cd_off + cd_size - 1)
    out, p = {}, 0
    for _ in range(n):
        f = struct.unpack("<IHHHHHHIIIHHHHHII", cd[p:p + 46])
        if f[0] != 0x02014B50:
            raise SystemExit(f"bad central directory entry at {p}")
        nl, el, cl = f[10], f[11], f[12]
        out[cd[p + 46:p + 46 + nl].decode()] = (f[4], f[7], f[8], f[9], f[16])
        p += 46 + nl + el + cl
    return out


def main(argv):
    if len(argv) < 2 or any(n not in MEMBERS for n in argv[1:]):
        raise SystemExit(f"usage: {os.path.basename(__file__)} OUTDIR NAME... (NAME of {', '.join(MEMBERS)})")
    outdir = argv[0]
    os.makedirs(outdir, exist_ok=True)
    cd = directory()
    for name in argv[1:]:
        path, want = MEMBERS[name]
        method, crc, csize, size, off = cd[path]
        head = get(off, off + 29)
        nl, el = struct.unpack("<HH", head[26:30])
        raw = get(off + 30 + nl + el, off + 30 + nl + el + csize - 1)
        data = zlib.decompressobj(-15).decompress(raw) if method == 8 else raw
        if len(data) != size or zlib.crc32(data) != crc:
            raise SystemExit(f"{path}: size or CRC-32 differs from the archive's directory")
        got = hashlib.sha256(data).hexdigest()
        if got != want:
            raise SystemExit(f"{path}: sha256 {got}, pinned {want}")
        with open(os.path.join(outdir, name), "wb") as f:
            f.write(data)
        print(f"{name}\t{path}\t{size}\t{got}")


if __name__ == "__main__":
    main(sys.argv[1:])
