#!/usr/bin/env python3
"""
parse_getorf.py — getorf FASTA → CSV
输出列: peptide_id, sequence, length_aa, source
（由 Snakemake script: 调用，自动注入 snakemake 对象）
"""
import csv, re

FA = snakemake.input.orf
OUT = snakemake.output.csv
SRC = snakemake.params.source

recs, h, buf = [], None, []
with open(FA, errors='ignore') as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if h is not None:
                recs.append((h, ''.join(buf)))
            h, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
if h is not None:
    recs.append((h, ''.join(buf)))

with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(["peptide_id", "sequence", "length_aa", "source"])
    n = 0
    for hdr, s in recs:
        s = s.upper().replace('*', '')
        if not s:
            continue
        pid = hdr.split()[0]
        w.writerow([pid, s, len(s), SRC])
        n += 1

print(f"[parse_getorf] {SRC}: {n}/{len(recs)} 条 → {OUT}", flush=True)
