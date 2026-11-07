#!/usr/bin/env python3
"""
noncoding_pipeline.py — 非编码区小肽：FASTA → CSV → CD-HIT NR100 → CSV
"""
import csv, os, subprocess

FA      = snakemake.input.orf
OUT_RAW = snakemake.output.raw
OUT_100 = snakemake.output.nr100
P       = snakemake.params
TMP     = P.tmpdir
SRC     = P.source
PRE     = P.prefix
THREADS = P.threads

os.makedirs(TMP, exist_ok=True)

# ── 解析 getorf FASTA ──
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

fa = os.path.join(TMP, f"{PRE}.nc_in.fa")
with open(fa, 'w') as f, open(OUT_RAW, 'w', newline='', encoding='utf-8') as out:
    w = csv.writer(out)
    w.writerow(["peptide_id", "sequence", "length_aa", "source"])
    n = 0
    for hdr, s in recs:
        s = s.upper().replace('*', '')
        if not s:
            continue
        pid = hdr.split()[0]
        f.write(f">{pid}\n{s}\n")
        w.writerow([pid, s, len(s), SRC])
        n += 1
print(f"[noncoding] {SRC}: 原始 {n} 条 → {OUT_RAW}", flush=True)

if n == 0:
    with open(OUT_100, 'w', newline='', encoding='utf-8') as f:
        csv.writer(f).writerow(["peptide_id","sequence","length_aa","source"])
    raise SystemExit(0)

# ── CD-HIT NR100 ──
out_fa = os.path.join(TMP, f"{PRE}.nc_nr100.fa")
cmd = (f"cd-hit -i {fa} -o {out_fa} -c {P.c100} -n {P.n100} -d 0 -M 0 -T {THREADS} "
       f"-s {P.s100} -aL {P.aL100} -g {P.g100}")
subprocess.run(cmd, shell=True, check=True)

m = 0
with open(out_fa) as f, open(OUT_100, 'w', newline='', encoding='utf-8') as out:
    w = csv.writer(out)
    w.writerow(["peptide_id", "sequence", "length_aa", "source"])
    hid, buf = None, []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if hid is not None:
                s = ''.join(buf); w.writerow([hid, s, len(s), SRC]); m += 1
            hid, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if hid is not None:
        s = ''.join(buf); w.writerow([hid, s, len(s), SRC]); m += 1
print(f"[noncoding] {SRC}: NR100 {n} → {m} 条", flush=True)
