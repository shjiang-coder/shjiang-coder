#!/usr/bin/env python3
"""
rna_extract.py — 按 cmscan 坐标从基因组提取小RNA序列
输出: sequences/{en}/all_types.csv  (gene_id,sequence) + 按类型拆分 {type}.csv
"""
import csv, os, subprocess, shlex

CLS     = snakemake.input.cls
GENOME  = snakemake.input.genome
OUT     = snakemake.output.seqs
OUTDIR  = snakemake.params.outdir
CNAME   = snakemake.params.cname
MINLEN  = snakemake.params.minlen
THREADS = snakemake.threads

os.makedirs(OUTDIR, exist_ok=True)

# ── 读分类结果 ──
rows = []
with open(CLS, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try:
            rows.append({"gene_id": r["gene_id"], "type": r["rna_type"],
                         "seq_id": r["seq_id"], "start": int(r["start"]),
                         "end": int(r["end"]), "strand": r["strand"]})
        except (ValueError, KeyError):
            continue
print(f"[rna_extract] {CNAME}: {len(rows)} 条待提取", flush=True)

if not rows:
    for d in (OUT,):
        with open(d, 'w', encoding='utf-8') as f:
            f.write("gene_id,sequence\n")
    raise SystemExit(0)

# ── 建 faidx 索引 ──
if not os.path.exists(GENOME + ".fai"):
    subprocess.run(f"samtools faidx {shlex.quote(GENOME)}", shell=True, check=True)

# ── 批量提取 ──
regions = []
for r in rows:
    a, b = r["start"], r["end"]
    if a > b: a, b = b, a
    regions.append(f"{r['seq_id']}:{a}-{b}")

reg_file = os.path.join(OUTDIR, ".regions.txt")
with open(reg_file, 'w') as f:
    f.write('\n'.join(regions) + '\n')

proc = subprocess.run(
    f"samtools faidx -r {shlex.quote(reg_file)} {shlex.quote(GENOME)}",
    shell=True, capture_output=True, text=True)

seqs, cur, buf = {}, None, []
for line in proc.stdout.split('\n'):
    if line.startswith('>'):
        if cur: seqs[cur] = ''.join(buf)
        cur, buf = line[1:].split()[0], []
    elif line.strip():
        buf.append(line.strip())
if cur: seqs[cur] = ''.join(buf)

COMP = str.maketrans("ACTGNactgn", "TGACNtgacn")

# ── 输出 ──
by_type, n_ok = {}, 0
with open(OUT, 'w', newline='', encoding='utf-8') as f_out:
    w = csv.writer(f_out)
    w.writerow(["gene_id", "sequence"])
    for r, reg in zip(rows, regions):
        s = seqs.get(reg, "")
        if r["strand"] == "-":
            s = s.translate(COMP)[::-1]
        s = s.upper()
        if len(s) < MINLEN:
            continue
        w.writerow([r["gene_id"], s])
        by_type.setdefault(r["type"], []).append((r["gene_id"], s))
        n_ok += 1

for t, items in by_type.items():
    with open(os.path.join(OUTDIR, f"{t}.csv"), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(["gene_id", "sequence"])
        w.writerows(items)

os.remove(reg_file)
print(f"[rna_extract] {CNAME}: 成功 {n_ok}/{len(rows)} 条 → {OUT}", flush=True)
for t in sorted(by_type):
    print(f"    {t:14s} {len(by_type[t]):>6}", flush=True)
