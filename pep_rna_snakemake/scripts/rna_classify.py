#!/usr/bin/env python3
"""
rna_classify.py — cmscan tblout → 8 类分类
输出: classify/{en}/all_types.csv   列: gene_id,rna_type,rfam_id,seq_id,start,end,strand,score,evalue
     同时按类型拆成 8 个文件: classify/{en}/{type}.csv（gene_id,rna_type,rfam_id,pos）
"""
import csv, os, re

TBLOUT = snakemake.input.tblout
OUT    = snakemake.output.cls
LOOKUP = snakemake.params.lookup
TYPES  = snakemake.params.types.split(",")
CNAME  = snakemake.params.cname
OUTDIR = snakemake.params.outdir

os.makedirs(OUTDIR, exist_ok=True)

# ── Rfam accession → 类型 ──
lookup = {}
with open(LOOKUP) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            lookup[p[0]] = p[2]

def rna_type(rfam_acc, family):
    t = lookup.get(rfam_acc, "")
    if t in TYPES:
        return t
    # 兜底启发式
    fl = (family or "").lower()
    if rfam_acc == "RF00005" or "trna" in fl:      return "tRNA"
    if "mir" in fl:                                return "miRNA"
    if "rrna" in fl or fl.startswith(("5s","5.8s","18s","28s","ssu","lsu")): return "rRNA"
    if "snorna" in fl or fl.startswith("snor"):    return "snoRNA"
    if re.match(r'^u\d', fl) or "snrna" in fl:     return "snRNA"
    if "srp" in fl:                                return "SRP"
    if "ribozyme" in fl:                           return "ribozyme"
    return "other_ncRNA"

# ── 解析 tblout（cmscan --tblout 标准格式，无 idx 列）──
# 列: target(家族) accession query(scaffold) acc mdl mdlfrom mdlto seqfrom seqto strand trunc pass gc bias score evalue inc desc
rows = []
with open(TBLOUT, errors='ignore') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) < 16:
            continue
        family, rf_acc, seq_id = p[0], p[1].split('/')[0], p[2]
        try:
            sf, st, strand = int(p[7]), int(p[8]), p[9]
            score, evalue = float(p[14]), float(p[15])
        except (ValueError, IndexError):
            continue
        t = rna_type(rf_acc, family)
        rows.append({"gene_id": f"{CNAME}|{t}|{rf_acc}|{seq_id}:{sf}-{st}",
                     "rna_type": t, "rfam_id": rf_acc, "family": family,
                     "seq_id": seq_id, "start": sf, "end": st, "strand": strand,
                     "score": round(score,2), "evalue": f"{evalue:.2e}"})

# ── 写总表 ──
with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=["gene_id","rna_type","rfam_id","family",
                                      "seq_id","start","end","strand","score","evalue"])
    w.writeheader(); w.writerows(rows)

# ── 按类型拆分（每个类型一个文件，无表头，格式: gene_id,type,rfam,pos）──
by_type = {t: [] for t in TYPES}
for r in rows:
    by_type.setdefault(r["rna_type"], []).append(r)
for t, items in by_type.items():
    with open(os.path.join(OUTDIR, f"{t}.csv"), 'w', encoding='utf-8') as f:
        for r in items:
            f.write(f"{r['gene_id']},{r['rna_type']},{r['rfam_id']},{r['seq_id']}:{r['start']}-{r['end']}\n")

print(f"[rna_classify] {CNAME}: {len(rows)} hits", flush=True)
for t in TYPES:
    print(f"    {t:14s} {len(by_type.get(t,[])):>6}", flush=True)
