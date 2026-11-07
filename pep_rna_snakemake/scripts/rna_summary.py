#!/usr/bin/env python3
"""
rna_summary.py — 小RNA 汇总总表
输出: _rna_summary.csv  列: 中文名,英文名,tRNA,miRNA,rRNA,snoRNA,snRNA,SRP,ribozyme,other_ncRNA,合计
"""
import csv, os

OUT   = snakemake.output.csv
TYPES = snakemake.params.types.split(",")
SPEC  = snakemake.config["species"]
SEQD  = os.path.join(snakemake.config["outdir"], "smallrna", "sequences")


def cnt(path):
    if not os.path.exists(path):
        return 0
    with open(path, encoding='utf-8') as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)   # 去表头


rows, totals = [], {t: 0 for t in TYPES}
for cname, info in SPEC.items():
    en = info["en"]
    vals = [cnt(os.path.join(SEQD, en, f"{t}.csv")) for t in TYPES]
    rows.append([cname, en] + vals + [sum(vals)])
    for t, v in zip(TYPES, vals):
        totals[t] += v
rows.append(["合计", "-"] + [totals[t] for t in TYPES] + [sum(totals.values())])

with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(["中文名", "英文名"] + TYPES + ["合计"])
    w.writerows(rows)

print(f"[rna_summary] → {OUT}", flush=True)
hdr = ["物种"] + TYPES + ["合计"]
print("  " + "  ".join(f"{h:<12}" for h in hdr), flush=True)
for r in rows:
    print("  " + "  ".join(f"{str(x):<12}" for x in r), flush=True)
