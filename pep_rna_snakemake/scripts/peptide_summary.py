#!/usr/bin/env python3
"""
peptide_summary.py — 小肽汇总总表
输出列: 中文名,英文名,编码区原始,编码区NR100,编码区NR90,非编码区原始,非编码区NR100,合计NR100
"""
import csv, os

OUT     = snakemake.output.csv
NC_EN   = snakemake.params.nc_enabled
SPECIES = snakemake.config["species"]

PEP_CSV = os.path.join(snakemake.config["outdir"], "peptide", "csv")
PEP_NC  = os.path.join(snakemake.config["outdir"], "peptide", "noncoding")


def cnt(path):
    """CSV 数据行数（不含表头）"""
    if not os.path.exists(path):
        return 0
    with open(path, encoding='utf-8') as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)


rows = []
g = {k: 0 for k in ["raw", "n100", "n90", "nc", "nc100", "total100"]}
for cname, info in SPECIES.items():
    en = info["en"]
    raw   = cnt(os.path.join(PEP_CSV, f"{en}_peptides.csv"))
    n100  = cnt(os.path.join(PEP_CSV, f"{en}_peptides_nr100.csv"))
    n90   = cnt(os.path.join(PEP_CSV, f"{en}_peptides_nr90.csv"))
    nc    = cnt(os.path.join(PEP_NC,  f"{en}_noncoding.csv")) if NC_EN else 0
    nc100 = cnt(os.path.join(PEP_NC,  f"{en}_noncoding_nr100.csv")) if NC_EN else 0
    total100 = n100 + nc100
    rows.append([cname, en, raw, n100, n90, nc, nc100, total100])
    for k, v in zip(["raw","n100","n90","nc","nc100","total100"], [raw,n100,n90,nc,nc100,total100]):
        g[k] += v

rows.append(["合计", "-", g["raw"], g["n100"], g["n90"], g["nc"], g["nc100"], g["total100"]])

with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(["中文名","英文名","编码区原始","编码区NR100","编码区NR90",
                "非编码区原始","非编码区NR100","合计NR100"])
    w.writerows(rows)

print(f"[peptide_summary] → {OUT}")
for r in rows:
    print("  " + "  ".join(str(x) for x in r), flush=True)
