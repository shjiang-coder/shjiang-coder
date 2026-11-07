#!/usr/bin/env python3
"""
cdhit_pipeline.py — 编码区小肽 CD-HIT 去冗余
CSV → FASTA → cd-hit(NR100) → CSV
                → cd-hit(NR90)  → CSV
"""
import csv, os, subprocess

IN_CSV  = snakemake.input.csv
OUT_100 = snakemake.output.nr100
OUT_90  = snakemake.output.nr90
P       = snakemake.params
TMP     = P.tmpdir
SRC     = P.source
PRE     = P.prefix
THREADS = P.threads

os.makedirs(TMP, exist_ok=True)
fa = os.path.join(TMP, f"{PRE}.cdhit_in.fa")

# ── CSV → FASTA ──
n_in = 0
with open(IN_CSV, encoding='utf-8') as f, open(fa, 'w') as out:
    for r in csv.DictReader(f):
        pid = (r.get('peptide_id') or '').strip()
        seq = (r.get('sequence') or '').strip()
        if pid and seq:
            out.write(f">{pid}\n{seq}\n")
            n_in += 1
print(f"[cdhit] {SRC}: 输入 {n_in} 条", flush=True)

if n_in == 0:
    for o in (OUT_100, OUT_90):
        with open(o, 'w', newline='', encoding='utf-8') as f:
            csv.writer(f).writerow(["peptide_id","sequence","length_aa","source"])
    print(f"[cdhit] {SRC}: 空输入，写出空 CSV", flush=True)
    raise SystemExit(0)


def run_cdhit(out_fa, c, n, s, aL, g, tag):
    clstr = out_fa + ".clstr"
    cmd = (f"cd-hit -i {fa} -o {out_fa} -c {c} -n {n} -d 0 -M 0 -T {THREADS} "
           f"-s {s} -aL {aL} -g {g}")
    subprocess.run(cmd, shell=True, check=True)
    n_out = sum(1 for line in open(out_fa) if line.startswith('>'))
    print(f"[cdhit] {SRC} {tag}: {n_in} → {n_out} 条", flush=True)
    return out_fa


def fa_to_csv(fa_path, out_csv):
    n = 0
    with open(fa_path) as f, open(out_csv, 'w', newline='', encoding='utf-8') as out:
        w = csv.writer(out)
        w.writerow(["peptide_id", "sequence", "length_aa", "source"])
        hid, buf = None, []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if hid is not None:
                    s = ''.join(buf)
                    w.writerow([hid, s, len(s), SRC]); n += 1
                hid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
        if hid is not None:
            s = ''.join(buf)
            w.writerow([hid, s, len(s), SRC]); n += 1
    return n


fa100 = run_cdhit(os.path.join(TMP, f"{PRE}.nr100.fa"),
                  P.c100, P.n100, P.s100, P.aL100, P.g100, "NR100")
fa_to_csv(fa100, OUT_100)

fa90 = run_cdhit(os.path.join(TMP, f"{PRE}.nr90.fa"),
                 P.c90, P.n90, P.s90, P.aL90, P.g90, "NR90")
fa_to_csv(fa90, OUT_90)

print(f"[cdhit] {SRC}: ✅ 完成", flush=True)
