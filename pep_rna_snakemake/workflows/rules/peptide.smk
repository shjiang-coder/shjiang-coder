# ============================================================
# peptide.smk — 小肽模块规则
#  genome → getorf(编码区/非编码区) → FASTA → CSV → CD-HIT(NR100/NR90) → CSV → 汇总
#  说明: 所有文件名统一用英文名 {en}（避免中文路径在部分工具下出问题）；
#        中文名通过 params.source 写入 CSV 的 source 列。
# ============================================================

PE = PEP
_C  = PE.get("cdhit", {})
NR100 = _C.get("nr100", {"c": 1.0, "n": 5, "s": 1.0, "aL": 1.0, "g": 1})
NR90  = _C.get("nr90",  {"c": 0.9, "n": 5, "s": 0.0, "aL": 0.0, "g": 1})
CDHIT_THREADS = _C.get("threads", 16)
MIN_AA = PE.get("length", {}).get("min_aa", 10)
MAX_AA = PE.get("length", {}).get("max_aa", 100)
CODING_FIND = PE.get("coding", {}).get("find", 1)
NC_FIND     = PE.get("noncoding", {}).get("find", 0)

ALL_EN = [e for e, _, _ in sp_iter()]


# ── 1) 编码区 getorf ──
rule peptide_getorf_coding:
    input:
        genome = lambda wc: SPECIES_BY_EN[wc.en]["genome"],
    output:
        orf  = os.path.join(PEP_ORF, "{en}.coding.orf.fa"),
        done = os.path.join(PEP_LOG, "{en}.coding.getorf.done"),
    params:
        minsize = MIN_AA * 3,
        maxsize = MAX_AA * 3,
        find = CODING_FIND,
    threads: 2
    log: os.path.join(PEP_LOG, "{en}.coding.getorf.log")
    conda: "../envs/peptide.yaml"
    shell:
        """
        mkdir -p {PEP_ORF} {PEP_LOG}
        getorf -sequence {input.genome} -outseq {output.orf} \
               -minsize {params.minsize} -maxsize {params.maxsize} \
               -find {params.find} -reverse T > {log} 2>&1
        touch {output.done}
        """


# ── 2) 编码区 FASTA → CSV ──
rule peptide_csv_coding:
    input:
        orf  = os.path.join(PEP_ORF, "{en}.coding.orf.fa"),
        done = os.path.join(PEP_LOG, "{en}.coding.getorf.done"),
    output:
        csv = os.path.join(PEP_CSV, "{en}_peptides.csv"),
    params:
        source = lambda wc: SPECIES_BY_EN[wc.en]["cname"],
    conda: "../envs/peptide.yaml"
    script: "../../scripts/parse_getorf.py"


# ── 3) 编码区 CD-HIT NR100 / NR90 ──
rule peptide_cdhit_coding:
    input:
        csv = os.path.join(PEP_CSV, "{en}_peptides.csv"),
    output:
        nr100 = os.path.join(PEP_CSV, "{en}_peptides_nr100.csv"),
        nr90  = os.path.join(PEP_CSV, "{en}_peptides_nr90.csv"),
    params:
        tmpdir = PEP_ORF,
        source = lambda wc: SPECIES_BY_EN[wc.en]["cname"],
        c100 = NR100["c"], n100 = NR100["n"], s100 = NR100["s"], aL100 = NR100["aL"], g100 = NR100["g"],
        c90  = NR90["c"],  n90  = NR90["n"],  s90  = NR90["s"],  aL90  = NR90["aL"],  g90  = NR90["g"],
        threads = CDHIT_THREADS,
        prefix = "{en}",
    threads: CDHIT_THREADS
    conda: "../envs/peptide.yaml"
    script: "../../scripts/cdhit_pipeline.py"


# ── 4) 非编码区（可选）──
if PEP_NC_EN:
    rule peptide_getorf_noncoding:
        input:
            genome = lambda wc: SPECIES_BY_EN[wc.en]["genome"],
        output:
            orf  = os.path.join(PEP_NC, "{en}.nc.orf.fa"),
            done = os.path.join(PEP_LOG, "{en}.nc.getorf.done"),
        params:
            minsize = MIN_AA * 3,
            maxsize = MAX_AA * 3,
            find = NC_FIND,
        threads: 2
        log: os.path.join(PEP_LOG, "{en}.nc.getorf.log")
        conda: "../envs/peptide.yaml"
        shell:
            """
            mkdir -p {PEP_NC} {PEP_LOG}
            getorf -sequence {input.genome} -outseq {output.orf} \
                   -minsize {params.minsize} -maxsize {params.maxsize} \
                   -find {params.find} -reverse T > {log} 2>&1
            touch {output.done}
            """

    rule peptide_noncoding_pipeline:
        input:
            orf  = os.path.join(PEP_NC, "{en}.nc.orf.fa"),
            done = os.path.join(PEP_LOG, "{en}.nc.getorf.done"),
        output:
            raw   = os.path.join(PEP_NC, "{en}_noncoding.csv"),
            nr100 = os.path.join(PEP_NC, "{en}_noncoding_nr100.csv"),
        params:
            tmpdir = PEP_NC,
            source = lambda wc: SPECIES_BY_EN[wc.en]["cname"],
            c100 = NR100["c"], n100 = NR100["n"], s100 = NR100["s"], aL100 = NR100["aL"], g100 = NR100["g"],
            threads = CDHIT_THREADS,
            prefix = "{en}",
        threads: CDHIT_THREADS
        conda: "../envs/peptide.yaml"
        script: "../../scripts/noncoding_pipeline.py"


# ── 5) 小肽汇总总表 ──
rule peptide_summary_total:
    input:
        coding = expand(os.path.join(PEP_CSV, "{en}_peptides.csv"), en=ALL_EN),
        nr100  = expand(os.path.join(PEP_CSV, "{en}_peptides_nr100.csv"), en=ALL_EN),
        nr90   = expand(os.path.join(PEP_CSV, "{en}_peptides_nr90.csv"), en=ALL_EN),
        noncoding = (expand(os.path.join(PEP_NC, "{en}_noncoding.csv"), en=ALL_EN)
                     if PEP_NC_EN else []),
        nc_nr100  = (expand(os.path.join(PEP_NC, "{en}_noncoding_nr100.csv"), en=ALL_EN)
                     if PEP_NC_EN else []),
    output:
        csv = os.path.join(PEP_CSV, "_peptide_summary_total.csv"),
    params:
        nc_enabled = PEP_NC_EN,
    conda: "../envs/peptide.yaml"
    script: "../../scripts/peptide_summary.py"
