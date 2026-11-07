# ============================================================
# smallrna.smk — 小RNA模块规则
#  genome → cmscan(Rfam) → tblout → 8类分类 → 序列提取 → 汇总
#  说明: 当 config 中 smallrna.enabled=false 时，整块规则不定义（流程自动跳过）
# ============================================================

if RNA_EN:   # ← 关键：关闭时整块跳过，不会读取缺失的配置项

    _CMS = RNA.get("cmscan", {})
    CMS_CPU   = _CMS.get("cpu", 16)
    CMS_RFAM  = _CMS.get("rfam_mode", True)
    CMS_GA    = _CMS.get("cut_ga", True)
    CMS_EXTRA = _CMS.get("extra", "")
    MIN_SEQ_LEN = RNA.get("min_seq_len", 10)
    RFAM_CM     = RNA.get("rfam_cm", "")
    RFAM_LOOKUP = RNA.get("rfam_lookup", "")

    # ── 1) cmscan ──
    rule rna_cmscan:
        input:
            genome = lambda wc: SPECIES_BY_EN[wc.en]["genome"],
        output:
            tblout = os.path.join(RNA_PRED, "{en}.tblout"),
            cmout  = os.path.join(RNA_PRED, "{en}.cmscan"),
            done   = os.path.join(RNA_LOG, "{en}.cmscan.done"),
        params:
            rfam = "--rfam" if CMS_RFAM else "",
            ga   = "--cut_ga" if CMS_GA else "",
            extra = CMS_EXTRA,
            cm   = RFAM_CM,
        threads: CMS_CPU
        log: os.path.join(RNA_LOG, "{en}.cmscan.log")
        conda: "../envs/rna.yaml"
        shell:
            """
            mkdir -p {RNA_PRED} {RNA_LOG}
            cmscan --cpu {threads} {params.rfam} {params.ga} {params.extra} \
                   --tblout {output.tblout} -o {output.cmout} \
                   {params.cm} {input.genome} > {log} 2>&1
            grep -q '^\\# \\[ok\\]' {output.tblout} || {{ echo "cmscan 未正常结束，请查看 {log}"; exit 1; }}
            touch {output.done}
            """

    # ── 2) 8 类分类 ──
    rule rna_classify:
        input:
            tblout = os.path.join(RNA_PRED, "{en}.tblout"),
            done   = os.path.join(RNA_LOG, "{en}.cmscan.done"),
        output:
            cls = os.path.join(RNA_CLS, "{en}", "all_types.csv"),
        params:
            lookup = RFAM_LOOKUP,
            types  = ",".join(RNA_TYPES),
            cname  = lambda wc: SPECIES_BY_EN[wc.en]["cname"],
            outdir = os.path.join(RNA_CLS, "{en}"),
        conda: "../envs/rna.yaml"
        script: "../../scripts/rna_classify.py"

    # ── 3) 按坐标提取序列 ──
    rule rna_extract:
        input:
            cls    = os.path.join(RNA_CLS, "{en}", "all_types.csv"),
            genome = lambda wc: SPECIES_BY_EN[wc.en]["genome"],
        output:
            seqs = os.path.join(RNA_SEQ, "{en}", "all_types.csv"),
        params:
            cname  = lambda wc: SPECIES_BY_EN[wc.en]["cname"],
            minlen = MIN_SEQ_LEN,
            outdir = os.path.join(RNA_SEQ, "{en}"),
        threads: 4
        conda: "../envs/rna.yaml"
        script: "../../scripts/rna_extract.py"

    # ── 4) 小RNA汇总总表 ──
    rule rna_summary:
        input:
            seqs = expand(os.path.join(RNA_SEQ, "{en}", "all_types.csv"), en=[e for e, _, _ in sp_iter()]),
        output:
            csv = os.path.join(RNA_SEQ, "_rna_summary.csv"),
        params:
            types = ",".join(RNA_TYPES),
        conda: "../envs/rna.yaml"
        script: "../../scripts/rna_summary.py"
