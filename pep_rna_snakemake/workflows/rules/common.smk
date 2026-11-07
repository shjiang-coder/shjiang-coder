# ============================================================
# common.smk — 公共函数与规则（输出文件清单、目录创建、汇总）
# ============================================================

PEP_CSV  = os.path.join(OUTPE, "csv")
PEP_ORF  = os.path.join(OUTPE, "orfs")
PEP_NC   = os.path.join(OUTPE, "noncoding")
PEP_LOG  = os.path.join(OUTPE, "logs")

RNA_PRED = os.path.join(OUTRN, "predict")
RNA_CLS  = os.path.join(OUTRN, "classify")
RNA_SEQ  = os.path.join(OUTRN, "sequences")
RNA_LOG  = os.path.join(OUTRN, "logs")


def sp_iter():
    """遍历物种：(en_key, 中文名, genome路径)"""
    for cname, info in SPECIES.items():
        yield info["en"], cname, info["genome"]


# 英文名 → 物种信息（供规则里 lambda wildcards 使用）
SPECIES_BY_EN = {
    info["en"]: {"cname": cname, "en": info["en"], "genome": info["genome"]}
    for cname, info in SPECIES.items()
}


def peptide_outputs():
    """小肽模块的最终输出"""
    outs = [os.path.join(PEP_CSV, "_peptide_summary_total.csv")]
    return outs


def smallrna_outputs():
    """小RNA模块的最终输出"""
    outs = [os.path.join(RNA_SEQ, "_rna_summary.csv")]
    for en, cname, _ in sp_iter():
        outs.append(os.path.join(RNA_SEQ, en, "all_types.csv"))
    return outs


rule create_dirs:
    """一次性建好所有输出子目录"""
    output:
        touch(PEP_CSV + "/.dirs_ok"),
        touch(RNA_SEQ + "/.dirs_ok"),
    run:
        for d in [PEP_CSV, PEP_ORF, PEP_NC, PEP_LOG, RNA_PRED, RNA_CLS, RNA_SEQ, RNA_LOG]:
            os.makedirs(d, exist_ok=True)
