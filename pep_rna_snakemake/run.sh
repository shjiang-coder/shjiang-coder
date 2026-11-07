#!/bin/bash
# ============================================================
# run.sh — 一键运行封装
# 用法:
#   bash run.sh          # 全流程
#   bash run.sh dry      # 干跑（只预览）
#   bash run.sh peptide  # 只跑小肽
#   bash run.sh rna      # 只跑小RNA
#   bash run.sh cores 48 # 指定核数
# ============================================================
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH=/pub/mambaforge/bin:/home/ug2451/.local/bin:$PATH

MODE="${1:-all}"
CORES="${2:-32}"
CONDA_FLAG=""
[ "${USE_CONDA:-1}" = "1" ] && CONDA_FLAG="--use-conda"

echo "============================================================"
echo " 小肽 + 小RNA 一体化流程"
echo " 模式: $MODE | 核数: $CORES | conda: ${USE_CONDA:-1}"
echo " 配置: config/config.yaml"
echo "============================================================"

case "$MODE" in
  dry)
      snakemake -n --cores "$CORES" $CONDA_FLAG
      ;;
  peptide)
      snakemake --cores "$CORES" $CONDA_FLAG --until peptide_summary_total --rerun-incomplete
      ;;
  rna)
      snakemake --cores "$CORES" $CONDA_FLAG \
          "$(python3 -c "
import yaml,sys
c=yaml.safe_load(open('config/config.yaml'))
print(c['outdir'].rstrip('/')+'/smallrna/sequences/_rna_summary.csv')
")" --rerun-incomplete
      ;;
  all|*)
      snakemake --cores "$CORES" $CONDA_FLAG --rerun-incomplete --keep-going
      ;;
esac

echo ""
echo "============================================================"
echo " ✅ 完成"
OUTDIR=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['outdir'])" 2>/dev/null || echo "<outdir>")
echo " 小肽汇总:    $OUTDIR/peptide/csv/_peptide_summary_total.csv"
echo " 小RNA汇总:   $OUTDIR/smallrna/sequences/_rna_summary.csv"
echo "============================================================"
