# 小肽 + 小RNA 一体化 Snakemake 自动化流程


把之前散装的「小肽预测（getorf → CD-HIT → CSV）」和「小RNA 注释（cmscan → 8类 → 序列）」两套流程，合并为**一个 Snakemake 工作流**，一条命令跑完，直接输出 CSV 表格。

---

## 一、目录结构

```
pep_rna_snakemake/
├── Snakefile                    主入口（含 include 与 all 规则）
├── config/
│   └── config.yaml              ★★★ 唯一需要您编辑的文件 ★★★
├── workflows/
│   ├── rules/
│   │   ├── common.smk           公共变量/函数（输出路径、物种映射）
│   │   ├── peptide.smk          小肽规则（5 条）
│   │   └── smallrna.smk         小RNA规则（4 条，enabled=false 时整块跳过）
│   └── envs/
│       ├── peptide.yaml         conda 环境：emboss(getorf) + cd-hit + samtools
│       └── rna.yaml             conda 环境：infernal(cmscan) + samtools
├── scripts/
│   ├── parse_getorf.py          getorf FASTA → CSV
│   ├── cdhit_pipeline.py        编码区 CD-HIT（NR100 + NR90）
│   ├── noncoding_pipeline.py    非编码区 CD-HIT（NR100）
│   ├── peptide_summary.py       小肽汇总总表
│   ├── rna_classify.py          cmscan tblout → 8 类
│   ├── rna_extract.py           按坐标提取小RNA序列
│   └── rna_summary.py           小RNA汇总总表
└── run.sh                       一键运行封装（含干跑/正式/断点续跑）
```

---

## 二、⭐ config.yaml 怎么编辑（重点）

打开 `config/config.yaml`，只需改 **4 个地方**：

### ① 输出目录（第 1 项）

```yaml
outdir: /home/ug2451/shjiang/project/pep_rna_results
```
所有结果都写到这里，**不存在会自动创建**。想换地方就改这一行。

### ② 物种列表（第 2 项）★最常改★

```yaml
species:
  佛手:                                    # ← 中文名（写进 CSV 的 source 列）
    en: foshou                             # ← 英文名（用于文件名，建议纯英文小写）
    genome: /绝对路径/xxx.fna              # ← 基因组或 CDS 的完整路径

  昆明山海棠:
    en: kunshanhaitang
    genome: /绝对路径/yyy.fna
```

**编辑要点：**
- **键（佛手）= 中文名**，会写进 CSV 的 `source` 列和汇总表的「中文名」列
- **`en` = 英文名**，只用于文件名（`foshou_peptides.csv`），建议纯小写英文
- **`genome` = 绝对路径**，可给基因组 FASTA，也可给 CDS FASTA（都可）
- **加物种**：照格式往下追加一段，缩进用 **2 个空格**（不能用 Tab）
- **减物种**：删掉或整段注释（在行首加 `#`）
- 想临时只跑一个物种做测试：把其他物种注释掉即可

### ③ 小肽模块（第 3 项）

```yaml
peptide:
  enabled: true              # ← 不想跑小肽就改 false
  length:
    min_aa: 10               # 肽段最短长度（氨基酸）
    max_aa: 100              # 肽段最长长度
  coding:                    # 编码区（ATG 起始）
    enabled: true
    find: 1                  # 1 = 只取 ATG 起始的 ORF
  noncoding:                 # 非编码区（stop-to-stop）
    enabled: true            # ← 不想跑非编码区改 false（能省一半时间）
    find: 0
  cdhit:
    threads: 16              # CD-HIT 线程
    nr100: {c: 1.0, n: 5, s: 1.0, aL: 1.0, g: 1}    # 精确去冗余
    nr90:  {c: 0.9, n: 5, s: 0.0, aL: 0.0, g: 1}    # 90% 聚类
```

**常见调整：**
| 想做什么 | 改哪里 |
|:---------|:-------|
| 只要 5-50 aa 的肽 | `min_aa: 5`, `max_aa: 50` |
| 不要非编码区（省时间） | `noncoding.enabled: false` |
| 要 NR80 | 加一段 `nr80: {c: 0.8, n: 5, ...}`（需同步改 peptide.smk，见「进阶」）|
| 只做精确去重不要聚类 | 把 `nr90` 里的 `c` 改成 `1.0` |
| CD-HIT 太慢 | `cdhit.threads` 调大（≤ CPU 核数）|

### ④ 小RNA模块（第 4 项）

```yaml
smallrna:
  enabled: true              # ← 不想跑小RNA就改 false
  rfam_cm: /path/Rfam.cm     # Rfam 协方差模型库（必须）
  rfam_lookup: /path/rfam_type_lookup.tsv   # Rfam ID→类型映射表（必须）
  cmscan:
    cpu: 16                  # 单物种线程（机器共享时调小，如 8）
    rfam_mode: true          # 目标库是 Rfam 时必须 true
    cut_ga: true             # 用 Rfam 官方 GA 阈值（推荐）
    extra: ""                # 附加参数，如 "-E 1e-3"
  rna_types:                 # 8 类，一般不用改
    - tRNA
    - miRNA
    - rRNA
    - snoRNA
    - snRNA
    - SRP
    - ribozyme
    - other_ncRNA
  min_seq_len: 10            # 提取序列最短长度
```

**⚠️ 小RNA 是最慢的一步**（全基因组 cmscan，单物种可能数小时到数天）。
- 只想先要小肽结果：`smallrna.enabled: false`
- 想加速：`cmscan.cpu` 调大 + 减少物种数；或先用 CDS 代替全基因组

---

## 三、运行方法

```bash
cd /home/ug2451/shjiang/project/pep_rna_snakemake

# ① 先干跑，确认要做什么（强烈建议）
snakemake -n --cores 32

# ② 正式运行（用系统已装工具）
snakemake --cores 32

# ③ 正式运行（推荐：自动创建 conda 环境，版本可复现）
snakemake --cores 32 --use-conda

# ④ 只跑小肽模块
snakemake --cores 32 --until peptide_summary_total

# ⑤ 中断后续跑（Snakemake 自动跳过已完成步骤）
snakemake --cores 32 --rerun-incomplete

# ⑥ 只重跑某个物种（删掉它的中间文件即可）
rm -rf <outdir>/peptide/csv/foshou_*
```

**或用封装脚本：**
```bash
bash run.sh dry      # 干跑
bash run.sh          # 全流程
bash run.sh peptide  # 只跑小肽
bash run.sh rna      # 只跑小RNA
```

---

## 四、输出文件说明

### 小肽（`{outdir}/peptide/`）

| 文件 | 内容 |
|:-----|:-----|
| `csv/{en}_peptides.csv` | 编码区原始预测 |
| `csv/{en}_peptides_nr100.csv` | 编码区 NR100 去冗余 |
| `csv/{en}_peptides_nr90.csv` | 编码区 NR90 聚类 |
| `noncoding/{en}_noncoding.csv` | 非编码区原始 |
| `noncoding/{en}_noncoding_nr100.csv` | 非编码区 NR100 |
| **`csv/_peptide_summary_total.csv`** | ★ **汇总总表** |

**小肽 CSV 格式（4 列）：**
```csv
peptide_id,sequence,length_aa,source
NC_079511.1_cds_XP_018447753.1_1_1,MKTTTITMKDFYGTRNSTVRILNREVRIKKEDQETIK...,60,佛手
```

**汇总总表格式：**
```csv
中文名,英文名,编码区原始,编码区NR100,编码区NR90,非编码区原始,非编码区NR100,合计NR100
佛手,foshou,5411454,4022327,3409000,,,4022327
```

### 小RNA（`{outdir}/smallrna/`）

| 文件 | 内容 |
|:-----|:-----|
| `predict/{en}.tblout` | cmscan 原始表格输出 |
| `classify/{en}/all_types.csv` | 全部命中 + 类型（10 列）|
| `classify/{en}/{type}.csv` | 按 8 类拆分 |
| **`sequences/{en}/{type}.csv`** | ★ **各类序列**（gene_id,sequence）|
| `sequences/{en}/all_types.csv` | 该物种全部序列 |
| **`_rna_summary.csv`** | ★ **汇总总表** |

**序列 CSV 格式：**
```csv
gene_id,sequence
佛手|tRNA|RF00005|MKYO02000002.1:160919-160846,GCAGATGTAGCCAAGTGGATTAAGGCAATGGATT...
```

**汇总总表格式：**
```csv
中文名,英文名,tRNA,miRNA,rRNA,snoRNA,snRNA,SRP,ribozyme,other_ncRNA,合计
佛手,foshou,849,168,584,25,96,14,2,1712,3450
```

---

## 五、验证状态

**已通过实测（2026-11-07）：**
- ✅ Snakemake 干跑：DAG 建成 27 个任务（3 物种）
- ✅ 小肽模块真实执行：测试基因组 496 KB
  - 编码区 4,544 → NR100 2,930 → NR90 2,632
  - 非编码区 19,232 → NR100 12,378
  - 汇总表输出正确

**⚠️ 小RNA 模块需实跑验证**（全基因组 cmscan 耗时较长，建议先用小基因组/单物种试跑）

---

## 六、注意事项

1. **`--configfile` 会与 `config/config.yaml` 合并**，不是替换。若要临时改配置，建议直接编辑 `config/config.yaml`，或用 `--config key=value` 覆盖单项。
2. **文件名统一用英文名 `{en}`**（避免中文工具兼容问题）；中文名写入 CSV 的 `source` / `中文名` 列。
3. **基因组需可读**；`rna_extract` 会自动 `samtools faidx` 建索引。
4. **断点续跑**：Snakemake 靠输出文件判断完成，中断后重跑自动跳过；若某步失败需删对应输出再跑。
5. **cmscan 完整性检查**：流程会校验 tblout 末尾是否有 `# [ok]`，没有则报错（防止截断结果被误用）。
6. 想加 **NR80/NR70** 等版本：在 `peptide.smk` 里照 `nr90` 的写法加一条 `run_cdhit` 调用即可。

---

## 七、环境依赖

| 模块 | 依赖 | 系统已装？ |
|:-----|:-----|:----------:|
| Snakemake | snaKemake ≥7 | ✅ 7.32.3 |
| 小肽 | emboss(getorf)、cd-hit、samtools、python3 | ✅ |
| 小RNA | infernal(cmscan)、samtools、python3 | ✅ |

用 `--use-conda` 可自动按 `workflows/envs/*.yaml` 创建隔离环境（推荐用于论文复现）。
