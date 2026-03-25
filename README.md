# genome_circos_plot

这是一个用于绘制基因组圈图（circos 风格）的 Python 脚本，支持 GC 含量、基因密度、重复序列密度、BGC、EggNOG/COG 注释、表达谱热图、可变剪切等多轨道数据的可视化，并支持通过配置文件精细控制布局与图例位置。

> 适用场景：细菌/真菌基因组多条 contig/scaffold/chromosome 的整体展示与多组学信息叠加。

---

## 功能特性

- 自动计算：GC 含量、基因密度、重复序列密度
- 支持可选数据叠加：
  - EggNOG/COG 注释
  - 表达谱热图（多条件）
  - BGC 预测结果
  - 可变剪切事件
- 轨道布局与图例位置完全可配置（`track_layout.conf.json`）
- 支持自动扩展外圈（表达轨道数增加时不压缩内圈）
- 输出 PNG 图像 + Excel 汇总表

---

## 依赖

- Python 3.8+
- `matplotlib`
- `pandas`

建议安装：
```bash
pip install matplotlib pandas
```

---

## 快速开始

### 最小运行（仅 FASTA + GFF3）
```bash
python3 genome_circos.py \
  --fasta genome_assembly.fasta \
  --gff3 annotation.gff3 \
  --output test_circos.png \
  --excel test_circos.xlsx
```

### 完整运行（带表达/剪切/COG/BGC）
```bash
python3 genome_circos.py \
  --fasta genome_assembly.fasta \
  --gff3 annotation.gff3 \
  --eggnog eggnog.tsv \
  --expr expression.tsv \
  --bgc bgc.tsv \
  --splice splicing.tsv \
  --conf track_layout.conf.json \
  --output test_circos.png \
  --excel test_circos.xlsx
```

### 查看帮助
```bash
python3 genome_circos.py -h
```

---

## 输入文件格式

### 1) 基因组序列（FASTA）
- `--fasta genome_assembly.fasta`
- 脚本会自动识别序列名中是否包含：
  - `chromosome` → 显示为 `C1, C2...`
  - `scaffold` → 显示为 `S1, S2...`
  - `contig` → 显示为 `c1, c2...`

### 2) 注释文件（GFF3）
- `--gff3 annotation.gff3`
- 需要包含基因坐标（gene/CDS）

### 3) EggNOG/COG 注释（可选）
- `--eggnog eggnog.tsv`
- 默认列名：
  - `query`（基因 ID）
  - `COG_category`（COG 类别）

### 4) 表达谱（可选）
- `--expr expression.tsv`
- 默认列名：
  - `gene_id` / `condition` / `log2FC` / `padj`

### 5) BGC 结果（可选）
- `--bgc bgc.tsv`
- 默认列名：
  - `contig` / `start` / `end` / `type` / `strand`

### 6) 可变剪切（可选）
- `--splice splicing.tsv`
- 默认列名：
  - `seqid` / `start` / `end` / `event`

> 如列名不同，可通过命令行参数自定义列名。

---

## 配置文件（track_layout.conf.json）

该文件控制：
- 各轨道的半径位置（r0/r1）
- 图注位置/大小
- 表达色条位置
- 染色体标签的样式
- 自动扩展外圈参数

已包含大量中文注释，直接编辑即可。

常用参数示例：
```json
"auto_expr_splice": true,
"expression_track_width": 0.04,
"splicing_width": 0.06,
"expr_gap": 0.02
```

```json
"image": {
  "figsize": [15, 15],
  "dpi": 300,
  "auto_size": true
}
```

---

## 输出文件

- 图像（PNG）：`--output`
- 汇总表（Excel）：`--excel`

Excel 内包含：
- GC 含量
- 基因密度
- 重复序列密度
- BGC 列表
- EggNOG/COG 注释
- 表达谱
- 可变剪切事件

---

## 常见问题

**Q: 表达轨道数量变多后图像会重叠怎么办？**
A: 已支持自动扩展外圈和自动增大画布大小，请确保 `auto_expr_splice = true`。

**Q: COG 没有注释的区域如何显示？**
A: 默认灰色填充，可在 conf 中调整 `empty_color`。

---

## 作者说明

本脚本为定制化可视化工具，适合基因组多轨道综合展示。欢迎根据需要扩展。
