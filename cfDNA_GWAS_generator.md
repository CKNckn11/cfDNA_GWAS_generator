# cfdna-gwas-generate 使用教程

`cfdna-gwas-generate` 是一个 cfDNA germline GWAS / germline VCF 流程生成器。它不会直接跑分析，而是根据 `config.yaml` 和样本表生成一个独立项目目录，里面包含可运行的 shell 脚本。

当前命令帮助：

```bash
cfdna-gwas-generate -h
```

输出：

```text
usage: cfdna-gwas-generate [-h] --config CONFIG [--outdir OUTDIR]

Generate a cfDNA germline GWAS shell pipeline.

options:
  -h, --help       show this help message and exit
  --config CONFIG  YAML config file.
  --outdir OUTDIR  Override project.outdir.
```

## 1. 激活环境

如果在 `biosoft` 环境：

```bash
conda activate biosoft
```

如果在独立 `cfDNA` 环境：

```bash
conda activate cfDNA
```

确认命令可用：

```bash
which cfdna-gwas-generate
cfdna-gwas-generate -h
```

如果 `biosoft` 环境里找不到命令，可以用完整路径：

```bash
/home/ckn/.local/bin/cfdna-gwas-generate -h
```

## 2. 输入类型

通过 `config.yaml` 里的 `input.type` 区分初始文件类型。

### BAM 输入

```yaml
input:
  type: bam
  sample_table: /path/to/samples.bam.tsv
```

BAM 样本表至少需要：

```text
sample_id	bam
```

推荐格式：

```text
sample_id	bam	bai	phenotype	age	sex
Lib-1005	/path/Lib-1005.bam	/path/Lib-1005.bam.bai	Breast Cancer	50	1
```

说明：

- `sample_id`：样本 ID，也可以用 `iid`
- `bam`：BAM 文件路径
- `bai`：BAM 索引路径
- `phenotype`：表型，也可以用 `diagnosis_type`
- `age`：年龄
- `sex`：性别，建议 `1=male, 2=female`

### FASTQ 输入

```yaml
input:
  type: fastq
  sample_table: /path/to/samples.fastq.tsv
```

FASTQ 样本表至少需要：

```text
sample_id	fastq1	fastq2
```

推荐格式：

```text
sample_id	fastq1	fastq2	phenotype	age	sex
Lib-001	/path/Lib-001_R1.fastq.gz	/path/Lib-001_R2.fastq.gz	Cancer	60	1
```

## 3. 配置文件示例

BAM 示例配置在：

```text
/400T/ckn/cfDNA_GWAS_generator/examples/config.bam.example.yaml
```

FASTQ 示例配置在：

```text
/400T/ckn/cfDNA_GWAS_generator/examples/config.fastq.example.yaml
```

Breast 项目当前配置在：

```text
/400T/ckn/CFDNA_Breast_GWAS/config.breast.bam.yaml
```

Breast 样本表在：

```text
/400T/ckn/CFDNA_Breast_GWAS/samples.breast.bam.tsv
```

## 4. 生成项目

进入任意目录都可以运行，例如：

```bash
cfdna-gwas-generate \
  --config /400T/ckn/CFDNA_Breast_GWAS/config.breast.bam.yaml
```

如果想临时覆盖 `config.yaml` 里的输出目录：

```bash
cfdna-gwas-generate \
  --config /400T/ckn/CFDNA_Breast_GWAS/config.breast.bam.yaml \
  --outdir /400T/ckn/CFDNA_Breast_GWAS_test
```

生成后项目目录里会有：

```text
00.shell/
00_metadata/
01_fastq_qc/
02_alignment/
03_markdup/
04_bqsr/
05_bam_qc/
06_gvcf/
07_joint_vcf/
08_plink_qc/
09_pca/
10_gwas/
11_plots/
loginfo/
config.sh
config.yaml
README.md
```

BAM 输入时，`01_fastq_qc`、`02_alignment`、`03_markdup`、`04_bqsr` 目录通常不会使用。

## 5. 运行流程

进入生成后的项目目录：

```bash
cd /400T/ckn/CFDNA_Breast_GWAS
```

查看可用步骤：

```bash
bash 00.shell/run_pipeline.sh help
```

BAM 输入流程：

```text
metadata -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot
```

FASTQ 输入流程：

```text
metadata -> fastq_qc -> align -> markdup -> bqsr -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot
```

## 6. 只得到 germline gVCF / VCF，不跑 GWAS

如果只需要单样本 gVCF：

```bash
cd /400T/ckn/CFDNA_Breast_GWAS
bash 00.shell/run_pipeline.sh metadata
bash 00.shell/run_pipeline.sh bam_qc
bash 00.shell/run_pipeline.sh gvcf
```

单样本 gVCF 输出：

```text
06_gvcf/<sample_id>.g.vcf.gz
06_gvcf/<sample_id>.g.vcf.gz.tbi
```

如果需要多样本 joint germline VCF：

```bash
bash 00.shell/run_pipeline.sh joint
```

joint VCF 输出：

```text
07_joint_vcf/cohort_chr1.vcf.gz
07_joint_vcf/cohort_chr2.vcf.gz
...
07_joint_vcf/cohort_chr22.vcf.gz
```

如果只要 gVCF / joint VCF，可以不跑：

```text
plink
pca
gwas
plot
```

## 7. 后台运行

推荐使用 `nohup` 后台运行重任务。

### 跑 gVCF

```bash
cd /400T/ckn/CFDNA_Breast_GWAS
nohup bash 00.shell/step2_haplotypecaller_gvcf.sh \
  > loginfo/step2_haplotypecaller_gvcf.nohup.log \
  2>&1 &
```

查看主日志：

```bash
tail -f loginfo/step2_haplotypecaller_gvcf.nohup.log
```

查看某个样本日志：

```bash
tail -f loginfo/haplotypecaller_Lib-1005.log
```

查看是否还在跑：

```bash
ps -ef | grep HaplotypeCaller | grep -v grep
```

### 跑 joint VCF

```bash
cd /400T/ckn/CFDNA_Breast_GWAS
nohup bash 00.shell/step3_joint_genotype_by_chr.sh \
  > loginfo/step3_joint_genotype_by_chr.nohup.log \
  2>&1 &
```

查看日志：

```bash
tail -f loginfo/step3_joint_genotype_by_chr.nohup.log
```

## 8. 检查输出是否完整

检查 gVCF 数量：

```bash
find 06_gvcf -name *.g.vcf.gz | wc -l
find 06_gvcf -name *.g.vcf.gz.tbi | wc -l
```

两个数量应该一致。

检查某个 gVCF 样本名：

```bash
bcftools query -l 06_gvcf/Lib-1005.g.vcf.gz
```

检查 joint VCF：

```bash
ls -lh 07_joint_vcf/cohort_chr*.vcf.gz
```

检查 joint VCF 样本数：

```bash
bcftools query -l 07_joint_vcf/cohort_chr1.vcf.gz | wc -l
```

## 9. 半成品处理

如果 HaplotypeCaller 中途被 `Ctrl+C` 或系统中断，可能留下不完整文件。当前 Breast 脚本已经做了保护：只有同时存在 `.g.vcf.gz` 和 `.g.vcf.gz.tbi` 才会跳过。

手动清理没有索引的半成品：

```bash
cd /400T/ckn/CFDNA_Breast_GWAS/06_gvcf
for f in *.g.vcf.gz; do
  [ -f ${f}.tbi ] || rm -f $f
done
```

## 10. Breast 当前项目

Breast 项目目录：

```text
/400T/ckn/CFDNA_Breast_GWAS
```

Breast 样本数：

```text
349 个 BAM+BAI 完整样本
```

当前配置：

```text
/400T/ckn/CFDNA_Breast_GWAS/config.breast.bam.yaml
```

当前样本表：

```text
/400T/ckn/CFDNA_Breast_GWAS/samples.breast.bam.tsv
```

注意：当前 Breast 表型计数是：

```text
Breast Cancer: 348
Benign: 1
```

所以如果做 cancer vs benign GWAS，统计上不合适；但如果只是生成 germline gVCF / joint VCF，没有问题。

## 11. 常见问题

### cfdna-gwas-generate: command not found

先确认环境：

```bash
conda activate biosoft
which cfdna-gwas-generate
```

如果找不到，用完整路径：

```bash
/home/ckn/.local/bin/cfdna-gwas-generate -h
```

或者使用 cfDNA 环境：

```bash
conda activate cfDNA
cfdna-gwas-generate -h
```

### Ctrl+C 后任务还在跑

检查：

```bash
ps -ef | grep HaplotypeCaller | grep -v grep
```

停止指定 PID：

```bash
kill <PID>
```

不要随便 `killall java`，可能会影响别人的任务。

### 只想生成脚本，不运行分析

只运行：

```bash
cfdna-gwas-generate --config config.yaml
```

生成器只写项目目录，不会启动 HaplotypeCaller、BWA、GATK 等重任务。
