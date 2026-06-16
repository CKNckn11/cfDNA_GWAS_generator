#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import stat
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise SystemExit("PyYAML is required: python3 -m pip install pyyaml") from exc

VERSION = "0.1.0"


def q(value):
    return "'" + str(value).replace("'", "'\"'\"'") + "'"


def get(cfg, path, default=None):
    cur = cfg
    for key in path.split("."):
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def write_executable(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def normalize_chromosomes(chroms):
    if not chroms:
        return [f"chr{i}" for i in range(1, 23)]
    return [str(c) for c in chroms]


def required_columns(input_type):
    if input_type == "bam":
        return {"bam"}
    if input_type == "fastq":
        return {"fastq1", "fastq2"}
    raise SystemExit(f"Unsupported input.type: {input_type}. Use bam or fastq.")


def validate_samples(sample_table, input_type):
    with open(sample_table, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cols = set(reader.fieldnames or [])
        missing = required_columns(input_type) - cols
        if "sample_id" not in cols and "iid" not in cols:
            missing.add("sample_id_or_iid")
        if missing:
            raise SystemExit(f"Sample table missing columns: {', '.join(sorted(missing))}")
        rows = [r for r in reader if (r.get("sample_id") or r.get("iid") or "").strip()]
    if not rows:
        raise SystemExit("Sample table has no samples.")
    return rows


def config_sh(cfg, project_dir, input_type):
    tools = get(cfg, "tools", {})
    res = get(cfg, "resources", {})
    plink_qc = get(cfg, "plink_qc", {})
    gwas = get(cfg, "gwas", {})
    covars = get(gwas, "covars", ["age", "sex"] + [f"PC{i}" for i in range(1, 11)])
    if isinstance(covars, list):
        covars = ",".join(map(str, covars))
    chroms = " ".join(normalize_chromosomes(get(cfg, "reference.chromosomes")))
    known = get(cfg, "known_sites", {})

    return f"""#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR={q(project_dir)}
INPUT_TYPE={q(input_type)}

SHELL_DIR="${{PROJECT_DIR}}/00.shell"
META_DIR="${{PROJECT_DIR}}/00_metadata"
FASTQ_QC_DIR="${{PROJECT_DIR}}/01_fastq_qc"
ALIGN_DIR="${{PROJECT_DIR}}/02_alignment"
MARKDUP_DIR="${{PROJECT_DIR}}/03_markdup"
BQSR_DIR="${{PROJECT_DIR}}/04_bqsr"
BAM_QC_DIR="${{PROJECT_DIR}}/05_bam_qc"
GVCF_DIR="${{PROJECT_DIR}}/06_gvcf"
JOINT_DIR="${{PROJECT_DIR}}/07_joint_vcf"
PLINK_DIR="${{PROJECT_DIR}}/08_plink_qc"
PCA_DIR="${{PROJECT_DIR}}/09_pca"
GWAS_DIR="${{PROJECT_DIR}}/10_gwas"
PLOT_DIR="${{PROJECT_DIR}}/11_plots"
LOG_DIR="${{PROJECT_DIR}}/loginfo"

REFERENCE={q(get(cfg, "reference.fasta"))}
CHROMOSOMES={q(chroms)}
DBSNP={q(known.get("dbsnp", ""))}
KNOWN_INDELS={q(known.get("known_indels", ""))}
MILLS={q(known.get("mills", ""))}

FASTQC={q(tools.get("fastqc", "fastqc"))}
BWA={q(tools.get("bwa", "bwa"))}
SAMTOOLS={q(tools.get("samtools", "samtools"))}
GATK_JAR={q(tools.get("gatk_jar", ""))}
BCFTOOLS={q(tools.get("bcftools", "bcftools"))}
PLINK2={q(tools.get("plink2", "plink2"))}
RSCRIPT={q(tools.get("rscript", "Rscript"))}

BWA_THREADS="${{BWA_THREADS:-{res.get("bwa_threads", 8)}}}"
HC_THREADS="${{HC_THREADS:-{res.get("hc_threads", 2)}}}"
MAX_JOBS="${{MAX_JOBS:-{res.get("max_jobs", 4)}}}"
JAVA_MEM_MARKDUP="${{JAVA_MEM_MARKDUP:-{res.get("java_mem_markdup", "16G")}}}"
JAVA_MEM_BQSR="${{JAVA_MEM_BQSR:-{res.get("java_mem_bqsr", "16G")}}}"
JAVA_MEM_HC="${{JAVA_MEM_HC:-{res.get("java_mem_hc", "8G")}}}"
JAVA_MEM_JOINT="${{JAVA_MEM_JOINT:-{res.get("java_mem_joint", "16G")}}}"

PHENO_NAME="${{PHENO_NAME:-{get(cfg, "phenotype.name", "case_control")}}}"
CASE_LABEL={q(get(cfg, "phenotype.case_label", "case"))}
CONTROL_LABEL={q(get(cfg, "phenotype.control_label", "control"))}
COVAR_NAMES="${{COVAR_NAMES:-{covars}}}"

PLINK_GENO={q(plink_qc.get("geno", 0.05))}
PLINK_MIND={q(plink_qc.get("mind", 0.05))}
PLINK_MAF={q(plink_qc.get("maf", 0.01))}
PLINK_HWE={q(plink_qc.get("hwe", "1e-6"))}

mkdir -p "$SHELL_DIR" "$META_DIR" "$FASTQ_QC_DIR" "$ALIGN_DIR" "$MARKDUP_DIR" "$BQSR_DIR" \\
  "$BAM_QC_DIR" "$GVCF_DIR" "$JOINT_DIR" "$PLINK_DIR" "$PCA_DIR" "$GWAS_DIR" "$PLOT_DIR" "$LOG_DIR"
"""


def prepare_metadata_py(input_type):
    return f"""#!/usr/bin/env python3
import csv
import os

project = os.environ["PROJECT_DIR"]
input_type = os.environ["INPUT_TYPE"]
case_label = os.environ["CASE_LABEL"]
control_label = os.environ["CONTROL_LABEL"]
meta = os.environ["META_DIR"]
sample_table = os.path.join(meta, "input_samples.tsv")

def clean(x, default="NA"):
    x = (x or "").strip()
    return x if x else default

def sample_id(row):
    return clean(row.get("sample_id") or row.get("iid") or row.get("library_id"), "")

def phenotype(row):
    return clean(row.get("phenotype") or row.get("diagnosis_type"))

with open(sample_table, newline="", encoding="utf-8-sig") as handle:
    rows = [r for r in csv.DictReader(handle, delimiter="\\t") if clean(r.get("sample_id"), "")] 

manifest_cols = ["sample_id", "phenotype", "age", "sex"]
if input_type == "bam":
    manifest_cols += ["bam", "bai"]
else:
    manifest_cols += ["fastq1", "fastq2", "planned_bam"]

with open(os.path.join(meta, "sample_manifest.tsv"), "w", newline="", encoding="utf-8") as out:
    w = csv.DictWriter(out, fieldnames=manifest_cols, delimiter="\\t", extrasaction="ignore", lineterminator="\\n")
    w.writeheader()
    for r in rows:
        sid = sample_id(r)
        rec = {{
            "sample_id": sid,
            "phenotype": phenotype(r),
            "age": clean(r.get("age")),
            "sex": clean(r.get("sex")),
            "bam": clean(r.get("bam"), ""),
            "bai": clean(r.get("bai"), ""),
            "fastq1": clean(r.get("fastq1"), ""),
            "fastq2": clean(r.get("fastq2"), ""),
            "planned_bam": os.path.join(project, "04_bqsr", sid + ".recal.bam"),
        }}
        w.writerow(rec)

if input_type == "fastq":
    with open(os.path.join(meta, "fastq_list.tsv"), "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\\t", lineterminator="\\n")
        w.writerow(["iid", "fastq1", "fastq2"])
        for r in rows:
            w.writerow([sample_id(r), clean(r.get("fastq1"), ""), clean(r.get("fastq2"), "")])

with open(os.path.join(meta, "bam_list.tsv"), "w", newline="", encoding="utf-8") as out:
    w = csv.writer(out, delimiter="\\t", lineterminator="\\n")
    w.writerow(["iid", "bam"])
    for r in rows:
        sid = sample_id(r)
        bam = clean(r.get("bam"), "") if input_type == "bam" else os.path.join(project, "04_bqsr", sid + ".recal.bam")
        w.writerow([sid, bam])

with open(os.path.join(meta, "phenotype.txt"), "w", newline="", encoding="utf-8") as out:
    w = csv.writer(out, delimiter="\\t", lineterminator="\\n")
    w.writerow(["FID", "IID", os.environ["PHENO_NAME"]])
    for r in rows:
        sid = sample_id(r)
        label = phenotype(r)
        pheno = "2" if label == case_label else "1" if label == control_label else "-9"
        w.writerow([sid, sid, pheno])

with open(os.path.join(meta, "covariates_base.txt"), "w", newline="", encoding="utf-8") as out:
    w = csv.writer(out, delimiter="\\t", lineterminator="\\n")
    w.writerow(["FID", "IID", "age", "sex"])
    for r in rows:
        sid = sample_id(r)
        w.writerow([sid, sid, clean(r.get("age")), clean(r.get("sex"))])

print("samples:", len(rows))
print("metadata:", meta)
"""


def common_merge_covariates_py():
    return """#!/usr/bin/env python3
import csv
import os

project = os.environ["PROJECT_DIR"]
base = os.path.join(project, "00_metadata", "covariates_base.txt")
pcs = os.path.join(project, "09_pca", "pca.eigenvec")
out = os.path.join(project, "09_pca", "covariates_with_pcs.txt")

cov = {}
with open(base, newline="") as handle:
    for r in csv.DictReader(handle, delimiter="\t"):
        cov[(r["FID"], r["IID"])] = r

with open(pcs, newline="") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))

with open(out, "w", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["FID", "IID", "age", "sex"] + [f"PC{i}" for i in range(1, 21)])
    for row in rows:
        if len(row) < 22:
            row = row[0].split()
        fid, iid = row[0], row[1]
        c = cov.get((fid, iid), {"age": "NA", "sex": "NA"})
        w.writerow([fid, iid, c["age"], c["sex"]] + row[2:22])
print(out)
"""


def step0():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
PROJECT_DIR="$PROJECT_DIR" INPUT_TYPE="$INPUT_TYPE" META_DIR="$META_DIR" CASE_LABEL="$CASE_LABEL" CONTROL_LABEL="$CONTROL_LABEL" PHENO_NAME="$PHENO_NAME" \
  python3 "${SHELL_DIR}/prepare_metadata.py"
"""


def fastq_steps():
    return {
        "step1_fastq_qc.sh": """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
tail -n +2 "${META_DIR}/fastq_list.tsv" | while IFS=$'\t' read -r iid fq1 fq2; do
  echo "[FastQC] $iid"
  mkdir -p "${FASTQ_QC_DIR}/${iid}"
  "$FASTQC" -o "${FASTQ_QC_DIR}/${iid}" "$fq1" "$fq2" > "${LOG_DIR}/fastqc_${iid}.log" 2>&1
done
""",
        "step2_bwa_align.sh": """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
tail -n +2 "${META_DIR}/fastq_list.tsv" | while IFS=$'\t' read -r iid fq1 fq2; do
  out="${ALIGN_DIR}/${iid}.sorted.bam"
  [ -s "$out" ] && { echo "[skip] $iid sorted BAM exists"; continue; }
  echo "[BWA] $iid"
  "$BWA" mem -t "$BWA_THREADS" -R "@RG\tID:${iid}\tSM:${iid}\tPL:ILLUMINA" "$REFERENCE" "$fq1" "$fq2" 2> "${LOG_DIR}/bwa_${iid}.log" | \
    "$SAMTOOLS" sort -@ "$BWA_THREADS" -o "$out" -
  "$SAMTOOLS" index "$out"
done
""",
        "step3_markdup.sh": """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
tail -n +2 "${META_DIR}/fastq_list.tsv" | while IFS=$'\t' read -r iid fq1 fq2; do
  in="${ALIGN_DIR}/${iid}.sorted.bam"
  out="${MARKDUP_DIR}/${iid}.markdup.bam"
  metrics="${MARKDUP_DIR}/${iid}.markdup.metrics.txt"
  [ -s "$out" ] && { echo "[skip] $iid markdup BAM exists"; continue; }
  echo "[MarkDuplicates] $iid"
  java -Xmx"$JAVA_MEM_MARKDUP" -jar "$GATK_JAR" MarkDuplicates -I "$in" -O "$out" -M "$metrics" > "${LOG_DIR}/markdup_${iid}.log" 2>&1
  "$SAMTOOLS" index "$out"
done
""",
        "step4_bqsr.sh": """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
known_args=()
[ -n "$DBSNP" ] && known_args+=(--known-sites "$DBSNP")
[ -n "$KNOWN_INDELS" ] && known_args+=(--known-sites "$KNOWN_INDELS")
[ -n "$MILLS" ] && known_args+=(--known-sites "$MILLS")
tail -n +2 "${META_DIR}/fastq_list.tsv" | while IFS=$'\t' read -r iid fq1 fq2; do
  in="${MARKDUP_DIR}/${iid}.markdup.bam"
  table="${BQSR_DIR}/${iid}.recal.table"
  out="${BQSR_DIR}/${iid}.recal.bam"
  [ -s "$out" ] && { echo "[skip] $iid BQSR BAM exists"; continue; }
  echo "[BQSR] $iid"
  java -Xmx"$JAVA_MEM_BQSR" -jar "$GATK_JAR" BaseRecalibrator -R "$REFERENCE" -I "$in" "${known_args[@]}" -O "$table" > "${LOG_DIR}/bqsr_${iid}.log" 2>&1
  java -Xmx"$JAVA_MEM_BQSR" -jar "$GATK_JAR" ApplyBQSR -R "$REFERENCE" -I "$in" --bqsr-recal-file "$table" -O "$out" >> "${LOG_DIR}/bqsr_${iid}.log" 2>&1
  "$SAMTOOLS" index "$out"
done
""",
    }


def bam_qc_script():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
tail -n +2 "${META_DIR}/bam_list.tsv" | while IFS=$'\t' read -r iid bam; do
  out_dir="${BAM_QC_DIR}/${iid}"
  mkdir -p "$out_dir"
  echo "[BAM QC] $iid"
  "$SAMTOOLS" quickcheck -v "$bam" > "${out_dir}/${iid}.quickcheck.txt" 2>&1 || true
  "$SAMTOOLS" flagstat "$bam" > "${out_dir}/${iid}.flagstat.txt"
  "$SAMTOOLS" stats "$bam" > "${out_dir}/${iid}.stats.txt"
done
"""


def hc_script():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_hc() {
  local iid="$1"
  local bam="$2"
  local out="${GVCF_DIR}/${iid}.g.vcf.gz"
  local log="${LOG_DIR}/haplotypecaller_${iid}.log"
  [ -s "$out" ] && { echo "[skip] $iid gVCF exists"; return 0; }
  echo "[HaplotypeCaller] $iid"
  java -Xmx"$JAVA_MEM_HC" -jar "$GATK_JAR" HaplotypeCaller \
    -R "$REFERENCE" -I "$bam" -O "$out" -ERC GVCF \
    --native-pair-hmm-threads "$HC_THREADS" > "$log" 2>&1
}
export -f run_hc
export GVCF_DIR LOG_DIR GATK_JAR REFERENCE JAVA_MEM_HC HC_THREADS

jobs=0
tail -n +2 "${META_DIR}/bam_list.tsv" | while IFS=$'\t' read -r iid bam; do
  run_hc "$iid" "$bam" &
  jobs=$((jobs + 1))
  if [ "$jobs" -ge "$MAX_JOBS" ]; then
    wait -n
    jobs=$((jobs - 1))
  fi
done
wait
"""


def joint_script():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
map_file="${JOINT_DIR}/sample_gvcf_map.tsv"
find "$GVCF_DIR" -name "*.g.vcf.gz" | sort | awk -F/ '{f=$NF; sub(/\\.g\\.vcf\\.gz$/,"",f); print f "\\t" $0}' > "$map_file"
[ -s "$map_file" ] || { echo "No gVCFs found in $GVCF_DIR"; exit 1; }

for interval in $CHROMOSOMES; do
  db="${JOINT_DIR}/genomicsdb_${interval}"
  out="${JOINT_DIR}/cohort_${interval}.vcf.gz"
  log="${LOG_DIR}/joint_${interval}.log"
  if [ ! -d "$db" ]; then
    echo "[GenomicsDBImport] $interval"
    java -Xmx"$JAVA_MEM_JOINT" -jar "$GATK_JAR" GenomicsDBImport --sample-name-map "$map_file" --genomicsdb-workspace-path "$db" -L "$interval" > "$log" 2>&1
  fi
  if [ ! -s "$out" ]; then
    echo "[GenotypeGVCFs] $interval"
    java -Xmx"$JAVA_MEM_JOINT" -jar "$GATK_JAR" GenotypeGVCFs -R "$REFERENCE" -V "gendb://${db}" -O "$out" >> "$log" 2>&1
  fi
done
"""


def plink_script():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
vcf_list="${JOINT_DIR}/autosome_vcfs.txt"
ls "${JOINT_DIR}"/cohort_chr*.vcf.gz | sort -V > "$vcf_list"
[ -s "$vcf_list" ] || { echo "No joint VCF files found."; exit 1; }
merged="${JOINT_DIR}/cohort_autosomes.vcf.gz"
if [ ! -s "$merged" ]; then
  "$BCFTOOLS" concat -a -Oz -o "$merged" -f "$vcf_list"
  "$BCFTOOLS" index -t "$merged"
fi
"$PLINK2" --vcf "$merged" --max-alleles 2 --snps-only just-acgt --make-pgen --out "${PLINK_DIR}/cohort_raw"
"$PLINK2" --pfile "${PLINK_DIR}/cohort_raw" --geno "$PLINK_GENO" --mind "$PLINK_MIND" --maf "$PLINK_MAF" --hwe "$PLINK_HWE" midp --make-bed --out "${PLINK_DIR}/cohort_qc"
"""


def pca_script():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
base="${PLINK_DIR}/cohort_qc"
[ -f "${base}.bed" ] || { echo "Missing ${base}.bed"; exit 1; }
"$PLINK2" --bfile "$base" --indep-pairwise 50 5 0.2 --out "${PCA_DIR}/prune"
"$PLINK2" --bfile "$base" --extract "${PCA_DIR}/prune.prune.in" --pca 20 --out "${PCA_DIR}/pca"
PROJECT_DIR="$PROJECT_DIR" python3 "${SHELL_DIR}/merge_covariates_with_pcs.py"
"""


def gwas_script():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
base="${PLINK_DIR}/cohort_qc"
pheno="${META_DIR}/phenotype.txt"
covar="${PCA_DIR}/covariates_with_pcs.txt"
[ -f "${base}.bed" ] || { echo "Missing ${base}.bed"; exit 1; }
[ -f "$pheno" ] || { echo "Missing $pheno"; exit 1; }
[ -f "$covar" ] || { echo "Missing $covar"; exit 1; }
"$PLINK2" --bfile "$base" --pheno "$pheno" --pheno-name "$PHENO_NAME" --covar "$covar" --covar-name "$COVAR_NAMES" --glm hide-covar firth-fallback --out "${GWAS_DIR}/${PHENO_NAME}"
"""


def plot_r():
    return """#!/usr/bin/env Rscript
project <- Sys.getenv("PROJECT_DIR", unset = "/tmp")
gwas_dir <- file.path(project, "10_gwas")
plot_dir <- file.path(project, "11_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
files <- list.files(gwas_dir, pattern = "\\\\.glm\\\\.", full.names = TRUE)
if (length(files) == 0) stop("No PLINK2 .glm result files found")
d <- read.table(files[1], header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
if (!"P" %in% names(d)) stop("No P column in GWAS result")
d <- d[!is.na(d$P) & d$P > 0, ]
chr_col <- if ("X.CHROM" %in% names(d)) "X.CHROM" else "#CHROM"
d$CHR_NUM <- suppressWarnings(as.integer(gsub("^chr", "", as.character(d[[chr_col]]))))
d <- d[!is.na(d$CHR_NUM) & d$CHR_NUM >= 1 & d$CHR_NUM <= 22, ]
d <- d[order(d$CHR_NUM, d$POS), ]
chr_len <- tapply(d$POS, d$CHR_NUM, max)
offset <- cumsum(c(0, head(chr_len, -1)))
names(offset) <- names(chr_len)
d$BPcum <- d$POS + offset[as.character(d$CHR_NUM)]
png(file.path(plot_dir, "manhattan.png"), width = 1800, height = 900, res = 150)
cols <- c("#4C78A8", "#F58518")
plot(d$BPcum, -log10(d$P), pch = 20, cex = 0.45, col = cols[(d$CHR_NUM %% 2) + 1], xaxt = "n", xlab = "Chromosome", ylab = "-log10(P)", main = "cfDNA germline association")
axis(1, at = tapply(d$BPcum, d$CHR_NUM, median), labels = names(chr_len), cex.axis = 0.75)
abline(h = -log10(5e-8), col = "red", lty = 2)
abline(h = -log10(1e-5), col = "darkgray", lty = 3)
dev.off()
png(file.path(plot_dir, "qqplot.png"), width = 900, height = 900, res = 150)
obs <- sort(d$P)
exp <- ppoints(length(obs))
plot(-log10(exp), -log10(obs), pch = 20, cex = 0.6, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = "QQ plot")
abline(0, 1, col = "red", lty = 2)
dev.off()
message("Plots saved to: ", plot_dir)
"""


def run_pipeline(input_type, steps):
    cases = ["  help|*) usage ;;"]
    usage_lines = ["Usage: bash 00.shell/run_pipeline.sh <mode>", "", "Modes:"]
    all_names = []
    for mode, script in steps:
        all_names.append(mode)
        usage_lines.append(f"  {mode}")
        cases.insert(-1, f"""  {mode})
    run_step {q(mode)} "bash 00.shell/{script}"
    ;;""")
    usage_lines.append("  all")
    all_cmds = "\n".join([f"    run_step {q(m)} \"bash 00.shell/{s}\"" for m, s in steps])
    cases.insert(-1, f"""  all)
{all_cmds}
    ;;""")
    usage = "\\n".join(usage_lines)
    return f"""#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${{SCRIPT_DIR}}/config.sh"
MODE="${{1:-help}}"

run_step() {{
  local name="$1"
  local cmd="$2"
  local log="${{LOG_DIR}}/${{name}}.$(date +%Y%m%d_%H%M%S).log"
  echo "=== [$name] start: $(date) ==="
  echo "command: $cmd"
  bash -lc "cd '$PROJECT_DIR' && $cmd" > "$log" 2>&1
  echo "=== [$name] done: $(date) ==="
  echo "log: $log"
}}

usage() {{
  printf '%b\\n' {q(usage)}
}}

case "$MODE" in
{os.linesep.join(cases)}
esac
"""


def readme(input_type):
    if input_type == "bam":
        order = "metadata -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot"
    else:
        order = "metadata -> fastq_qc -> align -> markdup -> bqsr -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot"
    return f"""# Generated cfDNA GWAS Pipeline

Input type: `{input_type}`

Run order:

```text
{order}
```

Run all steps:

```bash
bash 00.shell/run_pipeline.sh all
```

Run one step:

```bash
bash 00.shell/run_pipeline.sh gvcf
```

Logs are written to `loginfo/`. Edit `config.sh` only if paths or resources need to be changed after generation.
"""


def generate(cfg, config_path, outdir):
    input_type = get(cfg, "input.type")
    sample_table = get(cfg, "input.sample_table")
    if not sample_table:
        raise SystemExit("Missing input.sample_table in config.")
    rows = validate_samples(sample_table, input_type)
    project_dir = str(Path(outdir or get(cfg, "project.outdir")).resolve())
    project = Path(project_dir)
    shell_dir = project / "00.shell"

    project.mkdir(parents=True, exist_ok=True)
    (project / "00_metadata").mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, project / "config.yaml")
    shutil.copy2(sample_table, project / "00_metadata" / "input_samples.tsv")
    write_executable(project / "config.sh", config_sh(cfg, project_dir, input_type))
    write_executable(shell_dir / "prepare_metadata.py", prepare_metadata_py(input_type))
    write_executable(shell_dir / "merge_covariates_with_pcs.py", common_merge_covariates_py())
    write_executable(shell_dir / "step0_prepare_metadata.sh", step0())

    steps = [("metadata", "step0_prepare_metadata.sh")]
    if input_type == "fastq":
        for name, text in fastq_steps().items():
            write_executable(shell_dir / name, text)
        steps += [
            ("fastq_qc", "step1_fastq_qc.sh"),
            ("align", "step2_bwa_align.sh"),
            ("markdup", "step3_markdup.sh"),
            ("bqsr", "step4_bqsr.sh"),
        ]
        bam_qc_name = "step5_bam_qc.sh"
        hc_name = "step6_haplotypecaller_gvcf.sh"
        joint_name = "step7_joint_genotype_by_chr.sh"
        plink_name = "step8_make_plink_and_qc.sh"
        pca_name = "step9_pca.sh"
        gwas_name = "step10_run_gwas.sh"
        plot_name = "step11_plot_gwas.R"
    else:
        bam_qc_name = "step1_bam_qc.sh"
        hc_name = "step2_haplotypecaller_gvcf.sh"
        joint_name = "step3_joint_genotype_by_chr.sh"
        plink_name = "step4_make_plink_and_qc.sh"
        pca_name = "step5_pca.sh"
        gwas_name = "step6_run_gwas.sh"
        plot_name = "step7_plot_gwas.R"

    write_executable(shell_dir / bam_qc_name, bam_qc_script())
    write_executable(shell_dir / hc_name, hc_script())
    write_executable(shell_dir / joint_name, joint_script())
    write_executable(shell_dir / plink_name, plink_script())
    write_executable(shell_dir / pca_name, pca_script())
    write_executable(shell_dir / gwas_name, gwas_script())
    write_executable(shell_dir / plot_name, plot_r())

    steps += [
        ("bam_qc", bam_qc_name),
        ("gvcf", hc_name),
        ("joint", joint_name),
        ("plink", plink_name),
        ("pca", pca_name),
        ("gwas", gwas_name),
        ("plot", plot_name),
    ]
    write_executable(shell_dir / "run_pipeline.sh", run_pipeline(input_type, steps))
    write_text(project / "README.md", readme(input_type))
    print(f"Generated pipeline: {project}")
    print(f"Samples: {len(rows)}")
    print(f"Run: cd {project} && bash 00.shell/run_pipeline.sh all")


def main():
    parser = argparse.ArgumentParser(
        prog="cfdna-gwas-generate",
        description=(
            f"cfDNA GWAS Pipeline Generator (Version = {VERSION}): "
            "Generate cfDNA germline GWAS shell scripts from a YAML config."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""example:
  cfdna-gwas-generate --config examples/config.bam.example.yaml
  python generate_pipeline.py --config examples/config.bam.example.yaml
""",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"cfdna-gwas-generate {VERSION}",
        help="show the version of cfdna-gwas-generate and exit.",
    )
    parser.add_argument("--config", required=True, metavar="CONFIG", help="YAML config file.")
    parser.add_argument("--outdir", metavar="OUTDIR", help="Override project.outdir from the YAML config.")
    args = parser.parse_args()

    with open(args.config, encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    generate(cfg, args.config, args.outdir)


if __name__ == "__main__":
    main()
