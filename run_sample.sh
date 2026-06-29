#!/usr/bin/env bash
# Freyja SARS-CoV-2 lineage abundance pipeline
# Platform : Illumina paired-end
# Primer scheme: ARTIC v5.3.2 (default)
#
# Usage: bash run_sample.sh <SAMPLE_ID> <R1.fastq.gz> <R2.fastq.gz> [OUTPUT_DIR]
#
# Example:
#   bash run_sample.sh sample01 \
#       data/raw/sample01_R1.fastq.gz \
#       data/raw/sample01_R2.fastq.gz \
#       data/results

set -euo pipefail

# ── Arguments ────────────────────────────────────────────────────────────────
SAMPLE="${1:?Usage: $0 <SAMPLE_ID> <R1.fastq.gz> <R2.fastq.gz> [OUTPUT_DIR]}"
R1="${2:?Missing R1 FASTQ}"
R2="${3:?Missing R2 FASTQ}"
OUTDIR="${4:-data/results}"

# ── Bundled references (no download needed) ───────────────────────────────────
FREYJA_DATA="$(python -c "import freyja; import os; print(os.path.join(os.path.dirname(freyja.__file__), 'data'))")"
REF="${FREYJA_DATA}/NC_045512_Hu-1.fasta"
GFF="${FREYJA_DATA}/NC_045512_Hu-1.gff"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRIMERS="${SCRIPT_DIR}/data/bed/ARTIC_V5.3.2.bed"

# ── Settings ──────────────────────────────────────────────────────────────────
THREADS=4
MIN_QUALITY=20
MIN_DEPTH=10
MIN_FREQ=0.03        # freyja variant frequency threshold
MAX_DEPTH=600000     # mpileup cap (keeps memory sane on deep wastewater)
# ── Skip if already processed ─────────────────────────────────────────────────
FREYJA_OUT="${OUTDIR}/${SAMPLE}/${SAMPLE}.freyja.tsv"
if [[ -f "${FREYJA_OUT}" ]]; then
    echo "[SKIP] ${SAMPLE} — already processed (freyja output exists)"
    exit 0
fi

# ── Output layout ─────────────────────────────────────────────────────────────
SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
mkdir -p "${SAMPLE_DIR}"/{trimmed,aligned,variants}

LOG="${SAMPLE_DIR}/${SAMPLE}.log"
exec > >(tee -a "${LOG}") 2>&1
echo "[$(date)] Starting pipeline for sample: ${SAMPLE}"

# ── Step 1: Quality trim with fastp ───────────────────────────────────────────
echo "[$(date)] Step 1/6 — Adapter trimming (fastp)"
fastp \
    --in1 "${R1}" \
    --in2 "${R2}" \
    --out1 "${SAMPLE_DIR}/trimmed/${SAMPLE}_R1.fastq.gz" \
    --out2 "${SAMPLE_DIR}/trimmed/${SAMPLE}_R2.fastq.gz" \
    --json "${SAMPLE_DIR}/trimmed/${SAMPLE}_fastp.json" \
    --html "${SAMPLE_DIR}/trimmed/${SAMPLE}_fastp.html" \
    --thread "${THREADS}" \
    --qualified_quality_phred "${MIN_QUALITY}" \
    --length_required 50 \
    --detect_adapter_for_pe

# ── Step 2: Align to SARS-CoV-2 reference ────────────────────────────────────
echo "[$(date)] Step 2/6 — Alignment (minimap2)"
minimap2 \
    -ax sr \
    -t "${THREADS}" \
    "${REF}" \
    "${SAMPLE_DIR}/trimmed/${SAMPLE}_R1.fastq.gz" \
    "${SAMPLE_DIR}/trimmed/${SAMPLE}_R2.fastq.gz" \
| samtools sort -@ "${THREADS}" -o "${SAMPLE_DIR}/aligned/${SAMPLE}.sorted.bam"

samtools index "${SAMPLE_DIR}/aligned/${SAMPLE}.sorted.bam"

# ── Step 3: Primer trimming with ivar ────────────────────────────────────────
echo "[$(date)] Step 3/6 — Primer trimming (ivar trim)"
ivar trim \
    -i "${SAMPLE_DIR}/aligned/${SAMPLE}.sorted.bam" \
    -b "${PRIMERS}" \
    -p "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed" \
    -q "${MIN_QUALITY}" \
    -m 50 \
    -s 4

samtools sort -@ "${THREADS}" \
    -o "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.sorted.bam" \
    "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.bam"
samtools index "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.sorted.bam"

# Clean up unsorted trimmed BAM
rm -f "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.bam"

# ── Step 4: Alignment stats ───────────────────────────────────────────────────
echo "[$(date)] Step 4/6 — Alignment statistics"
samtools flagstat "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.sorted.bam" \
    > "${SAMPLE_DIR}/aligned/${SAMPLE}.flagstat.txt"
samtools coverage "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.sorted.bam" \
    > "${SAMPLE_DIR}/aligned/${SAMPLE}.coverage.txt"

# ── Step 5: Variant calling + depth (freyja variants) ────────────────────────
echo "[$(date)] Step 5/6 — Variant calling + depth (freyja variants)"
freyja variants \
    "${SAMPLE_DIR}/aligned/${SAMPLE}.trimmed.sorted.bam" \
    --variants "${SAMPLE_DIR}/variants/${SAMPLE}" \
    --depths "${SAMPLE_DIR}/variants/${SAMPLE}.depths.tsv" \
    --ref "${REF}" \
    --minq "${MIN_QUALITY}"

# ── Step 6: Freyja demix ──────────────────────────────────────────────────────
echo "[$(date)] Step 6/6 — Lineage demixing (freyja demix)"
freyja demix \
    "${SAMPLE_DIR}/variants/${SAMPLE}.tsv" \
    "${SAMPLE_DIR}/variants/${SAMPLE}.depths.tsv" \
    --output "${SAMPLE_DIR}/${SAMPLE}.freyja.tsv" \
    --confirmedonly \
    --eps 0.01 \
    --covcut 10

echo "[$(date)] Done. Results → ${SAMPLE_DIR}/${SAMPLE}.freyja.tsv"
echo ""
echo "── Lineage summary ──────────────────────────────────────────────────────"
cat "${SAMPLE_DIR}/${SAMPLE}.freyja.tsv"
