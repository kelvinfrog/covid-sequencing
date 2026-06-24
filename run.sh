#!/usr/bin/env bash
# Master entry point for the COVID sequencing Freyja pipeline.
#
# Just drop files and run:
#   data/accessions/   ← drop a new .txt file with SRA accession numbers
#   data/metadata/     ← drop a new .csv file with sample metadata
#
# Then run: bash run.sh
#
# The script will:
#   1. Read only the NEWEST .txt in data/accessions/ for this batch
#   2. Download any accessions not yet in data/raw/
#   3. Run the Freyja pipeline on new samples (already-done ones are skipped)
#   4. Use the NEWEST .csv in data/metadata/ for time-series plotting
#   5. Aggregate ALL results ever processed (cumulative) and produce plots

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="data/raw"
RESULTS_DIR="data/results"
ACCESSIONS_DIR="data/accessions"
METADATA_DIR="data/metadata"
THREADS=4

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Collect all unique accessions across all .txt files
# ─────────────────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════"
echo " Step 1 — Picking newest accession file"
echo "════════════════════════════════════════════════════════════════════"

# Pick only the most recently modified .txt file
ACC_FILE=$(
    find "${ACCESSIONS_DIR}" -maxdepth 1 -name "*.txt" -print0 \
    | xargs -0 stat -f "%m %N" 2>/dev/null \
    | sort -rn \
    | head -1 \
    | cut -d' ' -f2-
)

if [[ -z "${ACC_FILE}" ]]; then
    echo "ERROR: No .txt files found in ${ACCESSIONS_DIR}/"
    echo "       Drop a file with one SRA accession per line and re-run."
    exit 1
fi

echo "Using: $(basename "${ACC_FILE}")"

ACCESSIONS=()
while IFS= read -r acc; do
    [[ -z "$acc" ]] && continue
    ACCESSIONS+=("$acc")
done < <({ cat "${ACC_FILE}"; echo; } | sed 's/\r//' | grep -vE '^\s*(#|$)')

echo "Accessions in this batch: ${#ACCESSIONS[@]}"

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Download any accessions not yet in data/raw/
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════════"
echo " Step 2 — Downloading missing samples"
echo "════════════════════════════════════════════════════════════════════"

mkdir -p "${RAW_DIR}"
FAILED_DOWNLOADS=()
DOWNLOADED=0
ALREADY_HAVE=0

for ACC in "${ACCESSIONS[@]}"; do
    R1="${RAW_DIR}/${ACC}_R1.fastq.gz"
    R2="${RAW_DIR}/${ACC}_R2.fastq.gz"

    if [[ -f "${R1}" && -f "${R2}" ]]; then
        echo "[EXISTS]   ${ACC}"
        ALREADY_HAVE=$((ALREADY_HAVE + 1))
        continue
    fi

    echo "[DOWNLOAD] ${ACC}..."
    if prefetch "${ACC}" --output-directory "${RAW_DIR}" --progress \
    && fasterq-dump "${RAW_DIR}/${ACC}" \
           --outdir "${RAW_DIR}" \
           --split-files \
           --threads "${THREADS}" \
           --progress; then
        gzip -f "${RAW_DIR}/${ACC}_1.fastq" "${RAW_DIR}/${ACC}_2.fastq"
        mv "${RAW_DIR}/${ACC}_1.fastq.gz" "${R1}"
        mv "${RAW_DIR}/${ACC}_2.fastq.gz" "${R2}"
        rm -rf "${RAW_DIR}/${ACC}"
        echo "  Done."
        DOWNLOADED=$((DOWNLOADED + 1))
    else
        echo "  ERROR: Download failed for ${ACC}"
        FAILED_DOWNLOADS+=("${ACC}")
        rm -rf "${RAW_DIR}/${ACC}"
    fi
done

echo ""
echo "Download summary: ${DOWNLOADED} new, ${ALREADY_HAVE} already present, ${#FAILED_DOWNLOADS[@]} failed"
if [[ ${#FAILED_DOWNLOADS[@]} -gt 0 ]]; then
    echo "Failed: ${FAILED_DOWNLOADS[*]}"
    echo "Re-run to retry failed downloads."
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Auto-detect most recently modified metadata CSV
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════════"
echo " Step 3 — Detecting metadata"
echo "════════════════════════════════════════════════════════════════════"

METADATA_FILE=""
if [[ -d "${METADATA_DIR}" ]]; then
    # Sort by modification time, pick newest
    METADATA_FILE=$(
        find "${METADATA_DIR}" -maxdepth 1 -name "*.csv" -print0 \
        | xargs -0 stat -f "%m %N" 2>/dev/null \
        | sort -rn \
        | head -1 \
        | cut -d' ' -f2-
    )
fi

if [[ -n "${METADATA_FILE}" ]]; then
    echo "Using: ${METADATA_FILE}"
else
    echo "No .csv found in ${METADATA_DIR}/ — time-series plot will be skipped."
    echo "Drop a metadata CSV there to enable it."
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Run pipeline + aggregate + plot
# ─────────────────────────────────────────────────────────────────────────────
echo ""
bash "${SCRIPT_DIR}/run_all.sh" "${RAW_DIR}" "${RESULTS_DIR}" "${METADATA_FILE}"
