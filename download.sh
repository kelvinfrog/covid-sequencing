#!/usr/bin/env bash
# Download SRA accessions and prepare paired FASTQs for the Freyja pipeline.
#
# Usage: bash download.sh <accessions.txt> [RAW_DIR]
#
# accessions.txt: plain text file, one SRA accession per line (e.g. SRR12345678)
#
# Example:
#   bash download.sh data/accessions.txt data/raw

set -euo pipefail

ACCESSION_FILE="${1:?Usage: $0 <accessions.txt> [RAW_DIR]}"
RAW_DIR="${2:-data/raw}"
THREADS=4

if [[ ! -f "${ACCESSION_FILE}" ]]; then
    echo "Accession file not found: ${ACCESSION_FILE}"
    exit 1
fi

mkdir -p "${RAW_DIR}"

# Strip blank lines and comments
ACCESSIONS=()
while IFS= read -r line; do
    line="${line//[$'\r\n']/}"   # strip carriage returns
    [[ -z "${line}" || "${line}" == \#* ]] && continue
    ACCESSIONS+=("${line}")
done < "${ACCESSION_FILE}"

TOTAL=${#ACCESSIONS[@]}
echo "Found ${TOTAL} accession(s) to download"
echo ""

FAILED=()

for i in "${!ACCESSIONS[@]}"; do
    ACC="${ACCESSIONS[$i]}"
    NUM=$((i + 1))
    echo "────────────────────────────────────────────────────────────────────"
    echo "[${NUM}/${TOTAL}] ${ACC}"
    echo "────────────────────────────────────────────────────────────────────"

    R1="${RAW_DIR}/${ACC}_R1.fastq.gz"
    R2="${RAW_DIR}/${ACC}_R2.fastq.gz"

    # Skip if already downloaded
    if [[ -f "${R1}" && -f "${R2}" ]]; then
        echo "  Already exists — skipping"
        continue
    fi

    # Step 1: prefetch (downloads .sra file)
    echo "  Downloading .sra..."
    if ! prefetch "${ACC}" --output-directory "${RAW_DIR}" --progress; then
        echo "  ERROR: prefetch failed for ${ACC}"
        FAILED+=("${ACC}")
        continue
    fi

    # Step 2: convert to FASTQ
    echo "  Converting to FASTQ..."
    if ! fasterq-dump "${RAW_DIR}/${ACC}" \
            --outdir "${RAW_DIR}" \
            --split-files \
            --threads "${THREADS}" \
            --progress; then
        echo "  ERROR: fasterq-dump failed for ${ACC}"
        FAILED+=("${ACC}")
        rm -rf "${RAW_DIR}/${ACC}"
        continue
    fi

    # Step 3: gzip
    echo "  Compressing..."
    gzip -f "${RAW_DIR}/${ACC}_1.fastq" "${RAW_DIR}/${ACC}_2.fastq"

    # Step 4: rename to _R1/_R2 convention expected by run_all.sh
    mv "${RAW_DIR}/${ACC}_1.fastq.gz" "${R1}"
    mv "${RAW_DIR}/${ACC}_2.fastq.gz" "${R2}"

    # Step 5: clean up intermediate .sra file
    rm -rf "${RAW_DIR}/${ACC}"

    echo "  Done → $(basename ${R1}), $(basename ${R2})"
done

echo ""
echo "════════════════════════════════════════════════════════════════════"
echo "Download complete: $((TOTAL - ${#FAILED[@]}))/${TOTAL} succeeded"

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo "Failed accessions:"
    for ACC in "${FAILED[@]}"; do
        echo "  - ${ACC}"
    done
    echo ""
    echo "Re-run this script to retry failed downloads (successful ones are skipped)."
    exit 1
fi

echo ""
echo "All FASTQs are in: ${RAW_DIR}"
echo "Run the pipeline with: bash run_all.sh"
