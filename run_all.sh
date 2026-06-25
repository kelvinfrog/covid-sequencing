#!/usr/bin/env bash
# Batch runner — processes samples and produces per-batch + all-batches aggregates.
#
# Usage: bash run_all.sh <RAW_DIR> <RESULTS_DIR> <METADATA_FILE> <BATCH_NAME> <ACC_FILE>

set -euo pipefail

RAW_DIR="${1:-data/raw}"
RESULTS_DIR="${2:-data/results}"
METADATA_FILE="${3:-}"
BATCH_NAME="${4:-unknown_batch}"
ACC_FILE="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AGGREGATE_ROOT="${RESULTS_DIR}/_aggregate"

# ── Helper: aggregate + plot whatever .freyja.tsv files are in STAGING_DIR ───
aggregate_and_plot() {
    local STAGING_DIR="$1"
    local OUT_DIR="$2"
    local LABEL="$3"

    local COUNT
    COUNT=$(find "${STAGING_DIR}" -maxdepth 1 -name "*.freyja.tsv" | wc -l | tr -d ' ')
    if [[ "${COUNT}" -eq 0 ]]; then
        echo "  [${LABEL}] No processed samples found — skipping"
        return
    fi

    echo "  [${LABEL}] Aggregating ${COUNT} sample(s)..."
    mkdir -p "${OUT_DIR}"

    freyja aggregate \
        "${STAGING_DIR}/" \
        --output "${OUT_DIR}/aggregated.tsv" \
        --ext freyja.tsv

    freyja plot \
        "${OUT_DIR}/aggregated.tsv" \
        --output "${OUT_DIR}/lineage_plot.pdf" \
        --mincov 0

    echo "  [${LABEL}] lineage_plot.pdf ✓"

    # Time-series plot — only if metadata with valid dates is available
    if [[ -n "${METADATA_FILE}" && -f "${METADATA_FILE}" ]]; then
        local TIMES_CSV="${OUT_DIR}/_times_metadata.csv"
        local FILTERED_AGG="${OUT_DIR}/_aggregated_timed.tsv"
        rm -f "${TIMES_CSV}" "${FILTERED_AGG}"

        python3 - <<PYEOF
import csv, sys
from datetime import datetime

metadata_file = "${METADATA_FILE}"
times_csv     = "${TIMES_CSV}"
agg_file      = "${OUT_DIR}/aggregated.tsv"
filtered_agg  = "${FILTERED_AGG}"

with open(metadata_file) as f:
    reader = csv.DictReader(f)
    cols = set(reader.fieldnames or [])
    if 'Run' in cols:
        id_col = 'Run'
    elif 'sample_id' in cols:
        id_col = 'sample_id'
    else:
        sys.exit(0)
    if 'collection_date' not in cols:
        sys.exit(0)
    rows = list(reader)

# Keep only rows with a valid YYYY-MM-DD date
valid_rows = []
for row in rows:
    try:
        datetime.strptime(row['collection_date'].strip(), '%Y-%m-%d')
        valid_rows.append(row)
    except ValueError:
        pass

# Keep only rows whose sample is in this aggregate's results
try:
    with open(agg_file) as f:
        f.readline()
        agg_samples = {line.split('\t')[0] for line in f if line.strip()}
    valid_rows = [r for r in valid_rows if r[id_col] + '.tsv' in agg_samples]
except FileNotFoundError:
    pass

if not valid_rows:
    sys.exit(0)

# Write times CSV
with open(times_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Sample', 'sample_collection_datetime'])
    for row in valid_rows:
        writer.writerow([row[id_col] + '.tsv', row['collection_date']])

# Write filtered aggregated TSV (freyja plot --times errors on unknown samples)
dated = {r[id_col] + '.tsv' for r in valid_rows}
with open(agg_file) as fin, open(filtered_agg, 'w') as fout:
    fout.write(fin.readline())
    for line in fin:
        if line.split('\t')[0] in dated:
            fout.write(line)

print(f"Time-series: {len(valid_rows)} sample(s) with valid dates")
PYEOF

        if [[ -f "${TIMES_CSV}" && -f "${FILTERED_AGG}" ]]; then
            freyja plot \
                "${FILTERED_AGG}" \
                --times "${TIMES_CSV}" \
                --interval MS \
                --output "${OUT_DIR}/lineage_timeseries.pdf" \
                --mincov 0
            echo "  [${LABEL}] lineage_timeseries.pdf ✓"
        fi
    fi

    echo "  [${LABEL}] → ${OUT_DIR}/"
}

# ── Run pipeline on current batch samples only ────────────────────────────────
if [[ -z "${ACC_FILE}" || ! -f "${ACC_FILE}" ]]; then
    echo "ERROR: No accession file provided or file not found: ${ACC_FILE}"
    exit 1
fi

BATCH_SAMPLES=()
while IFS= read -r acc; do
    acc=$(echo "${acc}" | sed 's/\r//' | tr -d '[:space:]')
    [[ -z "${acc}" || "${acc}" == \#* ]] && continue
    BATCH_SAMPLES+=("${acc}")
done < <({ cat "${ACC_FILE}"; echo; })

echo "Processing ${#BATCH_SAMPLES[@]} sample(s) from batch: ${BATCH_NAME}"

for SAMPLE in "${BATCH_SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[SKIP] ${SAMPLE} — FASTQ not found in ${RAW_DIR}"
        continue
    fi
    echo ""
    echo "════════════════════════════════════════════════════════════════════"
    echo " Processing: ${SAMPLE}"
    echo "════════════════════════════════════════════════════════════════════"
    bash "${SCRIPT_DIR}/run_sample.sh" "${SAMPLE}" "${R1}" "${R2}" "${RESULTS_DIR}"
done

# ── Per-batch aggregate ───────────────────────────────────────────────────────
echo ""
echo "[$(date)] Aggregating — batch: ${BATCH_NAME}"

BATCH_STAGING="${AGGREGATE_ROOT}/${BATCH_NAME}"
mkdir -p "${BATCH_STAGING}"
rm -f "${BATCH_STAGING}"/*.freyja.tsv

# Read accessions for this batch and stage only their freyja outputs
if [[ -n "${ACC_FILE}" && -f "${ACC_FILE}" ]]; then
    while IFS= read -r acc; do
        acc=$(echo "${acc}" | sed 's/\r//' | tr -d '[:space:]')
        [[ -z "${acc}" || "${acc}" == \#* ]] && continue
        TSV="${RESULTS_DIR}/${acc}/${acc}.freyja.tsv"
        if [[ -f "${TSV}" ]]; then
            cp "${TSV}" "${BATCH_STAGING}/"
        else
            echo "  [batch] WARNING: no result for ${acc} — not yet processed"
        fi
    done < <({ cat "${ACC_FILE}"; echo; })
else
    echo "  [batch] No accession file provided — skipping per-batch aggregate"
fi

aggregate_and_plot "${BATCH_STAGING}" "${BATCH_STAGING}" "${BATCH_NAME}"

# ── All-batches aggregate ─────────────────────────────────────────────────────
echo ""
echo "[$(date)] Aggregating — all batches"

ALL_STAGING="${AGGREGATE_ROOT}/all_batches"
mkdir -p "${ALL_STAGING}"
rm -f "${ALL_STAGING}"/*.freyja.tsv

find "${RESULTS_DIR}" -name "*.freyja.tsv" \
    -not -path "${AGGREGATE_ROOT}/*" \
    -exec cp {} "${ALL_STAGING}/" \;

aggregate_and_plot "${ALL_STAGING}" "${ALL_STAGING}" "all_batches"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "All done."
echo "  Per-batch results : ${AGGREGATE_ROOT}/${BATCH_NAME}/"
echo "  All-batches results: ${AGGREGATE_ROOT}/all_batches/"
