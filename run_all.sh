#!/usr/bin/env bash
# Batch runner — processes all paired FASTQ files in data/raw/
#
# Expects files named:  <SAMPLE>_R1.fastq.gz  and  <SAMPLE>_R2.fastq.gz
#
# Usage: bash run_all.sh [RAW_DIR] [RESULTS_DIR]
#
# Example:
#   bash run_all.sh data/raw data/results

set -euo pipefail

RAW_DIR="${1:-data/raw}"
RESULTS_DIR="${2:-data/results}"
METADATA_FILE="${3:-data/metadata/metadata.csv}"

# ── Discover samples from R1 files ───────────────────────────────────────────
R1_FILES=()
while IFS= read -r f; do
    R1_FILES+=("$f")
done < <(find "${RAW_DIR}" -name "*_R1.fastq.gz" | sort)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "No *_R1.fastq.gz files found in ${RAW_DIR}"
    exit 1
fi

echo "Found ${#R1_FILES[@]} sample(s) to process"

# ── Process each sample ───────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for R1 in "${R1_FILES[@]}"; do
    SAMPLE=$(basename "${R1}" _R1.fastq.gz)
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "${R2}" ]]; then
        echo "[SKIP] ${SAMPLE} — R2 file not found: ${R2}"
        continue
    fi

    echo ""
    echo "════════════════════════════════════════════════════════════════════"
    echo " Processing: ${SAMPLE}"
    echo "════════════════════════════════════════════════════════════════════"
    bash "${SCRIPT_DIR}/run_sample.sh" "${SAMPLE}" "${R1}" "${R2}" "${RESULTS_DIR}"
done

# ── Aggregate all freyja outputs ──────────────────────────────────────────────
# freyja aggregate only looks in the top-level directory, so collect all
# per-sample freyja.tsv files into one flat staging directory first.
echo ""
echo "[$(date)] Aggregating all samples..."
AGGREGATE_DIR="${RESULTS_DIR}/_aggregate"
mkdir -p "${AGGREGATE_DIR}"
find "${RESULTS_DIR}" -name "*.freyja.tsv" -not -path "${AGGREGATE_DIR}/*" \
    -exec cp {} "${AGGREGATE_DIR}/" \;

freyja aggregate \
    "${AGGREGATE_DIR}/" \
    --output "${RESULTS_DIR}/all_samples_aggregated.tsv" \
    --ext freyja.tsv

# ── Summary plot (all samples) ────────────────────────────────────────────────
echo "[$(date)] Generating summary plot..."
freyja plot \
    "${RESULTS_DIR}/all_samples_aggregated.tsv" \
    --output "${RESULTS_DIR}/lineage_plot.pdf" \
    --mincov 0

# ── Time-series plot (only if metadata file exists) ───────────────────────────
if [[ -f "${METADATA_FILE}" ]]; then
    echo "[$(date)] Generating time-series plot..."

    # Build freyja-format times CSV from master metadata:
    #   master:  sample_id, collection_date, ...
    #   freyja:  Sample (= sample_id + ".tsv"), sample_collection_datetime
    TIMES_CSV="${RESULTS_DIR}/_times_metadata.csv"
    rm -f "${TIMES_CSV}"
    python3 - <<PYEOF
import csv, sys, os
from datetime import datetime

metadata_file = "${METADATA_FILE}"
out_file      = "${TIMES_CSV}"

with open(metadata_file) as f:
    reader = csv.DictReader(f)
    cols = set(reader.fieldnames or [])

    # Accept both SRA Run Table format (Run) and our custom format (sample_id)
    if 'Run' in cols:
        id_col = 'Run'
    elif 'sample_id' in cols:
        id_col = 'sample_id'
    else:
        print("WARNING: metadata has no 'Run' or 'sample_id' column. Skipping time-series plot.", file=sys.stderr)
        sys.exit(0)

    if 'collection_date' not in cols:
        print("WARNING: metadata has no 'collection_date' column. Skipping time-series plot.", file=sys.stderr)
        sys.exit(0)

    rows = list(reader)

valid_rows = []
skipped = []
for row in rows:
    date_str = row['collection_date'].strip()
    try:
        datetime.strptime(date_str, '%Y-%m-%d')
        valid_rows.append(row)
    except ValueError:
        skipped.append(row[id_col])

if skipped:
    print(f"WARNING: Skipping {len(skipped)} sample(s) with missing/placeholder dates: {skipped}", file=sys.stderr)

if not valid_rows:
    print("WARNING: No samples with valid collection dates — skipping time-series plot.", file=sys.stderr)
    print("         Fill in collection_date (YYYY-MM-DD) in data/metadata/metadata.csv to enable it.", file=sys.stderr)
    sys.exit(0)

# Only include samples that are actually in the aggregated results
agg_file = "${RESULTS_DIR}/all_samples_aggregated.tsv"
try:
    with open(agg_file) as f:
        header = f.readline()
        agg_samples = {line.split('\t')[0] for line in f if line.strip()}
except FileNotFoundError:
    agg_samples = None

if agg_samples is not None:
    before = len(valid_rows)
    valid_rows = [r for r in valid_rows if r[id_col] + '.tsv' in agg_samples]
    dropped = before - len(valid_rows)
    if dropped:
        print(f"INFO: {dropped} metadata row(s) skipped — not in aggregated results (no FASTQ or not yet processed).", file=sys.stderr)

if not valid_rows:
    print("WARNING: No metadata rows match processed samples — skipping time-series plot.", file=sys.stderr)
    sys.exit(0)

with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Sample', 'sample_collection_datetime'])
    for row in valid_rows:
        writer.writerow([row[id_col] + '.tsv', row['collection_date']])

print(f"Written {len(valid_rows)} rows to {out_file}")

# Write a filtered aggregated TSV containing only samples that have metadata.
# freyja plot --times errors if the aggregated file has samples not in times_df.
dated_samples = {r[id_col] + '.tsv' for r in valid_rows}
filtered_agg = agg_file.replace('.tsv', '_timed.tsv')
try:
    with open(agg_file) as fin, open(filtered_agg, 'w') as fout:
        header_line = fin.readline()
        fout.write(header_line)
        written = 0
        for line in fin:
            sample = line.split('\t')[0]
            if sample in dated_samples:
                fout.write(line)
                written += 1
    print(f"Filtered aggregated TSV: {written} samples → {filtered_agg}")
except FileNotFoundError:
    print(f"WARNING: Could not find {agg_file}", file=sys.stderr)
PYEOF

    FILTERED_AGG="${RESULTS_DIR}/all_samples_aggregated_timed.tsv"
    if [[ -f "${TIMES_CSV}" && -f "${FILTERED_AGG}" ]]; then
        freyja plot \
            "${FILTERED_AGG}" \
            --times "${TIMES_CSV}" \
            --interval MS \
            --output "${RESULTS_DIR}/lineage_timeseries.pdf" \
            --mincov 0
        echo "  Time-series plot : ${RESULTS_DIR}/lineage_timeseries.pdf"
    fi
else
    echo "[INFO] No metadata file found at ${METADATA_FILE} — skipping time-series plot."
    echo "       Create data/metadata/metadata.csv to enable it (see README or run_all.sh for format)."
fi

echo ""
echo "All done."
echo "  Aggregated table : ${RESULTS_DIR}/all_samples_aggregated.tsv"
echo "  Summary plot     : ${RESULTS_DIR}/lineage_plot.pdf"
