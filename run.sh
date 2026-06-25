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
#   2. Detect metadata and check primer schemes before any downloading
#   3. Download any accessions not yet in data/raw/
#   4. Run the Freyja pipeline on new samples (already-done ones are skipped)
#   5. Aggregate results and produce plots

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="data/raw"
RESULTS_DIR="data/results"
ACCESSIONS_DIR="data/accessions"
METADATA_DIR="data/metadata"
THREADS=4

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Pick newest accession file
# ─────────────────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════"
echo " Step 1 — Picking newest accession file"
echo "════════════════════════════════════════════════════════════════════"

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

BATCH_NAME="$(basename "${ACC_FILE}" .txt)"
echo "Using: $(basename "${ACC_FILE}")  →  batch name: ${BATCH_NAME}"

ACCESSIONS=()
while IFS= read -r acc; do
    [[ -z "$acc" ]] && continue
    ACCESSIONS+=("$acc")
done < <({ cat "${ACC_FILE}"; echo; } | sed 's/\r//' | grep -vE '^\s*(#|$)')

echo "Accessions in this batch: ${#ACCESSIONS[@]}"

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Detect metadata + check primer schemes
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════════"
echo " Step 2 — Detecting metadata & checking primer schemes"
echo "════════════════════════════════════════════════════════════════════"

METADATA_FILE=""
if [[ -d "${METADATA_DIR}" ]]; then
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
    echo ""

    # Build newline-separated accession list to pass into Python
    BATCH_ACC_LIST=$(printf '%s\n' "${ACCESSIONS[@]}")

    python3 - <<PYEOF
import csv, sys

metadata_file = "${METADATA_FILE}"
acc_set = set("""${BATCH_ACC_LIST}""".strip().split())
pipeline_scheme = "ARTIC V5.3.2"

primer_col = None
date_col   = None
scheme_by_acc = {}

with open(metadata_file) as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames or []
    for col in fields:
        cl = col.lower()
        if 'primer' in cl and 'scheme' in cl:
            primer_col = col
        if col in ('Collection_Date', 'collection_date'):
            date_col = col
    for row in reader:
        run_id = row.get('Run', '').strip()
        if run_id in acc_set:
            scheme = row.get(primer_col, 'not found').strip() if primer_col else 'not found'
            scheme_by_acc[run_id] = scheme

if not primer_col:
    print("  [primer check] No primer scheme column found in metadata — skipping check")
    sys.exit(0)

mismatches = []
for acc in sorted(acc_set):
    scheme = scheme_by_acc.get(acc, 'not in metadata')
    ok = pipeline_scheme.lower().replace(' ', '') in scheme.lower().replace(' ', '')
    status = "OK" if ok else "MISMATCH"
    print(f"  {acc}: {scheme}  [{status}]")
    if not ok:
        mismatches.append((acc, scheme))

BED_SOURCES = {
    "artic v4":    ("ARTIC v4.x",    "https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019/V4.1"),
    "artic v5":    ("ARTIC v5.x",    "already bundled — no download needed"),
    "qiaseq":      ("QIAseq DIRECT", "NOT publicly available — contact Qiagen for the BED file"),
    "midnight":    ("Midnight",       "https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019/Midnight-1200"),
}

def bed_hint(scheme):
    sl = scheme.lower().replace(' ', '')
    for key, (label, url) in BED_SOURCES.items():
        if key.replace(' ', '') in sl:
            return label, url
    return scheme, "Search GitHub for your primer scheme name + 'SARS-CoV-2 BED file'"

if mismatches:
    print()
    print("  *** WARNING: primer scheme mismatch detected ***")
    print(f"  Pipeline uses: {pipeline_scheme}")
    print()
    seen = {}
    for acc, scheme in mismatches:
        label, url = bed_hint(scheme)
        if scheme not in seen:
            seen[scheme] = (label, url)
        print(f"    {acc}: {scheme}")
    print()
    print("  To fix — download the correct BED file(s) and re-run:")
    for scheme, (label, url) in seen.items():
        print(f"    {scheme}:")
        print(f"      {url}")
        print(f"      Save as:  data/bed/<filename>.bed")
    print()
    print("  Then tell me the filename and I will update the pipeline to use it.")
    print("  After updating, delete the bad results and re-run:")
    print("    rm -rf data/results/<SAMPLE_ID>/")
    print("    bash run.sh")
else:
    print()
    print(f"  All samples match pipeline primer scheme ({pipeline_scheme}) — OK")
PYEOF

else
    echo "No .csv found in ${METADATA_DIR}/ — primer scheme check skipped."
    echo "Drop a metadata CSV there to enable it."
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Download missing samples
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════════"
echo " Step 3 — Downloading missing samples"
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
# Step 4: Run pipeline + aggregate + plot
# ─────────────────────────────────────────────────────────────────────────────
echo ""
bash "${SCRIPT_DIR}/run_all.sh" "${RAW_DIR}" "${RESULTS_DIR}" "${METADATA_FILE}" "${BATCH_NAME}" "${ACC_FILE}"
