#!/usr/bin/env python3
"""Pre-flight primer-scheme mismatch check against metadata. Ported from
run.sh Step 2 so Nextflow warns the same way bash does before processing
samples with the wrong primer BED. Advisory only — never blocks the run.

Usage: primer_scheme_check.py <metadata_file> <acc_file> <pipeline_scheme>
"""
import csv
import sys

metadata_file, acc_file, pipeline_scheme = sys.argv[1:4]

with open(acc_file) as f:
    acc_set = {
        line.strip() for line in f
        if line.strip() and not line.strip().startswith('#')
    }

primer_col = None
scheme_by_acc = {}

with open(metadata_file) as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames or []
    for col in fields:
        cl = col.lower()
        if 'primer' in cl and 'scheme' in cl:
            primer_col = col
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
    "artic v4": ("ARTIC v4.x", "https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019/V4.1"),
    "artic v5": ("ARTIC v5.x", "already bundled — no download needed"),
    "qiaseq": ("QIAseq DIRECT", "NOT publicly available — contact Qiagen for the BED file"),
    "midnight": ("Midnight", "https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019/Midnight-1200"),
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
        print("      Save as:  data/bed/<filename>.bed")
    print()
    print("  Then update the pipeline to use it, delete the bad results, and re-run:")
    print("    rm -rf data/results/<SAMPLE_ID>/")
    print("    nextflow run main.nf --acc_file <file> --batch_name <name>")
else:
    print()
    print(f"  All samples match pipeline primer scheme ({pipeline_scheme}) — OK")
