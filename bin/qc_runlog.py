#!/usr/bin/env python3
"""QC run log for a Freyja batch — ported from run_all.sh so the Nextflow
pipeline produces the same collaborator QC thresholds (breadth/10x/depth/resid).

Usage: qc_runlog.py <batch_name> <freyja_db_dir> <accessions_file> <out_file>

Expects <accession>.freyja.tsv and <accession>.depths.tsv for each accession
to be present in the current working directory (Nextflow stages them there).
"""
import sys
import os
from datetime import datetime

batch_name, freyja_db_dir, acc_file, out_file = sys.argv[1:5]

# Freyja barcode version — read FREYJA_UPDATE's own timestamp file (explicit,
# not inferred from a file mtime like the old bundled-copy approach).
barcode_path = os.path.join(freyja_db_dir, "last_barcode_update.txt")
try:
    with open(barcode_path) as f:
        barcode_version = f.read().strip()
except OSError:
    barcode_version = "unknown"

with open(acc_file) as f:
    accessions = [line.strip() for line in f if line.strip()]

# Quality thresholds (from collaborator, item 3 — same as run_all.sh)
MIN_COVERAGE_BREADTH = 60.0   # % of genome with any coverage
MIN_COVERAGE_10X     = 60.0   # % of genome at >=10x depth
MIN_MEAN_DEPTH       = 10.0   # mean read depth across genome


def parse_depths(depths_file):
    """Return (breadth_pct, breadth_10x_pct, mean_depth) from a freyja depths file.
    Depths file columns: chrom, pos, base, depth"""
    depths = []
    try:
        with open(depths_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    depths.append(int(parts[3]))
    except OSError:
        return None, None, None
    if not depths:
        return None, None, None
    total = len(depths)
    breadth = sum(1 for d in depths if d > 0) / total * 100
    breadth_10x = sum(1 for d in depths if d >= 10) / total * 100
    mean_depth = sum(depths) / total
    return breadth, breadth_10x, mean_depth


sample_results = []
for acc in accessions:
    tsv = f"{acc}.freyja.tsv"
    depths_tsv = f"{acc}.depths.tsv"
    result = {
        "sample": acc, "coverage": "n/a", "cov_10x": "n/a",
        "mean_depth": "n/a", "resid": "n/a", "status": "not processed"
    }
    if os.path.isfile(tsv):
        data = {}
        with open(tsv) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    data[parts[0]] = parts[1]
        cov_val = float(data.get('coverage', 0))
        resid_val = float(data.get('resid', 0))
        result["coverage"] = f"{cov_val:.1f}%"
        result["resid"] = f"{resid_val:.2f}"

        breadth, breadth_10x, mean_depth = parse_depths(depths_tsv)
        if breadth_10x is not None:
            result["cov_10x"] = f"{breadth_10x:.1f}%"
            result["mean_depth"] = f"{mean_depth:.1f}x"
        else:
            breadth_10x = cov_val  # fall back to breadth if depths missing
            mean_depth = 0.0

        fails = []
        if cov_val < MIN_COVERAGE_BREADTH:
            fails.append(f"breadth {cov_val:.0f}%<{MIN_COVERAGE_BREADTH:.0f}%")
        if breadth_10x < MIN_COVERAGE_10X:
            fails.append(f"10x {breadth_10x:.0f}%<{MIN_COVERAGE_10X:.0f}%")
        if mean_depth < MIN_MEAN_DEPTH:
            fails.append(f"depth {mean_depth:.1f}x<{MIN_MEAN_DEPTH:.0f}x")
        if resid_val > 10:
            fails.append(f"resid {resid_val:.1f}>10")

        result["status"] = "FAIL: " + ", ".join(fails) if fails else "PASS"
    sample_results.append(result)

with open(out_file, 'w') as f:
    f.write("=" * 78 + "\n")
    f.write(f"  RUN LOG — {batch_name}\n")
    f.write("=" * 78 + "\n")
    f.write(f"  Date           : {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"  Batch name     : {batch_name}\n")
    f.write("  Engine         : Nextflow\n")
    f.write(f"  Freyja barcodes: {barcode_version}\n")
    f.write("  Primer BED     : ARTIC_V5.3.2.bed\n")
    f.write(f"  Samples        : {len(accessions)}\n")
    f.write("\n")
    f.write(f"  QC thresholds  : breadth >={MIN_COVERAGE_BREADTH:.0f}%  |  10x breadth >={MIN_COVERAGE_10X:.0f}%  |  mean depth >={MIN_MEAN_DEPTH:.0f}x  |  resid <=10\n")
    f.write("\n")
    f.write("-" * 78 + "\n")
    f.write(f"  {'Sample':<18} {'Breadth':>8} {'10x Brd':>8} {'MeanDep':>8} {'Resid':>6}  Result\n")
    f.write("-" * 78 + "\n")
    for r in sample_results:
        f.write(f"  {r['sample']:<18} {r['coverage']:>8} {r['cov_10x']:>8} {r['mean_depth']:>8} {r['resid']:>6}  {r['status']}\n")
    f.write("-" * 78 + "\n")
    passed = [r for r in sample_results if r['status'] == 'PASS']
    failed = [r for r in sample_results if r['status'].startswith('FAIL')]
    f.write(f"\n  {len(passed)} PASS  |  {len(failed)} FAIL\n")
    if failed:
        f.write("\n  Samples below QC thresholds (do not report lineages):\n")
        for r in failed:
            f.write(f"    {r['sample']}: {r['status']}\n")
    f.write("\n")
    f.write(f"  Outputs: _aggregate/{batch_name}/\n")
    f.write("=" * 78 + "\n")

print(open(out_file).read())
