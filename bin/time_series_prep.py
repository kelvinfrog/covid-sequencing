#!/usr/bin/env python3
"""Build the sample->collection_date table + date-filtered aggregate needed
for `freyja plot --times`.

Usage: time_series_prep.py <metadata_file> <agg_file> <times_csv_out> <filtered_agg_out>

Exits 0 without creating the output files if no valid dated rows are found —
callers should check for the output files' existence before running
`freyja plot --times`.
"""
import csv
import sys
from datetime import datetime

metadata_file, agg_file, times_csv, filtered_agg = sys.argv[1:5]

with open(metadata_file) as f:
    reader = csv.DictReader(f)
    cols = set(reader.fieldnames or [])
    if 'Run' in cols:
        id_col = 'Run'
    elif 'sample_id' in cols:
        id_col = 'sample_id'
    else:
        sys.exit(0)
    # SraRunTable exports use 'Collection_Date'; older metadata.csv-style
    # files use lowercase 'collection_date' — accept either, same detection
    # pattern already used in annotate_site_date.py and primer_scheme_check.py.
    if 'Collection_Date' in cols:
        date_col = 'Collection_Date'
    elif 'collection_date' in cols:
        date_col = 'collection_date'
    else:
        sys.exit(0)
    rows = list(reader)

# Keep only rows with a valid YYYY-MM-DD date
valid_rows = []
for row in rows:
    try:
        datetime.strptime(row[date_col].strip(), '%Y-%m-%d')
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
        writer.writerow([row[id_col] + '.tsv', row[date_col]])

# Write filtered aggregated TSV (freyja plot --times errors on unknown samples)
dated = {r[id_col] + '.tsv' for r in valid_rows}
with open(agg_file) as fin, open(filtered_agg, 'w') as fout:
    fout.write(fin.readline())
    for line in fin:
        if line.split('\t')[0] in dated:
            fout.write(line)

print(f"Time-series: {len(valid_rows)} sample(s) with valid dates")
