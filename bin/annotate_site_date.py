#!/usr/bin/env python3
"""Add site + collection_date rows to a per-sample freyja.tsv.

Join: this accession's `Run` row in metadata -> ww_surv_system_sample_id ->
first 4 characters -> cdc_site_code.csv's sample_id -> site. Collection date
comes from the same metadata row's Collection_Date/collection_date column.

Adding two new row labels (rather than touching freyja's own five) matters
here: freyja's own `agg()` reads each file with pandas
`read_csv(..., index_col=0)` and concatenates on the index, so these rows
flow straight through into aggregated.tsv as new columns.

Rewrites the whole file rather than blindly appending — this must be
idempotent. Nextflow's `-resume` can re-execute this task against an
already-annotated file (e.g. a downstream script changes, forcing a rerun,
but this task's cached output was already annotated); appending again would
create duplicate 'site'/'collection_date' row labels, which breaks
`freyja aggregate`'s pandas concat with "Reindexing only valid with
uniquely valued Index objects". Stripping any pre-existing site/
collection_date rows before adding fresh ones avoids that regardless of how
many times this runs against the same file.

Usage: annotate_site_date.py <freyja_tsv> <accession> <cdc_site_code_csv> <metadata_csv>

Never fails the sample over a missing join — writes NO_MATCH and warns on
stderr instead, since one sample without cdc_site_code coverage shouldn't
take down the whole batch.
"""
import csv
import sys

freyja_tsv, accession, cdc_site_code_csv, metadata_csv = sys.argv[1:5]

site_by_id = {}
with open(cdc_site_code_csv, newline='') as f:
    for row in csv.DictReader(f):
        site_by_id[row['sample_id'].strip()] = row['site']

ww_id = None
collection_date = None
with open(metadata_csv, newline='') as f:
    reader = csv.DictReader(f)
    cols = reader.fieldnames or []
    date_col = 'Collection_Date' if 'Collection_Date' in cols else ('collection_date' if 'collection_date' in cols else None)
    for row in reader:
        if row.get('Run', '').strip() == accession:
            ww_id = row.get('ww_surv_system_sample_id', '').strip()
            collection_date = row.get(date_col, '').strip() if date_col else ''
            break

if not ww_id:
    print(f"  [annotate] WARNING: {accession} not found in metadata — site/collection_date set to NO_MATCH", file=sys.stderr)
    site = "NO_MATCH"
    collection_date = collection_date or "NO_MATCH"
else:
    prefix = ww_id[:4]
    site = site_by_id.get(prefix)
    if not site:
        print(f"  [annotate] WARNING: {accession} (ww_surv prefix {prefix}) has no match in cdc_site_code.csv", file=sys.stderr)
        site = "NO_MATCH"
    collection_date = collection_date or "NO_MATCH"

# Strip any pre-existing site/collection_date rows first, so re-running
# this against an already-annotated file (e.g. under -resume) is idempotent
# rather than appending duplicate row labels.
with open(freyja_tsv) as f:
    lines = [
        line for line in f
        if not (line.startswith('site\t') or line.startswith('collection_date\t'))
    ]

with open(freyja_tsv, 'w') as f:
    f.writelines(lines)
    if lines and not lines[-1].endswith('\n'):
        f.write('\n')
    f.write(f"site\t{site}\n")
    f.write(f"collection_date\t{collection_date}\n")

print(f"  [annotate] {accession}: site={site} collection_date={collection_date}")
