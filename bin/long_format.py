#!/usr/bin/env python3
"""Convert an aggregated.tsv (freyja aggregate output) into long format.

Each input row (one sample) becomes multiple output rows:
  - one row per entry in `summarized` (the grouped/constellation lineages)
  - one row per entry in `lineages` (the fine-grained individual calls)

These two sets are kept on SEPARATE rows (not paired/joined) — summarized
groups and individual lineages don't correspond 1-to-1 in freyja's output
(a group like XFG.X can contain several individual lineages), and the file
itself doesn't record which individual lineage belongs to which group. A
summarized-group row has lineage/lineage_abundance blank; an
individual-lineage row has summarized_lineage/summarized_lineage_abundance
blank. resid/coverage/site/collection_date repeat on every row for that
sample.

Also strips freyja's `NAME* (NAME.X)` display format down to just `NAME.X`
for both summarized and individual lineage names — entries with no
parentheses (e.g. "Other", "Omicron") are left as-is.

Usage: long_format.py <aggregated.tsv> <output.csv>
"""
import ast
import csv
import re
import sys

agg_file, out_file = sys.argv[1:3]

PAREN_RE = re.compile(r'\(([^)]+)\)')


def simplify(name):
    """'XFG* (XFG.X)' -> 'XFG.X'; 'Other' -> 'Other' (no parens, left as-is)."""
    m = PAREN_RE.search(name)
    return m.group(1) if m else name


rows_out = []

with open(agg_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        # First column has no header name in freyja's own aggregated.tsv
        # (the header row starts with a tab), so DictReader keys it as ''.
        # Strip the .tsv suffix for a clean ID.
        sample_id = row.get('', '').strip()
        sample_id = re.sub(r'\.tsv$', '', sample_id)

        resid = row.get('resid', '')
        coverage = row.get('coverage', '')
        site = row.get('site', '')
        collection_date = row.get('collection_date', '')

        # Summarized groups: a Python-literal list of (name, abundance) tuples
        summarized_raw = row.get('summarized', '').strip()
        if summarized_raw:
            try:
                summarized_list = ast.literal_eval(summarized_raw)
            except (ValueError, SyntaxError):
                summarized_list = []
            for name, abundance in summarized_list:
                rows_out.append({
                    'ID': sample_id,
                    'summarized_lineage': simplify(name),
                    'summarized_lineage_abundance': abundance,
                    'lineage': '',
                    'lineage_abundance': '',
                    'resid': resid,
                    'coverage': coverage,
                    'site': site,
                    'collection_date': collection_date,
                })

        # Individual lineages: space-separated parallel lists
        lineages = row.get('lineages', '').split()
        abundances = row.get('abundances', '').split()
        for name, abundance in zip(lineages, abundances):
            rows_out.append({
                'ID': sample_id,
                'summarized_lineage': '',
                'summarized_lineage_abundance': '',
                'lineage': simplify(name),
                'lineage_abundance': abundance,
                'resid': resid,
                'coverage': coverage,
                'site': site,
                'collection_date': collection_date,
            })

fieldnames = ['ID', 'summarized_lineage', 'summarized_lineage_abundance',
              'lineage', 'lineage_abundance', 'resid', 'coverage', 'site',
              'collection_date']

with open(out_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows_out)

print(f"Wrote {len(rows_out)} rows ({out_file})")
