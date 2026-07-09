#!/usr/bin/env python3
"""Export freyja's individual-lineage -> summarized-constellation mapping
as a CSV — the actual `mapDict` freyja builds internally in
`buildLineageMap()` (freyja/sample_deconv.py) from curated_lineages.json,
just written out for reference instead of staying buried in a Python dict
during a `demix` run.

Replicates freyja's exact construction logic: records are sorted by their
number of pango_descendants, LARGEST first, then the map is built in that
order — so if a lineage is claimed by both a broad group and a narrower
one, the narrower (more specific) group's label wins, matching freyja's
own documented intent ("more specific designations will overwrite
broader, ancestral ones").

Simplifies `who_name` from freyja's `NAME* (NAME.X)` display format down to
just `NAME.X` — same transform bin/long_format.py applies to
`summarized_lineage` — so this file's summarized_lineage column matches
long_aggregated.csv's and the two can be joined directly.

Usage: build_lineage_map.py <curated_lineages.json> <output.csv>
"""
import csv
import json
import re
import sys

lineages_json, out_file = sys.argv[1:3]

PAREN_RE = re.compile(r'\(([^)]+)\)')


def simplify(name):
    m = PAREN_RE.search(name)
    return m.group(1) if m else name

with open(lineages_json) as f:
    dat = json.load(f)

# Same sort freyja itself uses — broadest groups (most descendants) first,
# so narrower groups processed later correctly override them below.
dat = sorted(dat, key=lambda x: len(x.get('pango_descendants', [])), reverse=True)

map_dict = {}
for record in dat:
    who_name = record.get('who_name')
    if who_name is None:
        continue
    for descendant in record.get('pango_descendants', []):
        map_dict[descendant] = simplify(who_name)

with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['lineage', 'summarized_lineage'])
    for lineage in sorted(map_dict):
        writer.writerow([lineage, map_dict[lineage]])

print(f"Wrote {len(map_dict)} lineage mappings ({out_file})")
