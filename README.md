# COVID-19 Wastewater Sequencing Pipeline

Freyja-based pipeline for SARS-CoV-2 lineage abundance estimation from wastewater sequencing data.

**Tools:** Freyja · minimap2 · ivar · samtools · fastp · seqtk · SRA toolkit  
**Conda environment:** `freyja-env`

---

## Quick Start

1. Drop your SRA accession list (`.txt`) into `data/accessions/`
2. Drop your SRA metadata table (`.csv`) into `data/metadata/`
3. Run:

```bash
conda activate freyja-env
cd ~/Library/Mobile\ Documents/com~apple~CloudDocs/covid_sequencing
bash run.sh
```

That's it. The script handles everything else automatically.

---

## Directory Structure

```
covid_sequencing/
│
├── run.sh                        ← MASTER SCRIPT — the only one you need to run
│
├── data/
│   ├── accessions/               ← DROP your SRA accession .txt files here
│   ├── metadata/                 ← DROP your SRA metadata .csv files here
│   ├── bed/                      ← Primer scheme BED files (ARTIC v5.3.2 bundled)
│   ├── raw/                      ← FASTQs are downloaded here automatically
│   └── results/                  ← All pipeline outputs land here
│         ├── <SAMPLE>/               Per-sample folder
│         └── _aggregate/
│               ├── <batch_name>/     Per-batch aggregate (current batch only)
│               │     aggregated.tsv
│               │     lineage_plot.pdf
│               │     lineage_timeseries.pdf
│               └── all_batches/      Cumulative aggregate (every sample ever run)
│                     aggregated.tsv
│                     lineage_plot.pdf
│                     lineage_timeseries.pdf
│
├── run_sample.sh                 ← Single-sample pipeline (used internally)
├── run_all.sh                    ← Batch runner + aggregate + plot (used internally)
├── download.sh                   ← Standalone downloader (optional)
└── README.md                     ← This file
```

---

## What the Script Does

| Step | What happens |
|------|-------------|
| 1 | Picks the **newest** `.txt` in `data/accessions/` |
| 2 | Detects metadata and **checks primer schemes** before downloading — warns if mismatch |
| 3 | Downloads any accessions not already in `data/raw/` (skips existing) |
| 4 | Runs the Freyja pipeline on each new sample (skips already-processed ones) |
| 5 | Produces a **per-batch** aggregate (current batch only) |
| 6 | Produces an **all-batches** aggregate (every sample ever processed) |

### Pipeline steps per sample

```
subsample → fastp → minimap2 → samtools → ivar trim → freyja variants → freyja demix
    │          │         │          │           │              │                │
cap at      trim      align      sort &    remove         call            estimate
3M reads   reads     to ref     index    primers        variants         lineages
```

> Subsampling only happens if the sample has more than 3 million reads. Most wastewater samples need far fewer reads than labs sequence — capping at 3M reduces processing time significantly without affecting lineage calls.

---

## Input File Formats

### Accession list (drop into `data/accessions/`)

Plain text, one SRA accession per line. Download from NCBI SRA → select runs → **Accession List**.

```
SRR39226689
SRR39226698
SRR23879090
```

- Lines starting with `#` are ignored
- File can have any name — script picks the newest one automatically
- Each new file creates its own batch folder in `_aggregate/`

### Metadata (drop into `data/metadata/`)

Download from NCBI SRA → select runs → **Metadata (RunInfo Table)**.

The pipeline uses these columns (all others are ignored):

| Column | Example | Required for |
|--------|---------|-------------|
| `Run` | SRR39226689 | Always |
| `amplicon_PCR_primer_scheme` | ARTIC V5.3.2 | Primer scheme check |
| `collection_date` | 2026-06-07 | Time-series plot only |

File can have any name — script picks the newest one automatically.

---

## Understanding the Output

### Per-batch vs all-batches

Every run produces two sets of aggregate outputs:

- `_aggregate/<batch_name>/` — only the samples from the current accession file
- `_aggregate/all_batches/` — every sample ever processed, updated automatically

### Per-sample result: `data/results/<SAMPLE>/<SAMPLE>.freyja.tsv`

| Field | What it means |
|-------|--------------|
| `summarized` | Lineage groups and their estimated % in the sample |
| `lineages` | Individual lineages detected |
| `abundances` | Estimated proportion of each lineage (sums to ~1) |
| `resid` | Fit quality — how well the mixture explains the data |
| `coverage` | % of genome covered by at least 1 read |

### Interpreting `resid` (residual)

The residual measures how well freyja's best lineage mixture explains the variant frequencies observed. Lower is better.

| resid | Interpretation |
|-------|---------------|
| < 2 | Excellent — lineage calls are very reliable |
| 2–5 | Good — results are trustworthy |
| 5–10 | Borderline — treat with some caution |
| > 10 | Poor — possible data quality issue (wrong primers, contamination) |

### Interpreting coverage

The pipeline uses `--mincov 0` so **all samples appear in plots** regardless of coverage. For wastewater, **20–50% coverage is normal and expected** — RNA in wastewater is degraded and viral concentrations are low. The 60% freyja default was designed for clinical samples.

> Use `resid` as the primary quality indicator, not coverage, for wastewater data.

---

## Primer Scheme Check

Before downloading anything, the pipeline checks the `amplicon_PCR_primer_scheme` column in your metadata and warns if any sample uses a different primer scheme than the pipeline default (**ARTIC v5.3.2**).

Example warning output:
```
  SRR12345678: QIAseq DIRECT  [MISMATCH]

  *** WARNING: primer scheme mismatch detected ***
  Pipeline uses: ARTIC V5.3.2
  Affected samples:
    SRR12345678 uses QIAseq DIRECT
  ...
```

**If you see a mismatch:**
1. Get the correct BED file for that primer scheme
2. Drop it into `data/bed/`
3. Tell the pipeline maintainer the filename — they will update the pipeline to use it
4. Delete the bad results and re-run:
```bash
rm -rf data/results/<SAMPLE_ID>/
bash run.sh
```

**Known primer schemes and BED file sources:**

| Scheme | Source |
|--------|--------|
| ARTIC v5.3.2 | Bundled — no action needed |
| ARTIC v4.x | https://github.com/artic-network/primer-schemes |
| Midnight | https://github.com/artic-network/primer-schemes |
| QIAseq DIRECT | Not public — contact Qiagen |

---

## Routine Workflow (adding new batches)

Each time you have a new batch of samples:

1. Go to NCBI SRA, find your runs, download:
   - **Accession List** → save as `.txt` → drop in `data/accessions/`
   - **RunInfo Table** → save as `.csv` → drop in `data/metadata/`

2. Run:
```bash
conda activate freyja-env
bash run.sh
```

The script will check primer schemes, download only new samples, process only unprocessed ones, and update both the per-batch and all-batches aggregate plots.

---

## Keeping Freyja Up to Date

Freyja's lineage barcode database is updated as new variants emerge. Run this periodically (e.g. once a month):

```bash
conda activate freyja-env
freyja update
```

---

## Troubleshooting

**"No .txt files found in data/accessions/"**  
→ Drop an accession list `.txt` file into `data/accessions/` first.

**Primer scheme mismatch warning**  
→ The pipeline detected that your samples use a different primer scheme than ARTIC v5.3.2. Get the correct BED file from the sequencing lab or the source listed in the warning, drop it in `data/bed/`, and update the pipeline before processing.

**Sample has very high resid (> 10)**  
→ Most likely a primer scheme mismatch. Check the `amplicon_PCR_primer_scheme` column in your metadata. The pipeline will warn you automatically if a mismatch is detected before the run starts.

**Time-series plot not generated**  
→ Your metadata `.csv` must have a `collection_date` column in `YYYY-MM-DD` format. Samples without a valid date are excluded from the time-series.

**"freyja: command not found"**  
→ Activate the conda environment first: `conda activate freyja-env`

**Want to re-process a sample from scratch**  
→ Delete its results folder and re-run:
```bash
rm -rf data/results/<SAMPLE_ID>/
bash run.sh
```

---

## Tool Versions

| Tool | Version | Purpose |
|------|---------|---------|
| freyja | 2.0.3 | Lineage demixing |
| ivar | 1.4.4 | Primer trimming |
| minimap2 | 2.31 | Read alignment |
| samtools | 1.21 | BAM processing |
| fastp | 1.1.0 | Read quality trimming |
| seqtk | 1.5-r133 | Read subsampling |
| sra-tools | 3.4.1 | SRA download |
