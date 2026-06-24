# COVID-19 Wastewater Sequencing Pipeline

Freyja-based pipeline for SARS-CoV-2 lineage abundance estimation from wastewater sequencing data.

**Tools:** Freyja · minimap2 · ivar · samtools · fastp · SRA toolkit  
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
│   ├── raw/                      ← FASTQs are downloaded here automatically
│   └── results/                  ← All pipeline outputs land here
│         ├── <SAMPLE>/               Per-sample folder
│         ├── all_samples_aggregated.tsv
│         ├── lineage_plot.pdf         Bar chart (all samples)
│         └── lineage_timeseries.pdf   Time-series (requires dates in metadata)
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
| 2 | Downloads any accessions not already in `data/raw/` (skips existing) |
| 3 | Picks the **newest** `.csv` in `data/metadata/` |
| 4 | Runs the Freyja pipeline on each new sample (skips already-processed ones) |
| 5 | Aggregates **all** results ever processed (cumulative) and produces plots |

### Pipeline steps per sample

```
fastp → minimap2 → samtools → ivar trim → freyja variants → freyja demix
  │         │          │           │              │                │
trim      align      sort &    remove         call            estimate
reads     to ref     index    primers        variants         lineages
```

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

### Metadata (drop into `data/metadata/`)

Download from NCBI SRA → select runs → **Metadata (RunInfo Table)**.

The pipeline uses these two columns (all others are ignored):

| Column | Example | Required for |
|--------|---------|-------------|
| `Run` | SRR39226689 | Always |
| `collection_date` | 2026-06-07 | Time-series plot only |

File can have any name — script picks the newest one automatically.

---

## Understanding the Output

### Per-sample result: `data/results/<SAMPLE>/<SAMPLE>.freyja.tsv`

| Field | What it means |
|-------|--------------|
| `summarized` | Lineage groups and their estimated % in the sample |
| `lineages` | Individual lineages detected |
| `abundances` | Estimated proportion of each lineage (sums to ~1) |
| `resid` | Fit quality — how well the mixture explains the data |
| `coverage` | % of genome covered by at least 1 read |

### Interpreting `resid` (residual)

The residual measures how well freyja's best lineage mixture explains the variant frequencies actually observed in your sample. Lower is better.

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

The script will only download and process new samples, then regenerate the aggregate plots with all samples (old + new).

---

## Primer Scheme Note

This pipeline trims primers using the **ARTIC v4.1** scheme. If your samples were sequenced with a different kit (e.g. **QIAseq DIRECT** from Qiagen), the primer trimming step will discard most reads, resulting in artificially low coverage (~25%) and a high `resid`.

**How to check:** look at the `amplicon_PCR_primer_scheme` column in your SRA metadata. If it says `ARTIC` you're fine. If it says something else, check the `resid` value — a `resid > 10` on a low-coverage sample is a sign of a primer scheme mismatch.

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

**Download seems stuck / silent for several minutes**  
→ It's probably gzipping large FASTQ files. This is silent and can take 5–10 minutes for large files. Check `data/raw/` to confirm files are appearing.

**Sample has very high resid (> 10)**  
→ Check the primer scheme in your metadata. If it's not ARTIC v4.1, the primer trimming may be wrong for those samples.

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
| sra-tools | 3.4.1 | SRA download |
