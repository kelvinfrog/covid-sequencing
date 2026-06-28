# COVID-19 Wastewater Sequencing Pipeline

Freyja-based pipeline for SARS-CoV-2 lineage abundance estimation from wastewater sequencing data. Supports both a bash workflow and a Nextflow workflow.

**Tools:** Freyja · minimap2 · ivar · samtools · fastp · SRA Toolkit · Nextflow  
**Conda environment:** `freyja-env`

---

## Quick Start

1. Drop your SRA accession list (`.txt`) into `data/accessions/`
2. Drop your SRA metadata table (`.csv`) into `data/metadata/`
3. Run:

**Bash (sequential):**
```bash
conda activate freyja-env
bash run.sh
```

**Nextflow (parallel — ~3× faster):**
```bash
conda activate freyja-env
nextflow run main.nf --acc_file data/accessions/<your_file>.txt --batch_name <batch_name>
```

---

## Directory Structure

```
covid_sequencing/
│
├── run.sh                        ← Bash master script — entry point
├── run_all.sh                    ← Batch runner + aggregate + run log
├── run_sample.sh                 ← Single-sample pipeline (called internally)
│
├── main.nf                       ← Nextflow pipeline (parallel)
├── nextflow.config               ← Nextflow runtime settings
│
├── data/
│   ├── accessions/               ← DROP your SRA accession .txt files here
│   ├── metadata/                 ← DROP your SRA metadata .csv files here
│   ├── bed/                      ← Primer scheme BED files (ARTIC v5.3.2 bundled)
│   ├── raw/                      ← FASTQs downloaded here automatically
│   └── results/
│         ├── <SAMPLE>/           ← Per-sample outputs (freyja.tsv, logs)
│         └── _aggregate/
│               ├── <batch_name>/ ← Per-batch aggregate (current batch only)
│               │     aggregated.tsv
│               │     lineage_plot.pdf
│               │     lineage_timeseries.pdf
│               │     run_log.txt         ← Run log with barcode version + QC
│               └── all_batches/  ← Cumulative aggregate (all samples ever run)
│                     aggregated.tsv
│                     lineage_plot.pdf
│                     lineage_timeseries.pdf
│
└── CLAUDE.md                     ← Domain context for AI assistant
```

---

## What the Pipeline Does

| Step | What happens |
|------|-------------|
| 1 | Picks the **newest** `.txt` in `data/accessions/` as the current batch |
| 2 | Detects metadata and **checks primer schemes** before downloading |
| 3 | Downloads FASTQs from SRA — first **2 million reads** per sample |
| 4 | Runs the Freyja pipeline per sample (skips already-processed ones) |
| 5 | Produces a **per-batch** aggregate and plots |
| 6 | Produces an **all-batches** aggregate and plots |
| 7 | Writes a **run log** with barcode version, coverage, and resid per sample |

### Per-sample pipeline steps

```
fastq-dump → fastp → minimap2 → samtools sort/index → ivar trim → freyja variants → freyja demix
     │          │         │              │                  │              │                │
download     quality    align          sort &           remove          call           estimate
first 2M     trim       to ref         index           primers        variants         lineages
reads
```

> Downloads are capped at 2 million reads per sample (`--maxSpotId 2000000`). For wastewater target-capture sequencing, 2M reads is sufficient for reliable lineage calls and keeps file sizes manageable (~200–350 MB per sample vs 1–2 GB for full downloads).

---

## Bash vs Nextflow

Both workflows produce identical results. Choose based on your needs:

| | Bash (`run.sh`) | Nextflow (`main.nf`) |
|--|--|--|
| **Parallelism** | Sequential (one sample at a time) | Parallel (2 samples at a time by default) |
| **Speed (10 samples)** | ~2 hours | ~40 minutes |
| **Resume failed runs** | Skip guard (re-run `run.sh`) | `-resume` flag caches completed steps |
| **Best for** | Simplicity, debugging | Speed, larger batches |

**To increase Nextflow parallelism** on a Linux server with more cores, edit `nextflow.config`:
```groovy
maxForks = 8   // run 8 samples at a time
```

---

## Input File Formats

### Accession list (`data/accessions/`)

Plain text, one SRA accession per line. Download from NCBI SRA → select runs → **Accession List**.

```
SRR39226689
SRR39226698
SRR23879090
```

- Lines starting with `#` are ignored
- The script always picks the **newest** file automatically
- Each file creates its own batch folder in `_aggregate/`

### Metadata (`data/metadata/`)

Download from NCBI SRA → select runs → **Metadata (RunInfo Table)**.

The pipeline uses these columns:

| Column | Example | Used for |
|--------|---------|---------|
| `Run` | SRR39226689 | Matching accessions |
| `amplicon_PCR_primer_scheme` | ARTIC V5.3.2 | Primer scheme check |
| `collection_date` | 2026-06-07 | Time-series plot |

The script always picks the **newest** `.csv` automatically.

---

## Understanding the Output

### Run log: `_aggregate/<batch_name>/run_log.txt`

Generated automatically after every run. Records:

- Date and batch name
- Freyja barcode database version (date of `curated_lineages.json`)
- Primer BED file used
- Per-sample coverage and resid
- Any samples flagged for quality issues

Example:
```
====================================================================
  RUN LOG — SRR_Acc_List_8
====================================================================
  Date           : 2026-06-28 16:25
  Freyja barcodes: 2026-06-24
  Primer BED     : ARTIC_V5.3.2.bed
  Samples        : 10

  Sample               Coverage    Resid  Status
  SRR39224810             45.9%     1.51  OK
  SRR39239887             52.6%    10.23  WARNING: high resid
====================================================================
```

> The barcode version matters: re-running old samples with a new barcode database can change lineage calls. The run log makes week-over-week comparisons traceable.

### Per-sample result: `data/results/<SAMPLE>/<SAMPLE>.freyja.tsv`

| Field | What it means |
|-------|--------------|
| `summarized` | Lineage groups and their estimated proportion |
| `lineages` | Individual lineages detected |
| `abundances` | Estimated proportion of each lineage (sums to ~1) |
| `resid` | Fit quality — lower is better |
| `coverage` | % of genome covered by at least 1 read |

### Interpreting `resid` (residual)

| resid | Interpretation |
|-------|---------------|
| < 2 | Excellent |
| 2–5 | Good |
| 5–10 | Borderline — treat with some caution |
| > 10 | Poor — likely a primer scheme mismatch or data quality issue |

### Interpreting coverage

The pipeline uses `--mincov 0` so all samples appear in output regardless of coverage. For wastewater, **20–50% coverage is normal** — viral RNA is degraded and concentrations are low. The freyja default of 60% was designed for clinical samples.

> Use `resid` as the primary quality indicator, not coverage, for wastewater data.

---

## Primer Scheme Check

Before downloading, the pipeline checks the `amplicon_PCR_primer_scheme` column in your metadata. It warns if any sample uses a scheme other than the pipeline default (**ARTIC v5.3.2**).

**If you see a mismatch:**
1. Get the correct BED file for that primer scheme
2. Drop it into `data/bed/`
3. Update the pipeline to use the new BED file
4. Delete the affected results and re-run:
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

## Keeping Freyja Up to Date

Freyja's lineage barcode database updates as new variants emerge. Run this periodically (e.g. once a month, or before a new batch):

```bash
conda activate freyja-env
freyja update
```

The run log records which barcode version was used for each batch, so you can always trace what database produced a given result.

---

## Routine Workflow (adding new batches)

1. Go to NCBI SRA, find your runs, download:
   - **Accession List** → save as `.txt` → drop in `data/accessions/`
   - **RunInfo Table** → save as `.csv` → drop in `data/metadata/`

2. Run:
```bash
conda activate freyja-env
bash run.sh
# or
nextflow run main.nf --acc_file data/accessions/<file>.txt --batch_name <name>
```

The pipeline checks primer schemes, downloads only new samples, processes only unprocessed ones, and updates both per-batch and all-batches aggregates.

---

## Troubleshooting

**"No .txt files found in data/accessions/"**  
→ Drop an accession list `.txt` into `data/accessions/` first.

**Primer scheme mismatch warning**  
→ Get the correct BED file from the source listed in the warning, drop it in `data/bed/`, and update the pipeline before re-processing.

**Sample flagged with high resid (> 10)**  
→ Most likely a primer scheme mismatch. Check the `amplicon_PCR_primer_scheme` column in your metadata.

**Time-series plot not generated**  
→ Your metadata `.csv` must have a `collection_date` column in `YYYY-MM-DD` format.

**"freyja: command not found"**  
→ Activate the conda environment: `conda activate freyja-env`

**Re-process a sample from scratch:**
```bash
rm -rf data/results/<SAMPLE_ID>/
bash run.sh
```

**Nextflow resume (skip already-completed steps):**
```bash
nextflow run main.nf --acc_file data/accessions/<file>.txt --batch_name <name> -resume
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
| sra-tools | 3.4.1 | SRA download (2M read cap) |
| nextflow | 24.x | Parallel workflow manager |
