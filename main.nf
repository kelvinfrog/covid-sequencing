#!/usr/bin/env nextflow
// ─────────────────────────────────────────────────────────────────────────────
// main.nf — Nextflow version of the COVID wastewater sequencing pipeline
//
// Run with:
//   nextflow run main.nf --acc_file data/accessions/your_batch.txt \
//                        --batch_name your_batch_name
//
// Resume a failed run (skips already-completed steps automatically):
//   nextflow run main.nf --acc_file ... --batch_name ... -resume
//
// This mirrors run.sh + run_all.sh + run_sample.sh feature-for-feature:
// primer-scheme pre-check, per-sample QC run log, per-batch aggregate,
// cross-run "all_batches" aggregate, and the metadata-driven time-series
// plot. See bin/*.py for the ported QC/reporting logic.
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMETERS — override any of these from the command line
// e.g. --threads 8  or  --max_reads 3000000
//
// Note: `?:` is used everywhere below so a value passed on the CLI (e.g.
// --ref /some/path) is not silently clobbered by the default assignment.
// ─────────────────────────────────────────────────────────────────────────────
params.acc_file      = params.acc_file      ?: ""
params.raw_dir       = params.raw_dir       ?: "data/raw"
params.results_dir   = params.results_dir   ?: "data/results"
params.metadata_dir  = params.metadata_dir  ?: "data/metadata"
params.cdc_site_code = params.cdc_site_code ?: "${projectDir}/cdc_site_code.csv"
params.bed           = params.bed           ?: "${projectDir}/data/bed/ARTIC_V5.3.2.bed"
params.max_reads     = params.max_reads     ?: 2000000
params.threads       = params.threads       ?: 4
// minq for freyja variants' ivar-variants call only (CDPHE: -q 20). NOT used
// for fastp anymore — see trim_min_quality/trim_min_length below, which are
// deliberately decoupled since CDPHE itself uses different thresholds for
// read-trimming (seqyclean -qual 30 30 -minlen 70) vs. variant calling.
params.min_quality      = params.min_quality      ?: 20
params.trim_min_quality = params.trim_min_quality ?: 30
params.trim_min_length  = params.trim_min_length  ?: 70
params.batch_name    = params.batch_name    ?: "nextflow_batch"
// `ref`, `gff`, and `metadata_file` are deliberately NOT declared here. This
// Nextflow version treats `params.X = ...` as a one-time declaration — the
// FIRST assignment wins and any later one (e.g. inside workflow{}, where the
// dynamic resolution actually needs to run) is silently discarded. All three
// are assigned exactly once, inside workflow{} below.

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: PRIMER_CHECK
//
// Advisory pre-flight check — warns if metadata's primer scheme column
// doesn't match the pipeline's ARTIC v5.3.2 BED. Never blocks the run,
// same as run.sh Step 2.
// ─────────────────────────────────────────────────────────────────────────────
process PRIMER_CHECK {
    input:
    path metadata
    path acc_list

    script:
    """
    primer_scheme_check.py ${metadata} ${acc_list} "ARTIC V5.3.2"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: DOWNLOAD
//
// Downloads the first N reads from SRA using fastq-dump.
// The `when` block skips the download if BOTH R1 and R2 already exist in
// data/raw/ — matching run.sh's [EXISTS] check (checking only R1 would
// treat a partial/failed prior download as complete).
// ─────────────────────────────────────────────────────────────────────────────
process DOWNLOAD {
    publishDir params.raw_dir, mode: 'copy'

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_R1.fastq.gz"), path("${accession}_R2.fastq.gz")

    when:
    !(file("${params.raw_dir}/${accession}_R1.fastq.gz").exists() &&
      file("${params.raw_dir}/${accession}_R2.fastq.gz").exists())

    script:
    """
    fastq-dump ${accession} \
        --maxSpotId ${params.max_reads} \
        --split-files \
        --gzip \
        --outdir .
    mv ${accession}_1.fastq.gz ${accession}_R1.fastq.gz
    mv ${accession}_2.fastq.gz ${accession}_R2.fastq.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FASTP
//
// Quality trims reads and removes adapters.
//
// Only the small QC reports (json/html) are published to results_dir — the
// trimmed FASTQ itself is a large intermediate that only ALIGN needs next;
// keeping a permanent copy doubles disk usage for no downstream benefit.
// It's still staged for ALIGN via Nextflow's own work directory as normal.
// ─────────────────────────────────────────────────────────────────────────────
process FASTP {
    publishDir { "${params.results_dir}/${accession}/trimmed" }, mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple val(accession), path(r1), path(r2)

    output:
    tuple val(accession), path("${accession}_trimmed_R1.fastq.gz"), path("${accession}_trimmed_R2.fastq.gz"), emit: reads
    path "${accession}_fastp.json", emit: json
    path "${accession}_fastp.html", emit: html

    script:
    """
    fastp \
        --in1 ${r1} --in2 ${r2} \
        --out1 ${accession}_trimmed_R1.fastq.gz \
        --out2 ${accession}_trimmed_R2.fastq.gz \
        --json ${accession}_fastp.json \
        --html ${accession}_fastp.html \
        --thread ${params.threads} \
        --qualified_quality_phred ${params.trim_min_quality} \
        --length_required ${params.trim_min_length} \
        --detect_adapter_for_pe
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: ALIGN
//
// Aligns trimmed reads to SARS-CoV-2 reference, then sorts and indexes BAM.
// ─────────────────────────────────────────────────────────────────────────────
process ALIGN {
    publishDir { "${params.results_dir}/${accession}/aligned" }, mode: 'copy'

    input:
    tuple val(accession), path(r1), path(r2)

    output:
    tuple val(accession), path("${accession}.sorted.bam"), path("${accession}.sorted.bam.bai")

    script:
    """
    minimap2 -ax sr -t ${params.threads} ${params.ref} ${r1} ${r2} \
        | samtools sort -@ ${params.threads} -o ${accession}.sorted.bam
    samtools index ${accession}.sorted.bam
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: IVAR_TRIM
//
// Removes primer sequences using the ARTIC v5.3.2 BED file, then re-sorts.
// `-e` (keep reads with no primer found) only, matching CDPHE's
// SC2_illumina_pe_assembly ivar_trim task exactly — this is confirmed
// tiled-amplicon data (same assay type CDPHE's pipeline is built for), so we
// align to their convention rather than overriding -q/-m/-s off ivar's
// defaults (q20/m30/s4).
// ─────────────────────────────────────────────────────────────────────────────
process IVAR_TRIM {
    publishDir { "${params.results_dir}/${accession}/aligned" }, mode: 'copy'

    input:
    tuple val(accession), path(bam), path(bai)
    path bed

    output:
    tuple val(accession), path("${accession}.trimmed.sorted.bam"), path("${accession}.trimmed.sorted.bam.bai")

    script:
    """
    ivar trim -e -i ${bam} -b ${bed} -p ${accession}.trimmed
    samtools sort -@ ${params.threads} ${accession}.trimmed.bam \
        -o ${accession}.trimmed.sorted.bam
    samtools index ${accession}.trimmed.sorted.bam
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_VARIANTS
//
// Calls variants and measures depth across the genome. --annot matches
// CDPHE's `ivar variants -g {gff}` — adds AA/codon columns to variants.tsv;
// doesn't change depths.tsv or downstream demix/lineage results.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_VARIANTS {
    publishDir { "${params.results_dir}/${accession}/variants" }, mode: 'copy'

    input:
    tuple val(accession), path(bam), path(bai)

    output:
    tuple val(accession), path("${accession}.tsv"), path("${accession}.depths.tsv")

    script:
    """
    freyja variants ${bam} \
        --variants ${accession} \
        --depths ${accession}.depths.tsv \
        --ref ${params.ref} \
        --annot ${params.gff} \
        --minq ${params.min_quality}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_UPDATE
//
// Downloads fresh UShER barcodes + curated lineages ONCE per pipeline
// invocation (not once per sample, unlike CDPHE's per-sample `freyja update`
// — this avoids N redundant downloads and guarantees every sample in the
// batch demixes against the identical barcode version). ~5s, ~16MB; adds a
// real network dependency to every run in exchange for always-current
// lineage calls instead of whatever was bundled at freyja-env install time.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_UPDATE {
    output:
    path "freyja_db"

    script:
    """
    mkdir -p freyja_db
    freyja update --outdir freyja_db
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_DEMIX
//
// Estimates lineage abundances — the core freyja step. --barcodes/--meta
// point at FREYJA_UPDATE's fresh download rather than freyja-env's bundled
// (install-time) copy.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_DEMIX {
    publishDir { "${params.results_dir}/${accession}" }, mode: 'copy'

    input:
    tuple val(accession), path(variants), path(depths)
    path freyja_db

    output:
    tuple val(accession), path("${accession}.freyja.tsv")

    script:
    """
    freyja demix ${variants} ${depths} \
        --output ${accession}.freyja.tsv \
        --confirmedonly \
        --eps 0.01 \
        --covcut 10 \
        --barcodes ${freyja_db}/usher_barcodes.feather \
        --meta ${freyja_db}/curated_lineages.json
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: ANNOTATE_SITE_DATE
//
// Appends `site` + `collection_date` rows to the per-sample freyja.tsv, via
// the cdc_site_code.csv join: this accession's Run row in metadata ->
// ww_surv_system_sample_id -> first 4 digits -> cdc_site_code.csv's
// sample_id -> site. Collection date comes from the same metadata row.
//
// Appending (not rewriting) is deliberate: freyja's own aggregate command
// reads each file with pandas index_col=0 and concatenates on the index, so
// these extra rows flow straight through into aggregated.tsv as new columns
// without needing to touch freyja's own five rows. See bin/annotate_site_date.py.
// ─────────────────────────────────────────────────────────────────────────────
process ANNOTATE_SITE_DATE {
    publishDir { "${params.results_dir}/${accession}" }, mode: 'copy'

    input:
    tuple val(accession), path(freyja_tsv)
    path cdc_site_code
    path metadata

    output:
    tuple val(accession), path(freyja_tsv)

    script:
    """
    annotate_site_date.py ${freyja_tsv} ${accession} ${cdc_site_code} ${metadata}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: CLEANUP_INTERMEDIATES
//
// Deletes the raw FASTQ and BAM files for a sample once its freyja.tsv is
// confirmed generated — same reasoning as the trimmed-FASTQ cleanup in
// FASTP, one stage later. At 400 samples/month, raw+BAM footprint alone
// (~614MB/sample measured) would exceed a 200GB volume within a single
// month; the actual valuable output (freyja.tsv + aggregates) is only a
// few KB/sample. BAM reuse for automated reanalysis is rare in practice —
// Nextflow's own -resume cache (not these published copies) is what
// actually drives reanalysis — and on-demand regeneration for the
// occasional manual investigation (e.g. a QC-flagged sample) is cheap
// (~13 min/sample), so there's little reason to keep these long-term.
//
// Takes the demix result as input (not just the accession) purely for
// dependency ordering — this must not run until freyja.tsv actually exists.
// Deletes at the real published paths (rawDirAbs/resultsDirAbs), not
// relative ones — same reason as FREYJA_AGGREGATE_ALL's `find`: process
// scripts execute in an isolated work directory, so a relative path here
// would silently resolve wrong.
// ─────────────────────────────────────────────────────────────────────────────
process CLEANUP_INTERMEDIATES {
    input:
    tuple val(accession), path(freyja_tsv)
    val rawDirAbs
    val resultsDirAbs

    script:
    """
    rm -f "${rawDirAbs}/${accession}_R1.fastq.gz" "${rawDirAbs}/${accession}_R2.fastq.gz"
    rm -f "${resultsDirAbs}/${accession}/aligned/${accession}.sorted.bam" \
          "${resultsDirAbs}/${accession}/aligned/${accession}.sorted.bam.bai" \
          "${resultsDirAbs}/${accession}/aligned/${accession}.trimmed.sorted.bam" \
          "${resultsDirAbs}/${accession}/aligned/${accession}.trimmed.sorted.bam.bai"
    echo "[cleanup] ${accession}: removed raw FASTQ + BAMs"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: QC_RUNLOG
//
// Per-sample QC PASS/FAIL grid (breadth/10x-breadth/mean-depth/resid) —
// ported from run_all.sh's collaborator QC thresholds. See bin/qc_runlog.py.
// ─────────────────────────────────────────────────────────────────────────────
process QC_RUNLOG {
    publishDir "${params.results_dir}/_aggregate/${params.batch_name}", mode: 'copy'

    input:
    val accessions
    path freyja_tsvs
    path depths_tsvs
    path freyja_db

    output:
    path "run_log.txt"

    script:
    def accBlock = accessions.join('\n')
    """
    cat > accessions.txt <<'ACCLIST'
${accBlock}
ACCLIST
    qc_runlog.py "${params.batch_name}" ${freyja_db} accessions.txt run_log.txt
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_AGGREGATE
//
// Combines this batch's freyja.tsv files into one aggregated TSV.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_AGGREGATE {
    publishDir "${params.results_dir}/_aggregate/${params.batch_name}", mode: 'copy'

    input:
    path tsvs

    output:
    path "aggregated.tsv"

    script:
    """
    mkdir -p staging
    for f in ${tsvs}; do cp \$f staging/; done
    freyja aggregate staging/ --output aggregated.tsv --ext freyja.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_AGGREGATE_ALL
//
// Cross-run "all_batches" aggregate — rescans EVERY *.freyja.tsv under
// results_dir (not just this run's samples), matching run_all.sh's
// all-batches step exactly. `tsvs` isn't read directly; it's only there to
// force this to run after the current batch's samples are published.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_AGGREGATE_ALL {
    publishDir "${params.results_dir}/_aggregate/all_batches", mode: 'copy'

    input:
    path tsvs
    val resultsDirAbs

    output:
    path "aggregated.tsv"

    script:
    """
    mkdir -p staging
    find "${resultsDirAbs}" -name "*.freyja.tsv" \
        -not -path "${resultsDirAbs}/_aggregate/*" \
        -exec cp {} staging/ \\;
    freyja aggregate staging/ --output aggregated.tsv --ext freyja.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_PLOT / FREYJA_PLOT_ALL
//
// Generates the lineage bar chart PDF from an aggregated TSV.
// Nextflow DSL2 disallows invoking the same process twice in one workflow,
// so the per-batch and all-batches variants are two small process blocks
// rather than one parameterized process.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_PLOT {
    publishDir "${params.results_dir}/_aggregate/${params.batch_name}", mode: 'copy'

    input:
    path aggregated

    output:
    path "lineage_plot.pdf"

    script:
    """
    freyja plot ${aggregated} \
        --output lineage_plot.pdf \
        --mincov 0
    """
}

process FREYJA_PLOT_ALL {
    publishDir "${params.results_dir}/_aggregate/all_batches", mode: 'copy'

    input:
    path aggregated

    output:
    path "lineage_plot.pdf"

    script:
    """
    freyja plot ${aggregated} \
        --output lineage_plot.pdf \
        --mincov 0
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: LONG_FORMAT / LONG_FORMAT_ALL
//
// One row per summarized-group entry and one row per individual-lineage
// entry (kept separate, not paired — see bin/long_format.py for why).
// Split into two process blocks for the same reason as FREYJA_PLOT above.
// ─────────────────────────────────────────────────────────────────────────────
process LONG_FORMAT {
    publishDir "${params.results_dir}/_aggregate/${params.batch_name}", mode: 'copy'

    input:
    path aggregated

    output:
    path "long_aggregated.csv"

    script:
    """
    long_format.py ${aggregated} long_aggregated.csv
    """
}

process LONG_FORMAT_ALL {
    publishDir "${params.results_dir}/_aggregate/all_batches", mode: 'copy'

    input:
    path aggregated

    output:
    path "long_aggregated.csv"

    script:
    """
    long_format.py ${aggregated} long_aggregated.csv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: TIME_SERIES / TIME_SERIES_ALL
//
// Metadata-driven time-series plot — ported from run_all.sh's
// aggregate_and_plot() helper. Silently produces nothing if the metadata
// has no valid dated rows for this aggregate (same as bash). Split into two
// process blocks for the same reason as FREYJA_PLOT above.
// ─────────────────────────────────────────────────────────────────────────────
process TIME_SERIES {
    publishDir "${params.results_dir}/_aggregate/${params.batch_name}", mode: 'copy'

    input:
    path aggregated
    path metadata

    output:
    path "lineage_timeseries.pdf", optional: true
    path "_times_metadata.csv", optional: true
    path "_aggregated_timed.tsv", optional: true

    script:
    """
    time_series_prep.py ${metadata} ${aggregated} _times_metadata.csv _aggregated_timed.tsv
    if [[ -f _times_metadata.csv && -f _aggregated_timed.tsv ]]; then
        freyja plot _aggregated_timed.tsv \
            --times _times_metadata.csv \
            --interval MS \
            --output lineage_timeseries.pdf \
            --mincov 0
    fi
    """
}

process TIME_SERIES_ALL {
    publishDir "${params.results_dir}/_aggregate/all_batches", mode: 'copy'

    input:
    path aggregated
    path metadata

    output:
    path "lineage_timeseries.pdf", optional: true
    path "_times_metadata.csv", optional: true
    path "_aggregated_timed.tsv", optional: true

    script:
    """
    time_series_prep.py ${metadata} ${aggregated} _times_metadata.csv _aggregated_timed.tsv
    if [[ -f _times_metadata.csv && -f _aggregated_timed.tsv ]]; then
        freyja plot _aggregated_timed.tsv \
            --times _times_metadata.csv \
            --interval MS \
            --output lineage_timeseries.pdf \
            --mincov 0
    fi
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW
//
// Connects all processes together. Data flows top to bottom.
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // Auto-detect the newest metadata CSV, same as run.sh Step 2 — only if
    // the caller didn't pass --metadata_file explicitly. Lives here (not at
    // top level) because this Nextflow version's parser rejects bare
    // statements/if-blocks outside a process/workflow/function body.
    if (!params.metadata_file) {
        def mdDir = new File(params.metadata_dir)
        if (mdDir.exists()) {
            def csvs = mdDir.listFiles({ d, n -> n.toLowerCase().endsWith('.csv') } as FilenameFilter)
            if (csvs) {
                params.metadata_file = csvs.sort { a, b -> b.lastModified() <=> a.lastModified() }[0].path
            }
        }
    }

    // Resolve the Freyja-bundled reference dynamically (same as
    // run_sample.sh's `python -c "import freyja..."`), instead of hardcoding
    // a path tied to one Python minor version. condaEnvPath follows the same
    // FREYJA_CONDA_ENV override as nextflow.config's `conda` directive, so
    // this resolves correctly on any machine, not just the Mac.
    def condaEnvPath = System.getenv('FREYJA_CONDA_ENV') ?: '/opt/anaconda3/envs/freyja-env'
    def freyjaDataDir = ["${condaEnvPath}/bin/python3", '-c',
        'import freyja, os; print(os.path.join(os.path.dirname(freyja.__file__), "data"))'
    ].execute().text.trim()
    params.ref = params.ref ?: "${freyjaDataDir}/NC_045512_Hu-1.fasta"
    params.gff = params.gff ?: "${freyjaDataDir}/NC_045512_Hu-1.gff"

    // Absolute results_dir/raw_dir — needed by FREYJA_AGGREGATE_ALL's `find`
    // and CLEANUP_INTERMEDIATES' `rm`, since process scripts execute in an
    // isolated work directory, not the launch directory, so a relative path
    // here would resolve wrong.
    def resultsDirAbs = new File(params.results_dir).absolutePath
    def rawDirAbs = new File(params.raw_dir).absolutePath

    // 0. Pre-flight primer-scheme check (advisory only)
    if (params.metadata_file) {
        PRIMER_CHECK(file(params.metadata_file), file(params.acc_file))
    }

    // 1. Read accession file — one channel item per accession.
    //    Strip \r explicitly (Windows-exported accession lists), matching
    //    run.sh's `sed 's/\r//'` — otherwise a stray \r ends up baked into
    //    the accession value and every downstream file lookup misses.
    accessions = Channel
        .fromPath(params.acc_file)
        .splitText()
        .map { it.replaceAll('\r', '').trim() }
        .filter { it && !it.startsWith('#') }

    // 2. For accessions already downloaded, skip DOWNLOAD and pick up files directly.
    new_reads = DOWNLOAD(accessions)
    existing_reads = accessions
        .filter { file("${params.raw_dir}/${it}_R1.fastq.gz").exists() &&
                  file("${params.raw_dir}/${it}_R2.fastq.gz").exists() }
        .map    { acc -> tuple(acc,
                              file("${params.raw_dir}/${acc}_R1.fastq.gz"),
                              file("${params.raw_dir}/${acc}_R2.fastq.gz")) }
    all_reads = new_reads.mix(existing_reads)

    // Fresh Freyja barcodes — once per invocation, shared across every
    // sample in this batch (see FREYJA_UPDATE comment for why not per-sample).
    // .first() converts this into a genuine value channel — without it,
    // Nextflow would pair its single emission positionally with the
    // multi-sample variants/qc_inputs channels and only run for sample #1.
    freyja_db = FREYJA_UPDATE().first()

    // 3. Per-sample steps — each runs on all samples in parallel automatically
    FASTP(all_reads)
    trimmed     = FASTP.out.reads
    aligned     = ALIGN(trimmed)
    trimmed_bam = IVAR_TRIM(aligned, file(params.bed))
    variants    = FREYJA_VARIANTS(trimmed_bam)
    demixed     = FREYJA_DEMIX(variants, freyja_db)

    // Annotate each freyja.tsv with site + collection_date via the
    // cdc_site_code.csv join — only if both metadata and cdc_site_code.csv
    // are actually available; otherwise pass the unannotated result through
    // rather than failing the whole run over an optional enrichment step.
    if (params.metadata_file && file(params.cdc_site_code).exists()) {
        ANNOTATE_SITE_DATE(demixed, file(params.cdc_site_code), file(params.metadata_file))
        demixed_final = ANNOTATE_SITE_DATE.out
    } else {
        demixed_final = demixed
    }

    // Delete raw FASTQ + BAMs now that this sample's freyja.tsv exists —
    // see CLEANUP_INTERMEDIATES comment for why this is safe and why BAM
    // retention isn't worth the storage cost at production volume.
    CLEANUP_INTERMEDIATES(demixed_final, rawDirAbs, resultsDirAbs)

    // 4. Per-batch aggregate + plot
    all_tsvs         = demixed_final.map { acc, tsv -> tsv }.collect()
    batch_aggregated = FREYJA_AGGREGATE(all_tsvs)
    FREYJA_PLOT(batch_aggregated)
    LONG_FORMAT(batch_aggregated)

    // 5. QC run log (per batch) — needs both freyja.tsv and depths.tsv per sample
    qc_inputs   = demixed_final.join(variants).map { acc, freyja_tsv, var_tsv, depths_tsv -> tuple(acc, freyja_tsv, depths_tsv) }
    acc_names   = qc_inputs.map { acc, f, d -> acc }.collect()
    freyja_tsvs = qc_inputs.map { acc, f, d -> f }.collect()
    depths_tsvs = qc_inputs.map { acc, f, d -> d }.collect()
    QC_RUNLOG(acc_names, freyja_tsvs, depths_tsvs, freyja_db)

    // 6. Time-series (per batch)
    if (params.metadata_file) {
        TIME_SERIES(batch_aggregated, file(params.metadata_file))
    }

    // 7. All-batches aggregate + plot (rescans every *.freyja.tsv on disk)
    all_batches_aggregated = FREYJA_AGGREGATE_ALL(all_tsvs, resultsDirAbs)
    FREYJA_PLOT_ALL(all_batches_aggregated)
    LONG_FORMAT_ALL(all_batches_aggregated)

    // 8. Time-series (all-batches)
    if (params.metadata_file) {
        TIME_SERIES_ALL(all_batches_aggregated, file(params.metadata_file))
    }
}
