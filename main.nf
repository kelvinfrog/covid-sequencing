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
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMETERS — override any of these from the command line
// e.g. --threads 8  or  --max_reads 3000000
// ─────────────────────────────────────────────────────────────────────────────
params.acc_file    = ""
params.raw_dir     = "data/raw"
params.results_dir = "data/results"
params.bed         = "${projectDir}/data/bed/ARTIC_V5.3.2.bed"
params.ref         = "/opt/anaconda3/envs/freyja-env/lib/python3.13/site-packages/freyja/data/NC_045512_Hu-1.fasta"
params.max_reads   = 2000000
params.threads     = 4
params.min_quality = 20
params.batch_name  = "nextflow_batch"

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: DOWNLOAD
//
// Downloads the first 2M reads from SRA using fastq-dump.
// The `when` block skips the download if the files already exist in data/raw/.
// This replaces the [EXISTS] check in run.sh.
// ─────────────────────────────────────────────────────────────────────────────
process DOWNLOAD {
    publishDir params.raw_dir, mode: 'copy'

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_R1.fastq.gz"), path("${accession}_R2.fastq.gz")

    when:
    !file("${params.raw_dir}/${accession}_R1.fastq.gz").exists()

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
// Note on publishDir: the path uses a closure { } so that ${accession}
// is resolved at runtime (when we know the actual sample name), not at
// script definition time. This is required in Nextflow DSL2 for dynamic paths.
// ─────────────────────────────────────────────────────────────────────────────
process FASTP {
    publishDir { "${params.results_dir}/${accession}/trimmed" }, mode: 'copy'

    input:
    tuple val(accession), path(r1), path(r2)

    output:
    tuple val(accession), path("${accession}_trimmed_R1.fastq.gz"), path("${accession}_trimmed_R2.fastq.gz")

    script:
    """
    fastp \
        --in1 ${r1} --in2 ${r2} \
        --out1 ${accession}_trimmed_R1.fastq.gz \
        --out2 ${accession}_trimmed_R2.fastq.gz \
        --json ${accession}_fastp.json \
        --html ${accession}_fastp.html \
        --thread ${params.threads} \
        --qualified_quality_phred ${params.min_quality} \
        --length_required 50 \
        --detect_adapter_for_pe
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: ALIGN
//
// Aligns trimmed reads to SARS-CoV-2 reference, then sorts and indexes BAM.
// Combined into one process because they always run together.
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
//
// path(bed) — Nextflow stages the BED file into the working directory so
// the script can reference it by filename only, regardless of where it lives.
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
    ivar trim -i ${bam} -b ${bed} -p ${accession}.trimmed -e
    samtools sort -@ ${params.threads} ${accession}.trimmed.bam \
        -o ${accession}.trimmed.sorted.bam
    samtools index ${accession}.trimmed.sorted.bam
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_VARIANTS
//
// Calls variants and measures depth across the genome.
// Most CPU-intensive step — samtools mpileup scans every genome position.
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
        --minq ${params.min_quality}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_DEMIX
//
// Estimates lineage abundances — the core freyja step.
// Produces the .freyja.tsv result file for each sample.
// ─────────────────────────────────────────────────────────────────────────────
process FREYJA_DEMIX {
    publishDir { "${params.results_dir}/${accession}" }, mode: 'copy'

    input:
    tuple val(accession), path(variants), path(depths)

    output:
    tuple val(accession), path("${accession}.freyja.tsv")

    script:
    """
    freyja demix ${variants} ${depths} \
        --output ${accession}.freyja.tsv \
        --confirmedonly \
        --eps 0.01 \
        --covcut 10
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESS: FREYJA_AGGREGATE
//
// Combines all samples' freyja.tsv files into one aggregated TSV for plotting.
//
// Key concept: this process receives ALL samples at once, not one at a time.
// That's because in the workflow we use .collect() which waits for every
// sample to finish before sending them here together.
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
// PROCESS: FREYJA_PLOT
//
// Generates the lineage bar chart PDF from the aggregated TSV.
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

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW
//
// Connects all processes together. Data flows top to bottom.
// Each line feeds the output of one process into the input of the next.
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // 1. Read accession file — creates one channel item per accession
    accessions = Channel
        .fromPath(params.acc_file)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }

    // 2. For accessions already downloaded, skip DOWNLOAD and pick up files directly.
    //    new_reads    — freshly downloaded (DOWNLOAD's when: block handles the check)
    //    existing_reads — already in data/raw/, grabbed directly from disk
    //    all_reads    — merges both so the rest of the pipeline sees all samples
    new_reads = DOWNLOAD(accessions)
    existing_reads = accessions
        .filter  { file("${params.raw_dir}/${it}_R1.fastq.gz").exists() }
        .map     { acc -> tuple(acc,
                               file("${params.raw_dir}/${acc}_R1.fastq.gz"),
                               file("${params.raw_dir}/${acc}_R2.fastq.gz")) }
    all_reads = new_reads.mix(existing_reads)

    // 3. Per-sample steps — each runs on all samples in parallel automatically
    trimmed     = FASTP(all_reads)
    aligned     = ALIGN(trimmed)
    trimmed_bam = IVAR_TRIM(aligned, file(params.bed))
    variants    = FREYJA_VARIANTS(trimmed_bam)
    demixed     = FREYJA_DEMIX(variants)

    // 4. Aggregate — .collect() waits for ALL samples before aggregating
    //    .map drops the accession name — we only need the file paths here
    all_tsvs   = demixed.map { acc, tsv -> tsv }.collect()
    aggregated = FREYJA_AGGREGATE(all_tsvs)

    // 5. Plot
    FREYJA_PLOT(aggregated)
}
