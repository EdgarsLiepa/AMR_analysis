process FASTQC_RAW {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastqc_raw"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

process FASTP_DECONTAMINATE {
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/fastp:1.1.0--3a0b489903631b1c'
    publishDir "${params.outdir}/trimmed_reads", mode: 'symlink', pattern: "*.fastq.gz"
    publishDir "${params.outdir}/fastp_reports", mode: 'symlink', pattern: "*.{html,json,log}"
    
    cpus 4

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_clean_R*.fastq.gz"), emit: reads
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.log"), emit: log
    
    script:
    def prefix = meta.id
    def read1 = reads[0]
    def read2 = reads[1]
    """

    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${prefix}_clean_R1.fastq.gz \\
        --out2 ${prefix}_clean_R2.fastq.gz \\
        --thread ${task.cpus} \\
        --html ${prefix}_fastp.html \\
        --json ${prefix}_fastp.json \\
        --detect_adapter_for_pe \\
        --cut_front --cut_tail \\
        2>| >(tee ${prefix}.fastp.log >&2)

    """
}

process FASTQC_CLEAN {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastqc_clean"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}


process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    publishDir "${params.OUTPUT}", mode: 'copy'
    
    input:
    path(fastqc_files, stageAs: "?/*")
    
    output:
    path("multiqc_report.html"), emit: html
    path("multiqc_report_data"), emit: data
    
    script:
    """
    multiqc . -n multiqc_report.html
    """
}



process REMOVE_HOST_HOSTILE {
    tag "${meta.id}"
    publishDir "${params.outdir}/final_microbial_reads", mode: 'copy'

    container 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'

    cpus 16

    input:
    tuple val(meta), path(reads)
    path index_dir

    output:
    // Hostile appends 'clean' to the output filenames automatically
    tuple val(meta), path("*.clean_[12].fastq.gz"), emit: reads
    path "${meta.id}.hostile.log", emit: log    


    script:
    
    def r1 = reads[0]
    def r2 = reads[1]

    """

    # add the cache redirect!
    export HOSTILE_CACHE_DIR=\$PWD/hostile_cache

    hostile clean \\
        --fastq1 ${r1} \\
        --fastq2 ${r2} \\
        --aligner bowtie2 \\
        --index ${index_dir}/human-t2t-hla-argos985 \\
        --threads ${task.cpus} \\
        --output . \\
        2> ${meta.id}.hostile.log
    """
}

