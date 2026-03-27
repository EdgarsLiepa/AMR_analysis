#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.metaphlan_db = "DB/metaphlan_db"
params.resfinder_db = "DB/resfinder_db/all"

params.projPath = ""
params.outdir = "${params.projPath}/results/202x-xx-xx"

params.samplesheet = "${params.projPath}/data/202x-xx-xx/samplesheet_202x-xx-xx.tsv"

params.reads_dir = "${params.projPath}/data/202x-xx-xx/merged_reads"  // Directory containing all FASTQ files

// Define the folder containing the index files
hostile_idx_ch = Channel.fromPath('.local/share/hostile', type: 'dir').first()


//params.human_reference = '/path/to/human_genome.fasta'
    
// Fastp quality parameters
params.qualified_quality_phred = 20
params.unqualified_percent_limit = 40
params.min_length = 50
params.cut_window_size = 4
params.cut_mean_quality = 20



include { FASTQC_RAW; FASTP_DECONTAMINATE; FASTQC_CLEAN } from "../modules/preProc.nf" 

include { MULTIQC as MULTIQC_RAW } from "../modules/preProc.nf" addParams(OUTPUT: "${params.outdir}/multiqc_raw")
include { MULTIQC as MULTIQC_CLEAN } from "../modules/preProc.nf" addParams(OUTPUT: "${params.outdir}/multiqc_clean")
include { REMOVE_HOST_HOSTILE } from "../modules/preProc.nf" addParams(OUTPUT: "${params.outdir}/host_removed")


process ARG_KMA {
    tag "${meta.id}" 
    publishDir "${params.outdir}/KMA", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.res"), emit: res
    path "${meta.id}_ARG_counts.tsv", emit: counts    

    cpus 8
    memory '64 GB'
    time = '1d'

    script:
    
    """

    # -ipe: Paired-end reads
    # -t_db: KMA database prefix
    # -1t1: One read to one template (forces strict assignment, prevents double-counting)
    # -nc: No consensus building (saves memory/time)
    ~/tools/kma/kma -ipe ${reads[0]} ${reads[1]} -t_db ${params.resfinder_db} -o ${meta.id} -1t1 -nc -t ${task.cpus}
    
    # 3. Calculate Normalized Metrics (Depth)
    # Extract the Template Name (col 1) and Depth (col 9) from the KMA .res file. 
    # Depth is length-normalized, making it vastly superior to raw counts.
    echo "${meta.id}" > ${meta.id}_ARG_counts.tsv

    grep -v "^#" ${meta.id}.res \
        | awk 'BEGIN{OFS="\\t"} {print \$1, \$4, \$6, \$9}' \
        | sed '1i Gene\\tTemplate_Length\\tTemplate_Coverage\\tDepth' \
        > ${meta.id}_ARG_counts.tsv

    """

}

process ARG_MAPPED {
    tag "${meta.id}"        
    publishDir "${params.outdir}/Bowtie", mode: 'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_ARG_mapped.bam"), path("${meta.id}_ARG_mapped.bam.bai")

    cpus 6
    memory '64 GB'
    time = '1d'

    script:
    
    """
    module load 'bio/samtools-1.9'
    
    # --------------------------------------------------------------------------------
    # Filter, compress, and sort the alignments
    # --------------------------------------------------------------------------------
    # Step 1: samtools view -h   -> Output SAM format including the header
    # Step 2: awk                -> Filter alignments based on pair status:
    #                                  \$7!="=" : Keep discordant pairs (mate mapped to different contig)
    #                                  (\$7=="=" && and(\$2,0x40)) : For concordant pairs, ONLY keep Read 1 to avoid double-counting
    # Step 3: samtools view -Su  -> Convert back to uncompressed BAM format for speed
    # Step 4: samtools sort      -> Sort by coordinates and output final BAM file
    # --------------------------------------------------------------------------------

    # Filter, compress, and sort the alignments
    samtools view -h -@ ${task.cpus} ${bam} | \\
    awk '\$7!="=" || (\$7=="=" && and(\$2,0x40)) {print \$0}' | \\
    samtools view -Su -@ ${task.cpus} - | \\
    samtools sort -@ ${task.cpus} -o ${meta.id}_ARG_mapped.bam
    
    samtools index -@ ${task.cpus} ${meta.id}_ARG_mapped.bam

    """
}


process MERGE_ARG_COUNTS {
    tag "Merge_ARG_KMA"
    
    // 1. Environment Management: Guarantee pandas is available
    container 'quay.io/biocontainers/pandas:1.4.3'
    
    // 2. Resource Management: Prevent OOM errors on large outer joins
    cpus 1
    memory '8 GB'
    
    publishDir "${params.outdir}/Final_Results", mode: 'copy'

    input:
    path count_files

    output:
    path "Master_ARG_RawDepth_Matrix.tsv", emit: depth_matrix
    path "Master_ARG_TPM_Matrix.tsv",      emit: tpm_matrix
    
    script:
    // 3. Groovy Array Injection: Convert Nextflow path list to a valid Python list string
    def py_file_list = count_files.collect { "'${it}'" }.join(', ')

    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    input_files = [${py_file_list}]

    if not input_files:
        raise RuntimeError("No input files provided by Nextflow.")

    depth_dfs = []
    tpm_dfs = []

    for f in input_files:
        # Read the file and set 'Gene' as the index
        df = pd.read_csv(f, sep='\\t', index_col='Gene')
        
        # Extract sample name (e.g., 'LM0563_ARG_counts.tsv' -> 'LM0563')
        sample_name = os.path.basename(f).replace('_ARG_counts.tsv', '')
        
        # 1. Isolate Raw Depth
        depth_series = df[['Depth']].copy()
        depth_series.columns = [sample_name]
        depth_dfs.append(depth_series)
        
        # 2. Calculate TPM
        tpm_series = depth_series.copy()
        total_depth = tpm_series[sample_name].sum()
        
        if total_depth > 0:
            # TPM = (Gene Depth / Total Sample Depth) * 1,000,000
            tpm_series[sample_name] = (tpm_series[sample_name] / total_depth) * 1e6
            
        tpm_dfs.append(tpm_series)

    # Merge all samples into master dataframes
    master_depth = pd.concat(depth_dfs, axis=1, join='outer').fillna(0).round(2)
    master_tpm = pd.concat(tpm_dfs, axis=1, join='outer').fillna(0).round(2)

    # Ensure the index has a clean name
    master_depth.index.name = 'Gene'
    master_tpm.index.name = 'Gene'

    # Save both files
    master_depth.to_csv("Master_ARG_RawDepth_Matrix.tsv", sep='\\t')
    master_tpm.to_csv("Master_ARG_TPM_Matrix.tsv", sep='\\t')

    """
}

workflow {

    // Read samplesheet and construct file paths
    def samples = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .unique { row -> row.sample_id }
        .map { row ->
            def meta = [id: row.sample_id]
            def read1 = file("${params.reads_dir}/${row.sample_id}_R1.fastq.gz")
            def read2 = file("${params.reads_dir}/${row.sample_id}_R2.fastq.gz")
            [meta, [read1, read2]]
        }
    
    // QC on raw reads
    FASTQC_RAW(samples)
    
    // Aggregate raw QC
    MULTIQC_RAW(
        FASTQC_RAW.out.zip.map { meta, zip -> zip }.collect() 
    )
    
    // Decontaminate and quality filter
    FASTP_DECONTAMINATE(samples)
    
    // QC on clean reads
    FASTQC_CLEAN(FASTP_DECONTAMINATE.out.reads)
    
    // --- Aggregate clean QC with fastp reports ---
    // use .mix() to combine both channels into one, then .collect() them into a single list
    def clean_reports = FASTQC_CLEAN.out.zip.map { meta, zip -> zip }
    .mix( FASTP_DECONTAMINATE.out.json.map { meta, json -> json } )
    .collect()

    MULTIQC_CLEAN(clean_reports)

    REMOVE_HOST_HOSTILE(FASTP_DECONTAMINATE.out.reads, hostile_idx_ch)


    // Map to ResFinder, Filter, and Quantify using KMA
    ARG_KMA(REMOVE_HOST_HOSTILE.out.reads)

    // Merge everything using your robust Pandas script
    MERGE_ARG_COUNTS(ARG_KMA.out.counts.collect())

}
