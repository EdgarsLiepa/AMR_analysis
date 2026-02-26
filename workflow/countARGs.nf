#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.metaphlan_db = 
params.resfinder_db = "DB/resfinder_db/all"

params.projPath = 
params.outdir = 

params.samplesheet = 

params.reads_dir =   // Directory containing all FASTQ files

// Define the folder containing the index files
hostile_idx_ch = Channel.fromPath('', type: 'dir').first()


    
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



process ARG_BOWTIE {
    tag "${meta.id}" 
    publishDir "${params.outdir}/Bowtie", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_ARG.bam")

    cpus 6
    memory '64 GB'
    time = '1d'

    script:
    """
    module load 'bio/bowtie2/2.4.2'
    module load 'bio/samtools-1.9'

    # --------------------------------------------------------------------------------
    # Align microbial reads to the ResFinder database
    # --------------------------------------------------------------------------------
    # -x : Path to the ResFinder Bowtie2 index prefix
    # -1 / -2 : Forward and reverse input reads
    # -D 20 -R 3 -N 1 -L 20 -i S,1,0.5 : Custom search and seed parameters.
    #      (These heavily increase alignment sensitivity to catch highly mutated 
    #       or novel resistance genes that slightly differ from the reference)
    # --------------------------------------------------------------------------------


    bowtie2 -x ${params.resfinder_db} -1 ${reads[0]} -2 ${reads[1]} \
    -D 20 -R 3 -N 1 -L 20 -i S,1,0.5 \
    --threads ${task.cpus} | samtools view -@ ${task.cpus} -Sb - > ${meta.id}_ARG.bam
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

process ARG_COUNT {
    tag "${meta.id}"
    publishDir "${params.outdir}/Bowtie/Counts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.id}_ARG_counts"

    cpus 2
    memory '6 GB'

    script:
    """
    
    module load 'bio/samtools-1.9'
    
    echo -e ${meta.id} > ${meta.id}_ARG_counts
    samtools idxstats ${bam} | grep -v "*" | cut -f1,3 >> ${meta.id}_ARG_counts
    """
}


process MERGE_ARG_COUNTS {
    publishDir "${params.outdir}/Final_Results", mode: 'copy'

    input:
    // Takes a collected list of all individual count files
    path count_files

    output:
    path "Master_ARG_Count_Matrix.tsv"

    script:
    """
    # 1. Extract the 'Gene' column from the very first file to use as our row names
    first_file=\$(ls *_ARG_counts.tsv | head -n 1)
    cut -f1 \$first_file > __genes.tmp

    # 2. Extract ONLY the count columns (Column 2) from all files and save as temps
    for f in *_ARG_counts.tsv; do
        cut -f2 \$f > \${f}.col.tmp
    done

    # 3. Paste the genes and all the sample count columns together side-by-side
    paste __genes.tmp *.col.tmp > Master_ARG_Count_Matrix.tsv

    # 4. Clean up the temporary files
    rm *.tmp
    """
}

workflow {

    // Read samplesheet and construct file paths
    def samples = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
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
    // We use .mix() to combine both channels into one, then .collect() them into a single list
    def clean_reports = FASTQC_CLEAN.out.zip.map { meta, zip -> zip }
    .mix( FASTP_DECONTAMINATE.out.json.map { meta, json -> json } )
    .collect()

    MULTIQC_CLEAN(clean_reports)

    REMOVE_HOST_HOSTILE(FASTP_DECONTAMINATE.out.reads, hostile_idx_ch)


    // METAPHLAN(combined_read_pairs_ch)
    ARG_BOWTIE(REMOVE_HOST_HOSTILE.out.reads)
    ARG_MAPPED(ARG_BOWTIE.out)
    count_ch = ARG_COUNT(ARG_MAPPED.out)

    // Collects all individual files into a single list and merges them
    // MERGE_ARG_COUNTS( count_ch.collect() )

}
