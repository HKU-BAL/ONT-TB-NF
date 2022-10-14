#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.sample_name="TB"
params.read_fq="$launchDir/read.fq"
params.threads="32"
params.output_dir="$launchDir/out"

params.ref="$baseDir/data/Mycobacterium_tuberculosis_H37Rv_genome_v3.fasta"
params.ref_index="${params.ref}.fai"
params.gene_bed="$baseDir/data/AMR_gene_H37RV_v3.bed"
params.nanofilt_options=""

params.TB_script_dir="$baseDir/scripts"
params.C3_model_n="r941_prom_sup_g5014"
params.alignment_work_space="20G"

params.help = false

if( params.help ) {

log.info """

A Nextflow pipeline for Mycobacterium tuberculosis (TB) antibiotic resistance gene analysis with ONT data.

=============================================
Usage:
    nextflow run run_tb.nf --sample_name [SAMPLE_NAME] --read_fq [READ FQ] --threads [THREADS] --output_dir [OUTPUT DIR] 
Input:
    * --sample_name: sample name. Default [${params.sample_name}]
    * --read_fq: Path of read FQ file. Default [read.fq]
    * --threads: number of threads for running. Default [${params.threads}]
    * --output_dir: name of output directory. Default [out]
    * --nanofilt_options: read filtering option. Default [None]

more information are available at [Gtihub page](https://github.com/HKU-BAL/ONT-TB-NF)
"""
    exit 0
}

log.info """\

Analysis TB data
================================
sample name        : $params.sample_name
input FQ           : $params.read_fq
threads            : $params.threads
output             : $params.output_dir
reference          : $params.ref
gene bed           : $params.gene_bed
Clair3 model       : $params.C3_model_n
QC NanoFilt option : [$params.nanofilt_options]

"""

process run_QC {
    debug true
	publishDir "$params.output_dir/1_qc", mode: 'copy'

    input:
    path read_fq

    output:
    path "out.fq"
    path "*.html"

    """
    # run fastQC, check read base quality, read quality, GC content, etc.
    fastqc -f fastq -o . -t ${params.threads} ${read_fq}

    # using nanofilt to filter
    # e.g. --nanfil_options "-q 7 " to filter read with mean quailty < 7
    if [[ "${params.nanofilt_options}" != "" ]]; then
        echo " run nanofilt_options"
        NanoFilt ${params.nanofilt_options} ${read_fq} > out.fq
	else
        ln -s ${read_fq} out.fq
    fi

    """
}

process run_aln {
	publishDir "$params.output_dir/2_aln", mode: 'copy'

    input:
    path ref 
    path read_fq

    output:
    path "${params.sample_name}_ori.bam"
    path "${params.sample_name}_ori.bam.bai"

    """
	minimap2 -ax map-ont -t $params.threads -I params.alignment_work_space $ref $read_fq | samtools sort  -@ $params.threads | samtools view -F 2048 -b > ${params.sample_name}_ori.bam
	samtools index ${params.sample_name}_ori.bam
    """
}


process run_aln_filtering {
	publishDir "$params.output_dir/2_aln", mode: 'copy'

    input:
    path "${params.sample_name}_ori.bam"
    path "${params.sample_name}_ori.bam.bai"

    output:
    path "${params.sample_name}.bam"
    path "${params.sample_name}.bam.bai"

    """
	python ${params.TB_script_dir}/filter_bam_by_AS_ratio.py ${params.sample_name}_ori.bam 1.2 0 ${params.sample_name}.bam
    """
}

process run_check_gene_coverage {
	publishDir "$params.output_dir/2_aln", mode: 'copy'
	
    input:
    path "${params.sample_name}.bam"
    path "${params.sample_name}.bam.bai"

	output:
	path "*.summary.txt"
	path "${params.sample_name}.regions.bed.gz"
	path "mapping.stats"
	stdout

	"""
	mosdepth -F 3844 -t $params.threads -b ${params.gene_bed} "${params.sample_name}" "${params.sample_name}.bam"
	samtools flagstat "${params.sample_name}.bam" > mapping.stats
    echo "number of TB read found:"
	head -n 1 mapping.stats
    zcat ${params.sample_name}.regions.bed.gz | awk '{ sum += \$5; n++ } END { if (n > 0) print "average AMR gene coverage: " sum / n; }'
	#echo 'contig\t start\t end\t name\t coverage'
	#zcat ${params.sample_name}.regions.bed.gz
	"""
}



process run_variant_calling {
    containerOptions "--cpus=${params.threads}"
	publishDir "$params.output_dir/3_vc", mode: 'copy'

	input:
    path "ref.fa"
    path "ref.fa.fai"
    path "${params.sample_name}.bam"
    path "${params.sample_name}.bam.bai"

	output:
    path "clair3_out/merge_output.vcf.gz"
    path "clair3_out/merge_output.vcf.gz.tbi"
    path "clair3_out"
	

	"""
	/opt/bin/run_clair3.sh \
	--bam=${params.sample_name}.bam \
	--ref_fn=ref.fa \
	--threads=${params.threads} \
	--chunk_size=100000 \
	--platform="ont" \
	--sample_name=${params.sample_name} \
	--model_path="/opt/models/${params.C3_model_n}" \
	--output=clair3_out \
	--include_all_ctgs \
	--no_phasing_for_fa \
	--haploid_precise > log


	"""
}

process run_tb_profiler {
    containerOptions "--cpus=${params.threads}"
	publishDir "$params.output_dir/4_tb", mode: 'copy'

	input:
    path "${params.sample_name}.bam"
    path "${params.sample_name}.bam.bai"

	output:
    path "results"
	
	"""
	tb-profiler profile \
	--bam ${params.sample_name}.bam \
	--platform nanopore \
	--csv \
	--pdf \
	--threads ${params.threads} \
	--prefix ${params.sample_name} > log

	"""
    
}

process run_get_consensus {
	debug true
	publishDir "$params.output_dir/4_cns", mode: 'copy'
	
    input:
    path "${params.sample_name}.bam"
    path "${params.sample_name}.bam.bai"

	output:
	path "${params.sample_name}.fa" 
	path "flye"
	path "target.fa"

	"""
	samtools fasta "${params.sample_name}.bam" >  "${params.sample_name}.fa" 
	flye --nano-hq "${params.sample_name}.fa" --out-dir flye --threads ${params.threads} --meta --genome-size 4.4m || true

	if [ ! -f "flye/assembly.fasta" ]; then
		echo "no consensus generated, will use original FA"
		ln -s "${params.sample_name}.fa" target.fa 
	else
		echo "consensus found, will use the consensus"
		ln -s flye/assembly.fasta target.fa 
	fi
	echo "target read/consensus"
	grep ">" target.fa | wc -l
	"""
}



workflow {
    (read_fq, _) = run_QC(params.read_fq)
    run_aln(params.ref, read_fq)
    run_aln_filtering(run_aln.out)
    (_, _, _, o) = run_check_gene_coverage(run_aln_filtering.out)
	o.view{"$it"}
    run_variant_calling(params.ref, params.ref_index, run_aln_filtering.out)
	//(_, _, cns_f) = run_get_consensus(run_aln_filtering.out)
	run_tb_profiler(run_aln_filtering.out)
}

workflow.onComplete {
    print "================================"
    print "Finish TB analysis pipeline"
    print "[1] QC results at:                                        ${params.output_dir}/1_qc"
    print "[2] Aligment results at:                                  ${params.output_dir}/2_aln"
    print "    | TB aligned bam at:                                  ${params.output_dir}/2_aln/${params.sample_name}.bam"
    print "    | Alignment results statistics at:                    ${params.output_dir}/2_aln/mapping.stats"
    print "    | AMR genes coverages at:                             ${params.output_dir}/2_aln/${params.sample_name}.regions.bed.gz"
    print "[3] Variant calling results at:                           ${params.output_dir}/3_vc"
    print "    | variant calling result at:                          ${params.output_dir}/3_vc/clair3_out/merge_output.vcf.gz"
    print "[4] TB analysis report at:                                ${params.output_dir}/4_tb"
    print "    | TB pdf report at:                                   ${params.output_dir}/4_tb/results/${params.sample_name}.results.pdf"

    print """
================================
Command line: ${workflow.commandLine}
Completed at: ${workflow.complete}
Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
exit status : ${workflow.exitStatus}
"""
}
