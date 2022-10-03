#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.ref="$baseDir/data/ref.fa"
params.read_fq="$baseDir/data/read.fq"
params.sample_name="TB_1"
params.amplicon_bed="$baseDir/data/amplicon.bed"
params.TB_script_dir="/autofs/bal31/jhsu/home/data/TB/scripts"
params.C3_model_n="r941_prom_sup_g5014"

params.threads="32"

params.output_dir="$baseDir/out"

log.info """\

Analysis TB data
================================
sample name  : $params.sample_name
input FQ     : $params.read_fq
reference    : $params.ref
amplicon bed : $params.amplicon_bed
output       : $params.output_dir
threads      : $params.threads
Clair3 model : $params.C3_model_n

"""

process run_QC {
	publishDir "$params.output_dir/1_qc", mode: 'copy'

    input:
    path read_fq

    output:
    path "$params.sample_name"

    """
    NanoPlot -o $params.sample_name -p $params.sample_name --fastq $read_fq -t $params.threads --raw
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
	minimap2 -ax map-ont -t $params.threads -I 1000G $ref $read_fq | samtools sort  -@ $params.threads | samtools view -F 2048 -b > ${params.sample_name}_ori.bam
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
	mosdepth -F 3844 -t $params.threads -b ${params.amplicon_bed} "${params.sample_name}" "${params.sample_name}.bam"
	samtools flagstat "${params.sample_name}.bam" > mapping.stats
	head -n 2 mapping.stats
	#echo 'contig\t start\t end\t name\t coverage'
	#zcat ${params.sample_name}.regions.bed.gz
	"""
}

process run_variant_calling {
	debug true
	publishDir "$params.output_dir/3_vc"

	input:
    path "${params.sample_name}.bam"
    path "${params.sample_name}.bam.bai"

    output:
    val 0
	

	"""
docker run \
--mount type=bind,source=${params.ref},target=/ref.fa \
--mount type=bind,source=${params.ref}.fai,target=/ref.fa.fai \
--mount type=bind,source=${params.amplicon_bed},target=/tar.bed \
--mount type=bind,source=$params.output_dir/2_aln/${params.sample_name}.bam,target=/input.bam \
--mount type=bind,source=$params.output_dir/2_aln/${params.sample_name}.bam.bai,target=/input.bam.bai \
-v ${params.output_dir}/3_vc:/out \
hkubal/clair3:v0.1-r12 /opt/bin/run_clair3.sh \
--bam=/input.bam \
--ref_fn=/ref.fa \
--bed_fn=/tar.bed \
--threads=${params.threads} \
--chunk_size=100000 \
--platform="ont" \
--sample_name="${params.sample_name}" \
--model_path="/opt/models/${params.C3_model_n}" \
--output=/out \
--include_all_ctgs \
--no_phasing_for_fa \
--haploid_precise 


	"""
}

process run_get_consensus {
	debug true
	publishDir "$params.output_dir/4_cns", mode: 'copy'
	
    input:
    val x

	output:
	path "tar.vcf.gz" 
	path "tar.vcf.gz.tbi" 
	path "amplicon.region"
	path "target.fa"

	"""

    cp ${params.output_dir}/3_vc/merge_output.vcf.gz ori.vcf.gz
    tabix ori.vcf.gz
    bcftools view -f 'PASS,.' ori.vcf.gz > tmp.vcf
    bedtools intersect -a tmp.vcf -b ${params.amplicon_bed} -header > tmp1.vcf
    bcftools view tmp1.vcf -O z -o tar.vcf.gz
    tabix tar.vcf.gz
    awk '{split(\$0,a," "); printf "%s:%s-%s\\n", a[1],a[2],a[3] }' ${params.amplicon_bed} > amplicon.region
    samtools faidx ${params.ref} -r amplicon.region | bcftools consensus tar.vcf.gz >  target.fa

	echo "target read/consensus"
	grep ">" target.fa | wc -l
	"""
}


process run_rgi {
	debug true
	publishDir "$params.output_dir/5_rgi"
	
    input:
    path "target.fa"


	"""
docker run \
--mount type=bind,source=${params.output_dir}/4_cns/target.fa,target=/input.fa \
-v ${params.output_dir}/5_rgi:/out \
finlaymaguire/rgi:latest \
rgi main -i /input.fa -o /out/rgi -t contig -n ${params.threads} --clean --low_quality
    wc -l ${params.output_dir}/5_rgi/rgi.txt

	"""
}

workflow {
    run_QC(params.read_fq)
    run_aln(params.ref, params.read_fq)
    run_aln_filtering(run_aln.out)
    (_, _, _, o) = run_check_gene_coverage(run_aln_filtering.out)
	o.view{"$it"}
    run_variant_calling(run_aln_filtering.out)
	(_, _, _, cns_f) = run_get_consensus(run_variant_calling.out)
	run_rgi(cns_f)
}

workflow.onComplete {
    print "================================"
    print "Finish TB analysis pipeline"
    print "QC results at:                ${params.output_dir}/1_qc"
    print "Aligment results at:          ${params.output_dir}/2_aln"
    print "Variant calling results at:   ${params.output_dir}/3_vc"
    print "Consensus calling results at: ${params.output_dir}/4_cns"
    print "TB analysis report at:        ${params.output_dir}/5_rgi"
    print "TB's json report can be uploaded to [https://card.mcmaster.ca/analyze/externalrgi] for visulization"

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
