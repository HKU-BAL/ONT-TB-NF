# ONT-TB-NF

Comprehensive pipeline for detection of *Mycobacterium tuberculosis* (TB) antibiotic resistance gene from ONT adaptive sequencing and amplicon data. 

Pipeline including quality control, alignment, variant calling, consensus generation and antimicrobial resistance prediction. 

## Features

- One command pipeline from sequencing data to TB analysis report. 
- Tailor-made for Adaptive sequencing data and Amplicon data
- Basecalling with ONT's Guppy, the whole pipeline can start from fast5 files.


## Quickstart 

Install [Nextflow](https://www.nextflow.io/) by using the following command: 

    curl -s https://get.nextflow.io | bash 
    
Install required packages and activate envirment with [conda](https://conda.io/projects/conda/en/latest/index.html):

    docker pull hkubal/clair3:v0.1-r12
    docker pull quay.io/biocontainers/tb-profiler:4.3.0--pypyh5e36f6f_0
    conda create -n ont_tb samtools=1.15.1 minimap2=2.24 nanoplot=1.40.2 mosdepth=0.3.3 flye=2.9.1 nanofilt fastqc bedtools -c bioconda
    conda activate ont_tb
    # clone ONT-TB-NF
    git clone https://github.com/HKU-BAL/ONT-TB-NF.git
    cd ONT-TB-NF


Launch the pipeline execution with the following command: 

```
nextflow run_tb.nf --help
nextflow run_tb_amplicon.nf --help
```

### For adaptive sequencing

```
TB_NF_DIR={ONT-TB-NF PATH}
NF_S=${TB_NF_DIR}/run_tb.nf

SAMPLE_ID={NAME}
FQ={YOUR FQ FILE}
THREADS={THREAD}                  # threads number, e.g. 16
OUT_DIR={ABSOLUTE OUTPUT PATH}    # output path, abolute path required

nextflow run ${NF_S} \
--read_fq ${FQ} \
--sample_name ${SAMPLE_ID} \
--threads ${THREADS} \
--output_dir ${OUT_DIR}

```

### For amplicon sequencing

For Amplicon sequencing data, the user should provide the amplicon bed regions with `--amplicon_bed` option.

```
TB_NF_DIR={ONT-TB-NF PATH}
NF_S=${TB_NF_DIR}/run_tb_amplicon.nf

AMPLICON_BED={AMPLICON BED}       # your amplicon region 
FQ={YOUR FQ FILE}
SAMPLE_ID={NAME}
THREADS={THREAD}                  # threads number, e.g. 16
OUT_DIR={ABSOLUTE OUTPUT PATH}    # output path, abolute path required

nextflow run ${NF_S} \
--read_fq ${FQ} \
--sample_name ${SAMPLE_ID} \
--amplicon_bed ${GENE_BED} \
--threads ${THREADS} \
--output_dir ${OUT_DIR}

```

### For using Guppy for basecalling 

```

FAST5_DIR={Input FAST5 folders}
GUPPY_BASECALLER_PATH={Guppy basecaller path}                       # e.g. guppy_basecaller
GUPPY_CONFIG={Guppy config file path}                               # e.g. dna_r10.4_e8.1_sup.cfg

SAMPLE_ID={NAME}
TB_NF_DIR={ONT-TB-NF PATH}
NF_S=${TB_NF_DIR}/run_tb_amplicon.nf                                # e.g. run_tb.nf or run_tb_amplicon.nf
THREADS={THREAD}                                                    # threads number, e.g. 16
OUT_DIR={ABSOLUTE OUTPUT PATH}                                      # output path, abolute path required

nextflow run ${NF_S} --fast5_dir ${FAST5_DIR} --guppy_basecaller_path ${GUPPY_BASECALLER_PATH} --guppy_config_path ${GUPPY_CONFIG} --guppy_options "--device 'cuda:0'" --sample_name ${SAMPLE_ID} --threads ${THREADS} --output_dir ${OUT_DIR}

```

## Pipeline Summary

The ONT-TB-NF pipeline is dedecated for decting TB from ONT data. The data can be obtained from standard whole genome sequncing from MinION, from adaptive sequcing, or from amplicon sequcing (by amplify sepecifig regions in TB genome). 

For apply to WGS and adaptive sequcing (like from [readfish](https://www.nature.com/articles/s41587-020-00746-x) or [UNCALLED](https://www.nature.com/articles/s41587-020-0731-9)), please use the default `run_tb.nf` pipeline. 
For apply to Amplicon sequecing, please use the `run_tb_amplicon.nf` pipeline. 

In general, the ONT-TB-NF pipeline performs the following tasks:

- Basecalling (Guppy, guppy_basecaller)
- Sequncing quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Filtering and Trimming of read ([NanoFilt](https://github.com/wdecoster/nanofilt))
- Alignment ([minimap2](https://github.com/lh3/minimap2))
- Variant calling ([Clair3](https://github.com/HKU-BAL/Clair3))
- Antibiotic resistome finding ([TBProfiler](https://github.com/jodyphelan/TBProfiler))



## Pipeline results


Here is a brief description of output files created for each sample:
```
Basecalling results at:       {YOUR OUTPUR DIR}/0_bc
QC results at:                {YOUR OUTPUR DIR}/1_qc
Aligment results at:          {YOUR OUTPUR DIR}/2_aln
Variant calling results at:   {YOUR OUTPUR DIR}/3_vc
TB analysis report at:        {YOUR OUTPUR DIR}/4_tb
```

## Requirements 

* [Nextflow](https://www.nextflow.io) 22.04.5 (or later)
* Java 17 
* [Docker](https://www.docker.com/) 20.10.7 

