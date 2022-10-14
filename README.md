# ONT-TB-NF
A Nextflow pipeline for *Mycobacterium tuberculosis* (TB) antibiotic resistance gene analysis with ONT data.


## Pipeline Description


Comprehensive pipeline for detection of TB from ONT adaptive sequencing and amplicon data. 
Pipeline including quality control, alignment, variant calling, consensus generation and antimicrobial resistance prediction. 

## Quickstart 

Install [Nextflow](https://www.nextflow.io/) by using the following command: 

    curl -s https://get.nextflow.io | bash 
    
Install required packages and activate envirment with [conda](https://conda.io/projects/conda/en/latest/index.html):
    
    conda create -n ont_tb samtools=1.15.1 minimap2=2.24 nanoplot=1.40.2 mosdepth=0.3.3 flye=2.9.1 nanofilt fastqc bedtools -c bioconda
    conda activate ont_tb


Launch the pipeline execution with the following command: 


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

## Pipeline Summary

The ONT-TB-NF pipeline is dedecated for decting TB from ONT data. The data can be obtained from standard whole genome sequncing from MinION, from adaptive sequcing, or from amplicon sequcing (by amplify sepecifig regions in TB genome). For apply to WGS and adaptive sequcing (like from readfish or UNCALLED), please use the default `run_tb.nf` pipeline. For apply to amplicon sequecing, please use the `run_tb_amplicon.nf' pipeline. In general, the ONT-TB-NF pipeline performs the following tasks:

- Sequncing quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Filtering and Trimming of read ([NanoFilt](https://github.com/wdecoster/nanofilt))
- Alignment ([[minimap2](https://github.com/lh3/minimap2))
- Variant calling ([Clair3](https://github.com/HKU-BAL/Clair3))
- Antibiotic resistome prection ([TBProfiler](https://github.com/jodyphelan/TBProfiler))


--nanofilt_options "STR"    Define the NanoFilt option for sequnencing filterinf and trimming, e.g. "-q 7" for filtering average quality low than 7.
--

## Pipeline results


Here is a brief description of output files created for each sample:
```
QC results at:                {YOUR OUTPUR DIR}/1_qc
Aligment results at:          {YOUR OUTPUR DIR}/2_aln
Variant calling results at:   {YOUR OUTPUR DIR}/3_vc
TB analysis report at:        {YOUR OUTPUR DIR}/4_tb
```

## Requirements 

* [Nextflow](https://www.nextflow.io) 22.04.5 (or later)
* Java 17 
* [Docker](https://www.docker.com/) 20.10.7 

