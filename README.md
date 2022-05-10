# wi-STRs

A Nextflow pipeline used for STR variant calling.

## Execution of pipeline using Nextflow
```
git clone https://github.com/AndersenLab/wi-STRs.git

cd wi-STRs

nextflow STR_variants_calling.nf --ref_fa=bin/c_elegans.PRJNA13758.WS276.genomic.fa --ref_str=build_ref/ref_ce.hipstr_reference.bed --bamfiles=bam_test.txt 

```

## Required software packages that should be in users PATH

1. [nextflow-v19.07.0](https://www.nextflow.io/docs/latest/getstarted.html)
2. [BCFtools-v1.9](https://samtools.github.io/bcftools/bcftools.html)
3. [HipSTR-v0.7](https://github.com/tfwillems/HipSTR)
4. [bedtools-v2.29.2](https://bedtools.readthedocs.io/en/latest/)

## Pipeliine parameters

* --ref_fa

Reference genome in FASTA format from WormBase https://wormbase.org/

* --ref_str

Reference STR in BED format

* --bamfiles

We use a sample sheet as the input file, including full path to bam files, see the example in `bam_test.txt`.

* --email

Add your email with command

* --out

Add result folder name. The default is "STR_Results-*date*"

## Output

This pipeline will generate three vcf files:

STR_all_raw.vcf.gz :              The raw vcf includes all the STR variants across samples.

STR_all_filtered.vcf.gz :         The filitered vcf from STR_all_raw.vcf.gz using the script [HipSTR filter_vcf.py](https://github.com/tfwillems/HipSTR/blob/master/scripts/filter_vcf.py) and recommended settings.

STR_all_filtered_Fmiss01.vcf.gz:  The final vcf from STR_all_filtered.vcf.gz by filtering STR variants with equal or more than 10% missing data across all samples using BCFtools.


## Building the STR reference

STR reference can be built using the script build_ref/STRref.sh

```
cd build_ref

STRref.sh <prefix of output> <reference genome in FASTA format>

```

The script STRref.sh is built using [scripts](https://github.com/HipSTR-Tool/HipSTR-references/tree/master/scripts) and the [framework](https://github.com/HipSTR-Tool/HipSTR-references/blob/master/mouse/mouse_reference.md) on [HipSTR-references](https://github.com/HipSTR-Tool/HipSTR-references).

### Required software package

1. [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html)

