# Introduction

This repository contains a Snakemake workflow for performing a basic CUT&RUN analysis pipeline.

This workflow performs the following steps:

1.  Optional: Retrieval of publically available sequencing data from NCBI GEO/SRA.
2.  Optional: merge reads from the same sample that were sequenced on separate lanes or sequencing runs.
3.  Align raw reads to a reference genome. The workflow will automatically retrieve any publically available reference genome and build the bowtie2 index.
4.  Filter aligned reads based on alignment quality and discard multi-mapping reads. The filtering parameters can be customized.
5.  Optional: filter aligned reads to discard reads aligning to contigs or mitochondrial genome. This filtering can also be customized.
6.  Separate aligned reads based on fragment size.
7.  Optional: Perform peak calling using MACS2.
8.  Generate bigWig files containing z-score normalized read depth to allow data visualization in the genome browser.
9.  Optional: spike-in normalization for experiments where exogenous control chromatin or cells are added to each sample as an internal control. Spike-in normalization can also be performed using barcoded, recombinant spike-in nucleosomes bearing different histone modifications (e.g. EpiCypher SNAP-CUTANA).

# Running the workflow

Follow the steps below to run this workflow:

## Step 1: Software installation

Ensure that you have a conda-based Python3 distribution installed.
I recommend [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

Snakemake and Snakedeploy are used to deploy and run this workflow.
These can be installed using Mamba:

```         
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```

OR Conda:

```         
conda create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```

## Step 2: Deploy the workflow

Activate the conda environment:

```
conda activate snakemake
```

Create a new directory for your analysis and enter it:

```         
mkdir -p path/to/project-workdir
cd path/to/project-workdir
```

Deploy the workflow:

```         
snakedeploy deploy-workflow  https://github.com/tjgibson/NGS-workflow-CUTandRUN . --branch main 
```

This command will create all the files necessary for running this workflow.

## Step 3: Entering your sample information

The previous step should have created the file `config/units.tsv` in your local project directory.
This is a template that can be modified to contain the information about the samples you wish to analyze.
Before modifying it, take a look at the file to get a sense of how it is formatted.
The fields of this file are described below:

**sample_name:** The name you want to use for each individual sample.
I like to use some sort of naming convention like developmental-stage_genotype_treatment_rep1, where each piece of important information about the sample is separated by an underscore.
If a single piece of information, such as developmental stage, contains multiple words, separate them with a dash: L3-larva_rep1.

**unit_name:** If you sequenced the same sample across multiple lanes or sequencing runs, these will be listed as separate units.
For example, control_rep1_A and control_rep1_B.
If you only have a single fastq file for a sample, then the `unit_name` should be identical to the `sample_name`.

**fq1:** The file path to the location on your computer where the raw fastq file is stored.
These can be located anywhere on your computer, but I suggest placing them in a directory called `data`, `raw`, or `raw_data` within your analysis directory.
Fastq files should be compressed via gzip to save space.

**fq2:** The file path to the second fastq file.

**sra:** This workflow can perform automatic retrieval of publically available data from the Short Read Archive (SRA).
Enter the SRA accession number (e.g. SRR0000001) here.
If you provide both an SRA number and a local fastq file, the workflow will use the local file.

**read_format:** This should be set to `PE` for paired-end.
See note above regarding single-end data.

**call_peaks:** This field indicates if peak calling should be performed using MACS2 and should be set to `TRUE` or `FALSE`.

**input:** If you have control samples to be used during peak calling (e.g. IgG), this column will indicate which control samples should be used for which experimental samples.

**peak_type:** This specifies what type of peak calling will be performed by MACS2 and should be set to `narrow` (e.g. transcription factors, H3K27ac) or `broad` (e.g. H3K27me3, H3K9me3)

**sample_group:** Indicates the name to use for each group of replicates.
This name will be used to name bigWig and peak files for merged replicates.
For example, for samples control_rep1 and control_rep2, `sample_group` could be set to control.


## Step 4: Configuring the workflow

The workflow can be customized to suit your needs by editing the file `config/config.yaml`.
Below is a description of the different options that can be customized within the config file.

### Merging of technical replicates

If you have samples that were sequenced on multiple runs or sequencing lanes, change the following lines of the config.yaml file:

```         
mergeReads:
  activate: True
```

This will result in merging of the fastq files prior to alignment.

### Filtering of reads based on alignment to specific chromosomes

When working with Drosophila data (or other organisms with really good genome assemblies), I typically discard reads aligning to unplaced contigs or scaffolds.
For Drosophila, this entails only retaining reads aligning to the major chromosome arms (chr2L, chr2R, Chr3L, chr3R, chrX, and chrY).
This can be done with the following lines of the config file:

```         
filter_chroms: True
keep_chroms: config/keep_chroms.txt
```

You can modify the keep_chroms.txt file to reflect a different filtering strategy.
Make sure that the chromosome naming convention used in the file matches that of your chosen reference genome (see below).
For example, chr2L for UCSC vs. 2L for Ensembl.

### Setting the reference genome

This is where you set the reference genome to use for read alignment.
Give the genome a name and provide a link to the fasta file for the genome.
Reference genomes can be found through the [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html), [Ensembl](https://ensemblgenomes.org/), or organism-specific genome databases (e.g. [Flybase](https://ftp.flybase.net/genomes/)).
The pipeline will automatically retrieve the reference genome from the provided link and build the bowtie2 index required for alignment.
Modify the following lines as necessary for your desired reference genome:

```         
ref_genome:
  name: dm6
  link: https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
```

If you want to use a custom genome that is not publically available or a bowtie2 index that you have already built, create a folder called `resources/` in your working directory and copy your fasta file or the bowtie2 index files to this directory.
If using a custom fasta file, rename the file to be `ref_genome.fasta.gz`.
If using a custom bowtie2 index, use the prefix `genome` for the index files (e.g. `genome.1.bt2`).

### Setting the spike-in genome

This workflow can optionally perform spike-in normalization if your experiment included exogenous DNA that was added to each sample prior to immunoprecipitation. There are two options for spike-in normalization of CUT&RUN data:

1. Exogenous cells, chromatin or DNA from a different organism that are added to samples at some point during the experiment.
2. Barcoded, recombinant nucleosomes (e.g. EpiCypher SNAP-CUTANA spike-in)

For option 1, modify the following lines and provide a link to to a publically available genome for the spike-in organism:

```
# spike-in genome options ------------------------------------------------------------------
use_spikeIn: True
epiCypher_spikeIn: False
epiCypher_barcodes: config/epicypher_met_panel_nuc_barcodes.csv

spikeIn_genome:
  name: yeast_w303
  link: http://sgd-archive.yeastgenome.org/sequence/strains/W303/W303_SGD_2015_JRIU00000000/W303_SGD_2015_JRIU00000000.fsa.gz
```

For option 2, modify the lines as follows:

```
# spike-in genome options ------------------------------------------------------------------
use_spikeIn: True
epiCypher_spikeIn: True
epiCypher_barcodes: config/epicypher_met_panel_nuc_barcodes.csv
```

For option 2, you will need to provide a table describing each spike-in nucleosome and the DNA barcode associated with it. The pipeline will search for each barcode sequence in the raw fastq files and count the number of reads containing each barcode.

For either option, the percentage of spike-in reads will be used to calculate a scaling factor.
The values in the z-score normalized bigWig files will be multiplied by this scaling factor to produce spike-in normalized bigWig files.
The pipeline also outputs a table with the number of reference and spike-in reads, percentage of spike-in reads, and scaling factors.
If you wish to use a different strategy for spike-in normalization (e.g. skip initial z-score normalization), you should be able to use the information in this table to impliment alternate spike-in normalizations.


### modifying parameters for specific steps

The final section of the config file allows customization of the parameters for individual steps in the workflow:

```         
params:
  bowtie2_index: ""
  bowtie2_align: "-k 2 --very-sensitive --no-mixed --no-discordant -X 5000"
  filter_multireads: "-f bam -F '(mapping_quality >= 30) and ([XS] == null)'" 
  bigwigs_ind: "--binSize 10"
  bigwigs_merged: "--binSize 10"
  macs2_call_peaks_narrow: "-f BAMPE --keep-dup all -g dm --call-summits"
  macs2_call_peaks_broad: "-f BAMPE --keep-dup all -g dm"
  bigwigs_ind: "--binSize 10"
  bigwigs_merged: "--binSize 10"
  filter_peaks_by_replicates:
    min_overlap: 10
    replicate_threshold: 2

```

These can be modified as desired to alter these steps.
For example, by default the workflow will run bowtie2 with parameters `-k 2 --very-sensitive --no-mixed --no-discordant -X 5000`.
The parameters listed here in the config file will be passed directly to the corresponding step when it is run.
The parameters for `filter_peaks_by_replicates` indicates that peaks called by MACS2 will be filtered to only include peaks that were deteced in at least 2 replicates, considering a peak to be shared between replicates if the peak in each replicate overlap oneanother by at least 10 bp.

## Step 5: Running the workflow

You are now ready to run the workflow.
Before running the workflow, it is recommended to do a "dry run", in which Snakemake will check that all the required input and output files exist, and provide a summary of what the workflow will do.
To do a dry run, invoke the following command:

```         
snakemake -n
```

Examining the output of this dry run can help make sure that samples names and locations of raw files were entered correctly.

If everything looks good, run the workflow with the following command:

```         
snakemake --use-conda -c 1
```

The `-c 1` parameter tells snakemake to use a single core for running the workflow.
Depending on the number of processors/cores available on your computer, this number can be increased, which will speed up execution by allowing multiple steps to run in parallel and allowing steps such as read alignment to use multiple cores.

This workflow can also be run on a computing cluster or using cloud computing services.
See the relevant sections of the Snakemake documentation for more information: [Cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) [Cloud execution](https://snakemake.readthedocs.io/en/stable/executing/cloud.html)

# Navigating the output

Once finished, the workflow will have produced various output files that will be located in a folder named `results/`.
The results will be organized into the following subdirectories:

**aligned_reads:** BAM files containing aligned, sorted, filtered reads and BAI bam indices, separated by fragment size.
sample_small.bam will contain fragments \< 100bp, sample_large.bam will contain fragments \> 150bp, and sample_total.bam will contain all fragments.
A number of different intermediate BAM files will be generated as the workflow runs.
Upon completion, intermediate files will be deleted and only only the filtered, sorted BAM files will be retained.

**bigwigs:** z-score normalized bigWig files for individual replicates and merged replicates.
If spike-in normalization is indicated in the config file, spike-in normalized bigWigs will also appear here.

**peaks:** narrowPeak/broadPeak and bed files containing peaks called by MACS2.
This includes peaks called on individual replicates, merged replicates, and peaks filtered by replicates according to the criteria in the config file.

**scaling_factors:** Tables containing the number and percentage of spike-in reads, as well as the scaling factors. If a panel of barcoded nucleosomes is used as spike-in, an additional table will include the number or reads containing each barcode.

# Combining workflows for multiple data types

It may be desirable to process multiple data types (e.g. ChIP-seq, RNA-seq, ATAC-seq, CUT&RUN) and then integrate those results into downstream analysis.
This can be done in a single Snakemake workflow that contains multiple "modules" for the different data types.
For an example of how to build a complex workflow including multiple datasets as well as custom downstream analysis, see [this repository](https://github.com/tjgibson/S2_pioneers_manuscript) for an example.

# To-do

In the future, I hope to add some of the following features:

-   Add small, test fastq files to repo to allow for automated testing of rules

-   Create an HTML summary file containing read alignment and filtering statistics, number of peaks called, and other useful information about the workflow.

-   The optional step to take unaligned reads, BLAST them, and provide a summary of which organisms the unaligned reads likely originated from.

-   Add validation to ensure that unit table and config file are formatted correctly.
# smk-chipseq-pe-test
