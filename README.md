# VGEA
VGEA (Viral Genomes Easily Analyzed) is an RNA viral assembly toolkit.

VGEA was developed to aid in the analysis of next generation sequencing data. Users can do the following with this pipeline:

* Remove adapters, low quality bases/positions, and perform read-level QC 
* Align paired-end sequencing reads to the human reference genome.
* Extract unmapped/unaligned reads.
* Split bam files into forward and reverse reads. 
* Carry out de novo assembly of forward and reverse reads to generate contigs.
* Pre-process reads for quality and contamination. 
* Map reads to a reference tailored to the sample using corrected contigs supplemented by the user’s choice of reference sequences.
* Evaluate/assess the quality of genome assemblies.
* Collate results in a multiqc summary..

Dependencies: 

The VGEA pipeline requires the following dependencies:

* Snakemake (Köster et al., 2012)
* Fastp (Chen et al., 2018)
* BWA (Li and Durbin, 2009)
* Samtools (Li et al., 2009)
* IVA (Hunt et al., 2015)
* Shiver (Wymant et al., 2018)
* Seqtk ([Li, 2018](https://github.com/lh3/seqtk))
* Quast (Gurevich et al., 2013)
* Multiqc (Ewels et al., 2016)

VGEA was built on the Snakemake workflow management system and utilizes existing tools for each step: **fastp** (Chen et al., 2018) for read trimming and quality control, **bwa** (Li and Durbin, 2009) for mapping sequencing reads to the human reference genome, **samtools** (Li et al., 2009) for extracting unmapped reads and also for splitting bam files into fastq files, **iva** (Hunt et al., 2015) for de novo assembly to generate contigs, **shiver** (Wymant et al., 2018) to pre-process reads for quality and contamination, then map to a reference tailored to the sample using corrected contigs supplemented with the user’s choice of existing reference sequences and **quast** (Gurevich et al., 2013) to evaluate/assess the quality of genome assemblies. 
Finally, an interactive report of results are generated by **MultiQC** (Ewels et al., 2016).

## Installation 

### 1. Conda

This workflow requires conda to be installed and available on the system.

To do this install conda via the miniconda installers found [here](https://docs.conda.io/en/latest/miniconda.html) and instructions [here](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html).

Briefly: 

#### Linux

  To obtain the installer for linux use the following:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

  Then, install miniconda,

```
sh Miniconda3-latest-Linux-x86_64.sh
```

#### MacOS

  To obtain the installer for MacOS, you can [download](https://docs.conda.io/en/latest/miniconda.html) it manually or use wget:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

  Then, install miniconda,

```
sh Miniconda3-latest-MacOSX-x86_64.sh
```
### 2. Snakemake

Then install snakemake as follows:

    conda create -c bioconda -c conda-forge --name snakemake snakemake 

### 3. VGEA 

To complete installation of VGEA clone the directory and enter it.

```
git clone git@github.com:pauloluniyi/VGEA.git
cd VGEA
```

### 4. Data Dependency

Finally, you need to download a human reference genome. 
There is a convenience script provided to do this that can run as follows:

```
cd resources/
bash get_human_reference.sh 
cd ..
```

## Usage

There are 2 key files when running VGEA:

- A yaml config file containing key parameters e.g., `config/config.yaml`
 
- A 3 column tab-seperated sheet with a name and paths to r1 and r2 for each sample e.g., `./tests/integration/sample_table.tsv`

### Config

- `sample_table` is the path to the tsv created above with names and read locations for each sample (under headings id, r1, r2)

- `shiver_config_file` is the path to the config file for shiver, by default this is `config/shiver_config.sh`

- `human_reference_genome` is the path to the fasta file containing the human genome reference, if you used `get_human_reference.sh` this will be `resources/GRCh38_latest_genomic.fna`

- `viral_species` is one of the 4 pre-supported viral references (`HIV-1`, `SARS-CoV-2`, `LASV_L`, `LASV_S`) and is used to autofill the paths to various reference files used by default for these. If you supply your own viral reference/adapter/primer files then this doesn't need to be supplied

- `viral_reference_alignment` path to a reference alignment of the appropriate viral genomes for your samples (by default `resources/{viral_species}/MyRefAlignment.fasta`)

- `viral_reference_genome` path to a singular reference genome for your samples for QUAST based assembly assessment (by default `resources/{viral_species}/MyRefGenome.fasta`).

- `viral_reference_gene_features` path to a GFF3 containing gene features for the supplied viral reference genome (by default `resources/{viral_species}/MyRefFeatures.gff3`)

- `viral_sequencing_adapters` path to a fasta file containing the sequencing adapters used for your samples (by default `resources/{viral_species}/MyAdapters.fasta`)

- `viral_sequencing_primers` path to a fasta file containing the sequencing primers used for your samples (by default `resources/{viral_species}/MyRefAlignment.fasta`)

### Sample Table

A 3 column, tab-separated file with an *id* column containing sample names, *r1* with the path to the forward reads for that sample, and *r2* with the path to reverse reads for that sample.


    id	r1	r2
    test1	.tests/integration/test1_r1.fq.gz	.tests/integration/test1_r2.fq.gz
    test2	.tests/integration/test2_r1.fq.gz	.tests/integration/test2_r2.fq.gz

### Executing the workflow

The workflow will automatically install dependencies using conda/mamba if executed with `--use-conda` otherwise all dependencies listed at the start of the README must be manually installed via conda.

To run the workflow complete the above config and sample tables and execute:

```
snakemake --cores $n --use-conda --configfile $your_config
```

Where `$n` is the number of cores with which to execute the workflow and `$your_config` is the path to `config.yaml` you've created.

### Results

VGEA will output all results in the `results/` directory with a subfolder containing results for each sample and top-level folders collecting all log files, tool benchmarks, and an interactive multiQC html summary of results.

    results/
    ├── benchmarks                  # all benchmark files with hardware usage for each sample
    ├── logs                        # all log files for each rule and sample
    ├── multiqc                     # multiqc summary of quast, fastp, and de-hosting results
    │   ├── multiqc_data
    │   └── multiqc_report.html     # the interactive multiqc report
    ├── test1                       # all results for the sample "test1" 
    │   ├── test1_1.fastq           # fastp trimmed and de-hosted reads for the sample 
    │   ├── test1_2.fastq
    │   ├── test1.bam               # alignment against human reference 
    │   ├── test1.fasta             # final cleaned assembly from IVA and shiver
    │   ├── test1.fastp.html        # fastp report in html format
    │   ├── test1.fastp.json        # fastp report in json format
    │   ├── test1.flagstat          # dehosting mapping statistics
    │   ├── test1_iva               # folder containing all IVA assembly files
    │   ├── test1.quast_results     # folder containing all QUAST assembly assessment of test1.fast
    │   ├── test1_r1_trimmed.fq     # fastp trimmed reads
    │   ├── test1_r2_trimmed.fq     
    │   └── test1.sam               # sam file containing all the reads that didn't map to the human reference
    └── test2
        ├── shiver
        ├── test2_1.fastq
        ├── test2_2.fastq
        ├── test2.bam
        ├── test2.fasta
        ├── test2.fastp.html
        ├── test2.fastp.json
        ├── test2.flagstat
        ├── test2_iva
        ├── test2.quast_results
        ├── test2_r1_trimmed.fq
        ├── test2_r2_trimmed.fq
        └── test2.sam

## Containerized Singularity (Beta)

Alternatively, users can run the VGEA pipeline with all dependencies installed in a docker/singularity container.

This requires singularity and snakemake to be installed on the system but in theory provides a more reproducible version of the conda environments (not fully tested compared to just conda).

See [here](https://www.sylabs.io/docs/) for instructions to install Singularity.

Then the workflow can be run as normal with `--use-singularity` added e.g.,

```
snakemake --use-conda --use-singularity --configfile .test/integration/test_config.yaml -j 1
```

## Testing

To run a minimal integration test once snakemake and conda are installed:

```
snakemake --use-conda --configfile .test/integration/test_config.yaml -j 1
```

# Contributions

- Paul Eniola Oluniyi
- Gerry Tonkin-Hill (https://gtonkinhill.github.io/)
- Finlay Maguire (https://finlaymagui.re)
