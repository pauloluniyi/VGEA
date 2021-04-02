# VGEA
VGEA (Viral Genomes Easily Assembled) is an RNA viral assembly toolkit.

VGEA was developed to aid in the analysis of next generation sequencing data. Users can do the following with this pipeline:

* Merge fastq files into unaligned bam files.
* Split bam files into forward and reverse reads. 
* Carry out de novo assembly of forward and reverse reads to generate contigs.
* Pre-process reads for quality and contamination. 
* Map reads to a reference tailored to the sample using corrected contigs supplemented by the user’s choice of reference sequences.
* Evaluate/assess the quality of genome assemblies.

Dependencies: 

The VGEA pipeline requires the following dependencies:

* Python 3 (www.python.org)
* Picard (https://github.com/broadinstitute/picard)
* Snakemake (Köster et al., 2012)
* Samtools (Li et al., 2009)
* IVA (Hunt et al., 2015)
* Shiver (Wymant et al., 2018)
* Quast (Gurevich et al., 2013)

VGEA was built on the Snakemake workflow management system and utilizes existing tools for each step: **picard**(https://github.com/broadinstitute/picard) for converting fastq files to a bam file, **samtools** (Li et al., 2009) for splitting files, **iva** (Hunt et al., 2015) for de novo assembly to generate contigs, **shiver** (Wymant et al., 2018) to pre-process reads for quality and contamination, then map to a reference tailored to the sample using corrected contigs supplemented with the user’s choice of existing reference sequences and **quast** (Gurevich et al., 2013) to evaluate/assess the quality of genome assemblies.

## Quick start

Instructions to type in a shell

1. [Install](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) miniconda3

# Linux

  To obtain the installer for linux use the following:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

  Then, install miniconda,

```
sh Miniconda3-latest-Linux-x86_64.sh
```

# MacOS

  To obtain the installer for MacOS, you can [download](https://docs.conda.io/en/latest/miniconda.html) it manually or use wget:
```
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

  Then, install miniconda,

```
sh Miniconda3-latest-MacOSX-x86_64.sh
```

2. Create conda virtual environment

# Conda

For detailed installation of **conda**, follow the instructions in the following link: https://conda.io/projects/conda/en/latest/user-guide/install/index.html.

Following installation of conda; 

Create a `conda` environment and install the dependencies

```
conda create -n vgea
conda activate vgea
conda install -c bioconda picard
conda install -c bioconda iva
conda install -c bioconda shiver
conda install -c bioconda snakemake
conda install -c bioconda quast
```

Clone the directory

```
git clone git@github.com:pauloluniyi/VGEA.git
```

Run the snakemake pipeline assuming the input fastq files are in `/workingdir/`

```
snakemake -d /workingdir/
```

The location of the reference, primer and adapter fasta files for **shiver** and the reference and gene features files for **quast** can be controlled by adjusting the `config.yaml` in the VGEA working directory.

# Singularity

Alternatively, users can run the entire VGEA pipeline with all dependencies installed from the singularity container.

To install Singularity: See https://www.sylabs.io/docs/ for instructions to install Singularity.

* Building the container

If you are the administrator on your machine, you can build a local image of the VGEA container using the Singularity recipe file provided:

sudo singularity build vgea.simg Singularity

# Contributions

- Paul Eniola Oluniyi
- Gerry Tonkin-Hill (https://gtonkinhill.github.io/)
