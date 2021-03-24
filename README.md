# VGEA
VGEA (Viral Genomes Easily Assembled) is a pipeline for advanced assembly of RNA viral genomes from next-generation sequencing data

VGEA was developed to aid in the analysis of next generation sequencing data. Users can do the following with this pipeline:

❖ Split bam or fastq files into forward and reverse reads. ❖ Carry out de novo assembly of forward and reverse reads to generate contigs. ❖ Pre-process reads for quality and contamination. ❖ Map reads to a reference tailored to the sample using corrected contigs supplemented by the user’s choice of reference sequences.

Dependencies: 

The VGEA pipeline requires the following dependencies:

* Python 3 (www.python.org)
* Snakemake (Köster et al., 2012)
* Samtools (Li et al., 2009)
* IVA (Hunt et al., 2015)
* Shiver (Wymant et al., 2018)

This pipeline was built on the snakemake workflow management system (Koster and Rahmann 2012). Several tools are used to perform different tasks within the pipeline: Samtools (Li et al., 2009) for splitting of bam or fastq files into forward and reverse reads; IVA (Hunt et al., 2015) for de novo assembly to generate contigs; Shiver (Wymant et al., 2018) to pre-process reads for quality and contamination, then map to a reference tailored to the sample using corrected contigs supplemented with the user’s choice of existing reference sequences

# Conda

Create a `conda` environment and install the dependencies

```
conda create -n vgea
conda activate vgea
conda install -c bioconda iva
conda install -c bioconda shiver
conda install -c bioconda snakemake
```

Clone the directory

```
git clone git@github.com:pauloluniyi/VGEA.git
```

Run the snakemake pipeline assuming the input bam file is in `/workingdir/`

```
snakemake -d /workingdir/
```

The location of the reference, primer and adapter fasta files can be controlled by adjusting the `config.yaml` in the VGEA working directory.

# Singularity

Alternatively, users can run the entire VGEA pipeline with all dependencies installed from the singularity container.

To install Singularity: See https://www.sylabs.io/docs/ for instructions to install Singularity.

* Building the container

If you are the administrator on your machine, you can build a local image of the VGEA container using the Singularity recipe file provided:

sudo singularity build vgea.simg Singularity


A sample singularity image built using the test data provided is available at 10.5281/zenodo.3702287.

