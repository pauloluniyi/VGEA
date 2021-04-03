Bootstrap: debootstrap
OSVersion: xenial
MirrorURL: http://us.archive.ubuntu.com/ubuntu/
IncludeCmd: yes
Include: bash vim less man-db apt-utils tzdata

%help

    This is the container for the VGEA pipeline

%labels

    AUTHOR - Paul Eniola Oluniyi - pauleniolaoluniyi@gmail.com
    PhD researcher at ACEGID, Redeemer's University
    Version - v1.0

%setup

    # make sure the "pipeline/" folder exists
    mkdir -p ${SINGULARITY_ROOTFS}/pipeline

    # make sure the "/data" folder exists and in this folder, place your input files
    mkdir ${SINGULARITY_ROOTFS}/data

%environment

    PATH="/usr/local/anaconda/bin:$PATH"

%post

    echo "The post section is where you can install, and configure the container"

    mv /etc/apt/sources.list /etc/apt/sources.list.bak

    echo "deb http://us.archive.ubuntu.com/ubuntu/ xenial main restricted universe multiverse
    deb-src http://us.archive.ubuntu.com/ubuntu/ xenial main restricted universe multiverse
    deb http://us.archive.ubuntu.com/ubuntu/ xenial-security main restricted universe multiverse
    deb http://us.archive.ubuntu.com/ubuntu/ xenial-updates main restricted universe multiverse
    deb http://us.archive.ubuntu.com/ubuntu/ xenial-proposed main restricted universe multiverse
    deb http://us.archive.ubuntu.com/ubuntu/ xenial-backports main restricted universe multiverse
    deb-src http://us.archive.ubuntu.com/ubuntu/ xenial-security main restricted universe multiverse
    deb-src http://us.archive.ubuntu.com/ubuntu/ xenial-updates main restricted universe multiverse
    deb-src http://us.archive.ubuntu.com/ubuntu/ xenial-proposed main restricted universe multiverse
    deb-src http://us.archive.ubuntu.com/ubuntu/ xenial-backports main restricted universe multiverse" >> /etc/apt/sources.list


    # Install wget, git, make
    apt-get -y --force-yes update
    yes | apt-get install build-essential
    yes | apt-get install git
    yes | apt-get install wget
    yes | apt-get install autoconf autogen libtool
    yes | apt install bc

    # Install zlib, cmake
    yes | apt-get install zlib1g-dev libbz2-dev pkg-config cmake

    # Install curl, lzma
    yes | apt-get install curl liblzma-dev libncurses5-dev


    # Install Anaconda3
    if [ ! -d /usr/local/anaconda ]; then
         wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
            -O ~/anaconda.sh && \
         bash ~/anaconda.sh -b -p /usr/local/anaconda && \
         rm ~/anaconda.sh
    fi
    # set anaconda path
    export PATH="/usr/local/anaconda/bin:$PATH"
    
    # Install Python2
    conda install -c anaconda python2
    
    # Install Python3
    conda install -c anaconda python3

    # Install snakemake
    conda install -c bioconda -c conda-forge snakemake

    # Install fastaq
    yes | apt-get install python3-pip
    pip3 install pyfastaq

    # Install biopython
    conda install -c conda-forge biopython

    # Install blast
    conda install -c bioconda blast

    # Install samtools
    conda install -c bioconda samtools

    # Install mummer
    conda install -c bioconda mummer

    # Install fastp
    conda install -c bioconda fastp

    # Install mafft
    conda install -c bioconda mafft

    # Install smalt
    conda install -c bioconda smalt

    # Install bwa
    conda install -c bioconda bwa

    # Install unzip
    yes | apt-get install -y unzip zip

    # Install bowtie
    conda install -c bioconda bowtie2

    # Install kmc
    conda install -c bioconda kmc
    
    # Install Picard
    conda install -c bioconda picard
    
    # Install shiver
    conda install -c bioconda shiver
    
    # Install quast
    conda install -c bioconda quast

    # Install iva
    conda install -c bioconda iva

    # Install java
    yes | apt-get install default-jdk

    # add read/write access to data folder
    chmod -R o+rwx /data

%files

    # Copy input files
    {id}_1.fastq /data
    {id}_2.fastq /data
    MyRefAlignment.fasta /data
    MyAdapters.fasta /data
    MyPrimers.fasta /data
    quast_refseq.fasta /data
    quast_genefeatures.txt /data


    # Copy pipeline scripts
    
    Snakefile /pipeline

%runscript

    snakemake -q -s /pipeline/Snakefile
