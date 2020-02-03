Bootstrap: debootstrap
OSVersion: xenial
MirrorURL: http://us.archive.ubuntu.com/ubuntu/
IncludeCmd: yes
Include: bash vim less man-db apt-utils tzdata

%help

    This is the container for the VGEA pipeline

%labels

    AUTHOR Paul Oluniyi - oluniyip@run.edu.ng
    PhD researcher at ACEGID, Redeemer's University
    Version - v1.0.0

%setup

    # make sure the "pipeline/" folder exists
    mkdir -p ${SINGULARITY_ROOTFS}/pipeline
    # make sure the "/data" folder exists
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


    # Install snakemake
    conda install -c bioconda -c conda-forge snakemake

    # Install fastaq
    pip install fastaq

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

    # Install iva
    conda install -c bioconda iva

    # Install java
    yes | apt-get install default-jdk

    # Install Python2
    yes | apt install python-pip
    yes | apt install python2.7

    # Add the pipeline scripts to the $PATH and make sure they are executable
    chmod +x /pipeline/scripts/*.sh
    chmod +x /pipeline/scripts/*.jar
    chmod +x /pipeline/scripts/tools/*.py
    chmod +x /pipeline/scripts/tools/*.pyc
    echo 'export PATH=$PATH:/pipeline/scripts' >> $SINGULARITY_ENVIRONMENT
    echo 'export PATH=$PATH:/pipeline/scripts/tools' >> $SINGULARITY_ENVIRONMENT

    # add read/write access to data folder
    chmod -R o+rwx /data

%files

    # Copy test data files available at (https://drive.google.com/drive/folders/1Z-kMWj9fOEEu0lFjw_F4EBLl73JqklLE?usp=sharing)
    934.bam /data
    MyRefAlignment.fasta /data
    MyAdapters.fasta /data
    MyPrimers.fasta /data

    # Copy pipeline scripts
    scripts/ /pipeline
    Snakefile /pipeline

%runscript

    snakemake -q -s /pipeline/Snakefile
