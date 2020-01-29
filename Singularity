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

    # make sure the "pipeline/scripts" folder exists
    mkdir -p ${SINGULARITY_ROOTFS}/pipeline/scripts
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
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.10.0+-x64-linux.tar.gz
    echo 'PATH=$PATH:~/ncbi-blast-2.10.0+/bin/' >> $SINGULARITY_ENVIRONMENT

    # Install samtools
    conda install -c bioconda samtools

    # Install mummer
    conda install -c bioconda mummer

    # Install fastp
    git clone https://github.com/OpenGene/fastp.git
    cd fastp
    make && make install
    cd ..

    # Install mafft
    wget https://mafft.cbrc.jp/alignment/software/mafft-7.453-without-extensions-src.tgz
    tar -xzf mafft-7.453-without-extensions-src.tgz
    cd mafft-7.453-without-extensions/core/
    make clean
    make
    make install
    cd ~

    # Install java
    yes | apt-get install default-jdk

    # Install smalt
    wget https://sourceforge.net/projects/smalt/files/latest/download -O smalt.tgz
    tar -xzf smalt.tgz
    cd smalt-0.7.6/
    ./configure
    make
    make install
    cd ~

    # Install bwa
    git clone https://github.com/lh3/bwa.git
    cd bwa
    make
    echo 'PATH=$PATH:~/bwa/' >> $SINGULARITY_ENVIRONMENT
    
    # Install bowtie
    wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip/download -O bowtie2.zip 
    unzip bowtie2.zip
    echo 'PATH=$PATH:~/bowtie2-2.3.5.1-linux-x86_64/' >> $SINGULARITY_ENVIRONMENT

    # Install kmc 
    wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz
    tar -xzf KMC3.linux.tar
    echo 'PATH=$PATH:~/kmc/' >> $SINGULARITY_ENVIRONMENT
    echo 'PATH=$PATH:~/kmc_dump/' >> $SINGULARITY_ENVIRONMENT
    echo 'PATH=$PATH:~/kmc_tools/' >> $SINGULARITY_ENVIRONMENT

    # Install iva
    conda install -c bioconda iva

    
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

    # Copy data files
    data/OWOIF.bam
    data/MyRefAlignment.fasta
    data/MyAdapters.fasta
    data/MyPrimers.fasta

    # Copy pipeline scripts
    pipeline/scripts /pipeline
    pipeline/Snakefile /pipeline

%runscript
  
    snakemake -q -s /pipeline/Snakefile
    
