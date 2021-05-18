FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="5d84f5f70265dec1e783e19e1fdda3f098b3871d4f4d9916055c64e434a088ca"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/iva.yaml
#   prefix: /conda-envs/32d984c246da8df3871b0ad2191ea64b
#   name: iva
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - iva=1.0.11
RUN mkdir -p /conda-envs/32d984c246da8df3871b0ad2191ea64b
COPY workflow/envs/iva.yaml /conda-envs/32d984c246da8df3871b0ad2191ea64b/environment.yaml

# Conda environment:
#   source: workflow/envs/qc.yaml
#   prefix: /conda-envs/b1e51184c6742d9c5f7d3643dddf2d62
#   name: qc
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - multiqc=1.10.1
#     - quast=5.0.2
RUN mkdir -p /conda-envs/b1e51184c6742d9c5f7d3643dddf2d62
COPY workflow/envs/qc.yaml /conda-envs/b1e51184c6742d9c5f7d3643dddf2d62/environment.yaml

# Conda environment:
#   source: workflow/envs/read_processing.yaml
#   prefix: /conda-envs/b886d47ebcea709153b6f48c6b5f74a5
#   name: read_processing
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bwa=0.7.17
#     - fastp=0.20.1
#     - samtools=1.12
RUN mkdir -p /conda-envs/b886d47ebcea709153b6f48c6b5f74a5
COPY workflow/envs/read_processing.yaml /conda-envs/b886d47ebcea709153b6f48c6b5f74a5/environment.yaml

# Conda environment:
#   source: workflow/envs/seqtk.yaml
#   prefix: /conda-envs/6f11c32683b75717896b0694bdeeacbf
#   name: seqtk
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - seqtk=1.3
RUN mkdir -p /conda-envs/6f11c32683b75717896b0694bdeeacbf
COPY workflow/envs/seqtk.yaml /conda-envs/6f11c32683b75717896b0694bdeeacbf/environment.yaml

# Conda environment:
#   source: workflow/envs/shiver.yaml
#   prefix: /conda-envs/9bc42145b3a91c4c4a7d40b471672a9c
#   name: shiver
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - shiver=1.3.5
#     - bc
RUN mkdir -p /conda-envs/9bc42145b3a91c4c4a7d40b471672a9c
COPY workflow/envs/shiver.yaml /conda-envs/9bc42145b3a91c4c4a7d40b471672a9c/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/32d984c246da8df3871b0ad2191ea64b --file /conda-envs/32d984c246da8df3871b0ad2191ea64b/environment.yaml && \
    mamba env create --prefix /conda-envs/b1e51184c6742d9c5f7d3643dddf2d62 --file /conda-envs/b1e51184c6742d9c5f7d3643dddf2d62/environment.yaml && \
    mamba env create --prefix /conda-envs/b886d47ebcea709153b6f48c6b5f74a5 --file /conda-envs/b886d47ebcea709153b6f48c6b5f74a5/environment.yaml && \
    mamba env create --prefix /conda-envs/6f11c32683b75717896b0694bdeeacbf --file /conda-envs/6f11c32683b75717896b0694bdeeacbf/environment.yaml && \
    mamba env create --prefix /conda-envs/9bc42145b3a91c4c4a7d40b471672a9c --file /conda-envs/9bc42145b3a91c4c4a7d40b471672a9c/environment.yaml && \
    mamba clean --all -y
