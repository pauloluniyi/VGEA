rule iva_assembly:
    message:
        "IVA assembly of de-hosted trimmed reads: {wildcards.id}"
    conda:
        "../envs/iva.yaml"
    # container: "docker://quay.io/biocontainers/iva:1.0.11--py_0"
    threads: 4
    log:
        "results/logs/{id}/{id}_iva.log"
    benchmark:
        "results/benchmarks/iva/{id}.tsv"
    input:
        forward_read=rules.generate_dehosted_fastq.output.forward_read,
        reverse_read=rules.generate_dehosted_fastq.output.reverse_read,
    output:
        contigs="results/{id}/{id}_iva/contigs.fasta"
    params:
        output_folder="results/{id}/{id}_iva"
    shell:
        """
        rm -rf {params.output_folder} #to prevent snakemake pre-making the folder
        (iva --reads_fwd {input.forward_read} --reads_rev {input.reverse_read} --threads {threads} {params.output_folder}) > {log} 2>&1
        """


rule shiver_initialization:
    message:
        "Shiver run preparation: {wildcards.id}"
    conda:
        "../envs/shiver.yaml"
    # container: "docker://quay.io/biocontainers/shiver:1.3.5--py27_0"
    log:
        "results/logs/{id}/{id}_shiver_init.log"
    benchmark:
        "results/benchmarks/shiver_init/{id}.tsv"
    input:
        reference_alignment=config["viral_reference_alignment"],
        adapters=config["viral_sequencing_adapters"],
        primers=config["viral_sequencing_primers"],
        shiver_config=config["shiver_config_file"],
    output:
        initialization_directory=directory("results/{id}/shiver/{id}_shiver_init_dir"),
    params:
        shiver_init_dir_path="{id}_shiver_init_dir",
    shell:
        """
        cd results/{wildcards.id}/shiver
        (shiver_init.sh {params.shiver_init_dir_path} {input.shiver_config} {input.reference_alignment} {input.adapters} {input.primers}) > ../../../{log} 2>&1
        cd ../../../
        """


rule shiver_align_contigs_to_reference:
    message:
        "Shiver alignment of contigs to reference: {wildcards.id}"
    conda:
        "../envs/shiver.yaml"
    # container: "docker://quay.io/biocontainers/shiver:1.3.5--py27_0"
    log:
        "results/logs/{id}/{id}_shiver_align_contigs.log"
    benchmark:
        "results/benchmarks/shiver_align_contigs/{id}.tsv"
    input:
        initialization_directory=rules.shiver_initialization.output.initialization_directory,
        contigs_file=rules.iva_assembly.output.contigs,
        shiver_config=config["shiver_config_file"],
    output:
        blast_hits="results/{id}/shiver/{id}.blast",
        aligned_contigs_raw="results/{id}/shiver/{id}_raw_wRefs.fasta",
        aligned_contigs_cut="results/{id}/shiver/{id}_cut_wRefs.fasta",
    params:
        shiver_init_dir_path="{id}_shiver_init_dir",
    shell:
        """
        cd results/{wildcards.id}/shiver
        (shiver_align_contigs.sh {params.shiver_init_dir_path} {input.shiver_config} ../../../{input.contigs_file} {wildcards.id})  > ../../../{log} 2>&1
        cd ../../..
        """


rule shiver_read_remapping:
    message:
        "Shiver mapping reads to reference alignment: {wildcards.id}"
    conda:
        "../envs/shiver.yaml"
    # container: "docker://quay.io/biocontainers/shiver:1.3.5--py27_0"
    log:
        "results/logs/{id}/{id}_shiver_map_reads.log"
    benchmark:
        "results/benchmarks/shiver_map_reads/{id}.tsv"
    input:
        initialization_directory=rules.shiver_initialization.output.initialization_directory,
        contigs=rules.iva_assembly.output.contigs,
        blast_hits=rules.shiver_align_contigs_to_reference.output.blast_hits,
        aligned_contigs_cut=rules.shiver_align_contigs_to_reference.output.aligned_contigs_cut,
        forward_read=rules.generate_dehosted_fastq.output.forward_read,
        reverse_read=rules.generate_dehosted_fastq.output.reverse_read,
        shiver_config=config["shiver_config_file"],
    output:
        ref_seqs="results/{id}/shiver/{id}_ref.fasta",
        base_freqs="results/{id}/shiver/{id}_BaseFreqs.csv",
        base_freqs_global_aln="results/{id}/shiver/{id}_BaseFreqs_ForGlobalAln.csv",
        coords="results/{id}/shiver/{id}_coords.csv",
        insert_size_dist="results/{id}/shiver/{id}_InsertSizeCounts.csv",
        consensus_genome="results/{id}/shiver/{id}_remap_consensus_MinCov_15_30.fasta",
    params:
        shiver_init_dir_path="{id}_shiver_init_dir",
    shell:
        """
        cd results/{wildcards.id}/shiver
        (shiver_map_reads.sh {params.shiver_init_dir_path} {input.shiver_config} ../../../{input.contigs} {wildcards.id} ../../../{input.blast_hits} ../../../{input.aligned_contigs_cut} ../../../{input.forward_read} ../../../{input.reverse_read}) > ../../../{log} 2>&1
        cd ../../..
        """

rule tidy_shiver_output:
    message:
        "Cleaning shiver assembly for quast: {wildcards.id}"
    conda:
        "../envs/seqtk.yaml"
    log:
        "results/logs/{id}/{id}_shiver_tidy_contigs.log"
    benchmark:
        "results/benchmarks/shiver_tidy/{id}.tsv"
    input:
        consensus_genome="results/{id}/shiver/{id}_remap_consensus_MinCov_15_30.fasta",
    output:
        "results/{id}/{id}.shiver_assembly.fasta"
    shell:
        """
        seqtk seq -l0 {input} | head -n2 | sed '/>/!s/-//g' | sed 's/\\?/N/g' | sed 's/_remap_consensus//g' | seqtk seq -l80 > {output} 2> {log}
        """
