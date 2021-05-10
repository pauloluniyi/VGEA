rule assembly:
    message:
        "Assembly of forward and reverse reads"
    conda:
        "../envs/iva.yaml"
    # container: "docker://quay.io/biocontainers/iva:1.0.11--py_0"
    input:
        forward_read=rules.bamtoFastq.output.forward_read,
        reverse_read=rules.bamtoFastq.output.reverse_read,
    output:
        contigs="results/{id}/{id}_iva/contigs.fasta",
    threads: 8
    shell:
        "iva --reads_fwd {input.forward_read} --reads_rev {input.reverse_read} --threads {threads} {output}"


rule shiver_init:
    message:
        "Shiver initialization"
    conda:
        "../envs/shiver.yaml"
    # container: "docker://quay.io/biocontainers/shiver:1.3.5--py27_0"
    input:
        reference_alignment=config["viral_reference_genome"],
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
        shiver_init.sh {params.shiver_init_dir_path} {input.shiver_config} {input.reference_alignment} {input.adapters} {input.primers}
        cd ../../../
        """


rule align_contigs:
    message:
        "Aligning contigs"
    conda:
        "../envs/shiver.yaml"
    # container: "docker://quay.io/biocontainers/shiver:1.3.5--py27_0"
    input:
        initialization_directory=rules.shiver_init.output.initialization_directory,
        contigs_file=rules.assembly.output.contigs,
        shiver_config=config["shiver_config_file"],
    output:
        blast_hits="results/{id}/shiver/{id}.blast",
        aligned_contigs_raw="results/{id}/shiver/{id}_raw_wRefs.fasta",
        aligned_contigs_cut="results/{id}/shiver/{id}_cut_wRefs.fasta",
    shell:
        """
        cd results/{wildcards.id}/shiver
        shiver_align_contigs.sh {input.initialization_directory} {input.shiver_config} ../../../{input.contigs_file} {wildcards.id}
        cd ../../..
        """


rule map:
    message:
        "Mapping paired-end reads to reference genome"
    conda:
        "../envs/shiver.yaml"
    # container: "docker://quay.io/biocontainers/shiver:1.3.5--py27_0"
    input:
        initialization_directory=rules.shiver_init.output.initialization_directory,
        contigs=rules.assembly.output.contigs,
        blast_hits=rules.align_contigs.output.blast_hits,
        aligned_contigs_cut=rules.align_contigs.output.aligned_contigs_cut,
        forward_read=rules.bamtoFastq.output.forward_read,
        reverse_read=rules.bamtoFastq.output.reverse_read,
        shiver_config=config["shiver_config_file"],
    output:
        ref_seqs="results/{id}/shiver/{id}_ref.fasta",
        base_freqs="results/{id}/shiver/{id}_BaseFreqs.csv",
        base_freqs_global_aln="results/{id}/shiver/{id}_BaseFreqs_ForGlobalAln.csv",
        coords="results/{id}/shiver/{id}_coords.csv",
        insert_size_dist="results/{id}/shiver/{id}_InsertSizeCounts.csv",
        consensus_genome="results/{id}/shiver/{id}_remap_consensus_MinCov_10_30.fasta",
    shell:
        """
        cd results/{wildcards.id}/shiver
        shiver_map_reads.sh {input.initialization_directory} {input.shiver_config} ../../../{input.contigs} {wildcards.id} ../../../{input.blast_hits} ../../../{input.aligned_contigs_cut} ../../../{input.forward_read} ../../../{input.reverse_read}
        cd ../../..
        """
