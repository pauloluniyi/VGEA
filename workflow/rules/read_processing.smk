rule human_reference_index_creation:
    message:
        "Indexing the human reference genome"
    conda:
        "../envs/read_processing.yaml"
    # container: "docker://quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
    log:
        "results/logs/bwa_human_reference_index.log"
    benchmark:
        "results/benchmarks/human_reference_index.tsv"
    input:
        human_ref_genome=config["human_reference_genome"],
    output:
        multiext(
            config["human_reference_genome"], ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    shell:
        "bwa index {input} > {log} 2>&1" 


rule fastp_read_qc:
    message:
        "Trimming and performing read-level QC using FASTP: {wildcards.id}"
    conda:
        "../envs/read_processing.yaml"
    log:
        "results/logs/{id}/{id}_fastp.log"
    benchmark:
        "results/benchmarks/fastp/{id}.tsv"
    input:
        unpack(get_fastq)
    output:
        r1_trimmed="results/{id}/{id}_r1_trimmed.fq",
        r2_trimmed="results/{id}/{id}_r2_trimmed.fq",
        json_report="results/{id}/{id}.fastp.json",
        html_report="results/{id}/{id}.fastp.html"
    threads: 4
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1_trimmed} --out2 {output.r2_trimmed} --report_title {wildcards.id} --thread {threads} --json {output.json_report} --html {output.html_report} > {log} 2>&1" 


rule map_reads_to_human_genome:
    message:
        "Mapping reads to the human genome to remove human contaminants: {wildcards.id}"
    conda:
        "../envs/read_processing.yaml"
    # container: "docker://quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
    log:
        "results/logs/{id}/{id}_bwa_human_contamination_map.log"
    benchmark:
        "results/benchmarks/bwa_human/{id}.tsv"
    input:
        r1=rules.fastp_read_qc.output.r1_trimmed,
        r2=rules.fastp_read_qc.output.r2_trimmed,
        idx=rules.human_reference_index_creation.output,
    output:
        mapped_bam="results/{id}/{id}.sam",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
    threads: 4
    shell:
        "bwa mem -t {threads} {params.index} {input.r1} {input.r2} > {output} 2> {log}"


rule extract_dehosted_reads:
    message:
        "Extracting unmapped reads from bam file: {wildcards.id}"
    conda:
        "../envs/read_processing.yaml"
    # container: "docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1"
    log:
        "results/logs/{id}/{id}_samtools_extract_non_human.log"
    benchmark:
        "results/benchmarks/samtools_extract/{id}.tsv"
    input:
        mapped_reads=rules.map_reads_to_human_genome.output.mapped_bam,
    output:
        unmapped_bam="results/{id}/{id}.bam",
    shell:
        """
        samtools view -b -f12 {input} > {output.unmapped_bam} 2> {log}
        """


rule generate_dehosted_fastq:
    message:
        "Converting BAM file into fastq files of forward and reverse reads: {wildcards.id}"
    conda:
        "../envs/read_processing.yaml"
    # container: "docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1"
    log:
        "results/logs/{id}/{id}_bamtofastq.log"
    benchmark:
        "results/benchmarks/bamtofastq/{id}.tsv"
    input:
        unmapped_bam_file=rules.extract_dehosted_reads.output.unmapped_bam,
    output:
        forward_read="results/{id}/{id}_1.fastq",
        reverse_read="results/{id}/{id}_2.fastq",
        non_human_stats="results/{id}/{id}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output.non_human_stats} 2> {log}
        samtools fastq -N -1 {output.forward_read} -2 {output.reverse_read} {input} > {log} 2>&1
        """
