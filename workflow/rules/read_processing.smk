rule indexing:
 message: "Indexing the human reference genome"
 conda: "../envs/read_processing.yaml"
 #container: "docker://quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
 input:
  human_ref_genome = config['human_reference_genome']
 output:
    multiext(config['human_reference_genome'], ".amb", ".ann", ".bwt", ".pac", ".sa"),
 shell:
  "bwa index {input}"
  
rule map_to_human_genome:
 message: "Mapping reads to the human genome to remove human contaminants"
 conda: "../envs/read_processing.yaml"
 #container: "docker://quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
 input:
  unpack(get_fastq),
  idx=rules.indexing.output,
 output:
  mapped_bam = "results/{id}/{id}.sam"
 params:
  index=lambda w, input: os.path.splitext(input.idx[0])[0]
 threads: 4
 shell:
  "bwa mem -t {threads} {params.index} {input.r1} {input.r2} > {output}"

rule extract_unmapped_reads:
 message: "Extracting unmapped reads from bam file"
 conda: "../envs/read_processing.yaml"
 #container: "docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1"
 input:
  mapped_reads = rules.map_to_human_genome.output.mapped_bam
 output:
  unmapped_bam = "results/{id}/{id}.bam"
 shell:
  "samtools view -b -f12 {input} > {output}"
  
rule bamtoFastq:
 message: "Converting BAM file into fastq files of forward and reverse reads"
 conda: "../envs/read_processing.yaml"
 #container: "docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1"
 input:
  unmapped_bam_file = rules.extract_unmapped_reads.output.unmapped_bam
 output:
  forward_read = "results/{id}/{id}_1.fastq",
  reverse_read = "results/{id}/{id}_2.fastq"
 shell:
  "samtools fastq -N -1 {output[0]} -2 {output[1]} {input}"

