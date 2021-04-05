IDS, = glob_wildcards("{id}_R1.fastq"),
IDS, = glob_wildcards("{id}_R2.fastq")

wfbasedir = workflow.basedir
configfile: workflow.basedir + "/config.yaml"

rule all:
 input:
  mapped_bam = expand(["{id}.sam", id=IDS),
  unmapped_bam = expand(["{id}.bam", id=IDS),
  forward_reads = expand(["{id}_1.fastq"], id=IDS),
  reverse_reads = expand(["{id}_2.fastq"], id=IDS),
  contigs = expand(["{id}_iva"], id=IDS),
  initialization_directory = expand(["MyInitDir"]),
  blast_hits = expand(["{id}.blast"], id=IDS),
  aligned_contigs_raw = expand(["{id}_raw_wRefs.fasta"], id=IDS),
  aligned_contigs_cut = expand(["{id}_cut_wRefs.fasta"], id=IDS),
  bam_file = expand(["{id}.bam"], id=IDS),
  ref_seqs = expand(["{id}_ref.fasta"], id=IDS),
  base_freqs = expand(["{id}_BaseFreqs.csv"], id=IDS),
  base_freqs_global_aln = expand(["{id}_BaseFreqs_ForGlobalAln.csv"], id=IDS),
  coords = expand(["{id}_coords.csv"], id=IDS),
  insert_size_dist = expand(["{id}_InsertSizeCounts.csv"], id=IDS),
  quast_results = directory("quast_results")

rule indexing:
 message: "Indexing the human reference genome"
 input:
  human_ref_genome = config['hg19.fasta']
 shell:
  "bwa index {input}"
  
rule map_to_human_genome:
 message: "Mapping reads to the human genome to remove human contaminants"
 input:
  reads_1 = "{id}_R1.fastq",
  reads_2 = "{id}_R2.fastq",
  human_ref_genome = config['hg19.fasta']
 output:
  mapped_bam = "{id}.sam"
 shell:
  "bwa mem {input[2]} {input[0]} {input[1]} > {output}

rule extract_unmapped_reads:
 message: "Extracting unmapped reads from bam file"
 input:
  mapped_reads = rules.map_to_human_genome.output.mapped_bam
 output:
  unmapped_bam = "{id}.bam"
 shell:
  "samtools view -b -f12 {input} > {output}
  
rule bamtoFastq:
 message: "Converting BAM file into fastq files of forward and reverse reads"
 input:
  unmapped_bam_file = rules.extract_unmapped_reads.output.unmapped_bam
 output:
  forward_read = "{id}_1.fastq",
  reverse_read = "{id}_2.fastq"
 shell:
  "samtools fastq -N -1 {output[0]} -2 {output[1]} {input}"

rule assembly:
 message: "Assembly of forward and reverse reads"
 input:
  forward_read = rules.bamtoFastq.output.forward_read,
  reverse_read = rules.bamtoFastq.output.reverse_read
 output:
  contigs = directory("{id}_iva")
 shell:
  "iva -f {input[0]} -r {input[1]} {output}"

rule shiver_init:
 message: "Shiver initialization"
 input:
  Reference_alignment = config['Reference_alignment'],
  Adapters = config['Adapters'],
  Primers = config['Primers']
 output:
  initialization_directory = directory("MyInitDir")
 shell:
   "shiver_init.sh {output} {wfbasedir}/config.sh {input[0]} {input[1]} {input[2]}"

rule align_contigs:
 message: "Aligning contigs"
 input:
  initialization_directory = rules.shiver_init.output.initialization_directory,
  contigs_file = rules.assembly.output.contigs
 output:
  blast_hits = "{id}.blast",
  aligned_contigs_raw = "{id}_raw_wRefs.fasta",
  aligned_contigs_cut = "{id}_cut_wRefs.fasta"
 shell:
  "shiver_align_contigs.sh {input[0]} {wfbasedir}/config.sh {input[1]}/contigs.fasta 934"

#934 in the shell of rule align_contigs should be changed to the sample ID

rule map:
 message: "Mapping paired-end reads to reference genome"
 input:
  initialization_directory = rules.shiver_init.output.initialization_directory,
  contigs_file = rules.assembly.output.contigs,
  blast_hits = rules.align_contigs.output.blast_hits,
  aligned_contigs_cut = rules.align_contigs.output.aligned_contigs_cut,
  forward_read = rules.bamtoFastq.output.forward_read,
  reverse_read = rules.bamtoFastq.output.reverse_read
 output:
  ref_seqs = "{id}_ref.fasta",
  base_freqs = "{id}_BaseFreqs.csv",
  base_freqs_global_aln = "{id}_BaseFreqs_ForGlobalAln.csv",
  coords = "{id}_coords.csv",
  insert_size_dist = "{id}_InsertSizeCounts.csv"
 shell:
  "shiver_map_reads.sh {input[0]} {wfbasedir}/config.sh {input[1]}/contigs.fasta 934 \
 {input[2]} {input[3]} {input[4]} {input[5]}"
  
 #934 in the shell of rule map should be changed to the sample ID

 rule assembly_assessment/evaluation
  message: "Evaluate the quality of genome assembly"
  input:
   consensus_genome = "{id}_remap_consensus_MinCov_10_30.fasta",
   quast_ref_genome = config['quast_refseq.fasta'],
   gene_features = config['quast_genefeatures.txt']
  output:
   quast_results = directory("quast_results")
  shell:
   "python quast.py -r {input[1]} -g {input[2]} {input[0]}
   
