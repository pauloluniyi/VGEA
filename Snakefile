IDS, = glob_wildcards("{id}.bam")

rule all:
 input:
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
  min_cov = expand(["{id}_MinCov_X_Y.fasta"], id=IDS),
  min_cov_global_aln = expand(["{id}_MinCov_X_Y_ForGlobalAln.fasta"], id=IDS),
  coords = expand(["{id}_coords.csv"], id=IDS),
  insert_size_dist = expand(["{id}_InsertSizeCounts.csv"], id=IDS)

rule bamtoFastq:
 message: "Converting BAM file into fastq files of forward and reverse reads"
 input:
  expand(["{id}.bam"], id=IDS)
 output:
  forward_read = "{id}_1.fastq",
  reverse_read = "{id}_2.fastq"
 conda:
  "vgea.yml"
 shell:
  "samtools fastq -N -1 {output[0]} -2 {output[1]} {input}"

rule assembly:
 message: "Assembly of forward and reverse reads"
 input:
  forward_read = rules.bamtoFastq.output.forward_read,
  reverse_read = rules.bamtoFastq.output.reverse_read
 output:
  contigs = directory("{id}_iva")
 conda:
  "vgea.yml"
 shell:
  "iva -f {input[0]} -r {input[1]} {output}"

rule shiver_init:
 message: "Shiver initialization"
 input:
  Reference_alignment = "MyRefAlignment.fasta",
  Adapters = "MyAdapters.fasta",
  Primers = "MyPrimers.fasta"
 output:
  initialization_directory = directory("MyInitDir")
 conda:
  "vgea.yml"
 shell:
   "~/pipeline/scripts/shiver_init.sh {output} ~/pipeline/scripts/config.sh {input[0]} {input[1]} {input[2]}"

rule align_contigs:
 message: "Aligning contigs"
 input:
  initialization_directory = rules.shiver_init.output.initialization_directory,
  contigs_file = rules.assembly.output.contigs
 output:
  blast_hits = "{id}.blast",
  aligned_contigs_raw = "{id}_raw_wRefs.fasta",
  aligned_contigs_cut = "{id}_cut_wRefs.fasta"
 conda:
  "vgea.yml"
 shell:
  "~/pipeline/scripts/shiver_align_contigs.sh {input[0]} ~/pipeline/scripts/config.sh {input[1]}/contigs.fasta CV18"

#CV18 in the shell of rule align_contigs should be changed to the sample ID

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
  bam_file = "{id}.bam",
  ref_seqs = "{id}_ref.fasta",
  base_freqs = "{id}_BaseFreqs.csv",
  base_freqs_global_aln = "{id}_BaseFreqs_ForGlobalAln.csv",
  min_cov = "{id}_MinCov_X_Y.fasta",
  min_cov_global_aln = "{id}_MinCov_X_Y_ForGlobalAln.fasta",
  coords = "{id}_coords.csv",
  insert_size_dist = "{id}_InsertSizeCounts.csv"
 conda:
  "vgea.yml"
 shell:
  "~/pipeline/scripts/shiver_map_reads.sh {input[0]} ~/pipeline/scripts/config.sh {input[1]}/contigs.fasta CV18 \
 {input[2]} {input[3]} {input[4]} {input[5]}"

#CV18 in the shell of rule map should be changed to the sample ID
