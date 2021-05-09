from snakemake.utils import validate
import pandas as pd

wfbasedir = workflow.basedir
configfile: workflow.basedir + "/config.yaml"

##### complete species resource paths if following files defined with them #####
for config_resource in ['viral_reference_genome', 
                        'viral_sequencing_adapters',
                        'viral_sequencing_primers',
                        'viral_reference_gene_features']:
    config[config_resource] = config[config_resource].format(viral_species=config['viral_species'])

##### load sample sheets #####
sample_table = pd.read_csv(config["sample_table"], sep="\t").set_index("id", drop=False)
sample_table.index.names = ["id"]

IDS = sorted(sample_table['id'].drop_duplicates().values)

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = sample_table.loc[wildcards.id, ["r1", "r2"]].dropna()
    return {"r1": fastqs.r1, "r2": fastqs.r2}


rule all:
 input:
  mapped_bam = expand(["results/{id}/{id}.sam"], id=IDS),
  unmapped_bam = expand(["results/{id}/{id}.bam"], id=IDS),
  forward_reads = expand(["results/{id}/{id}_1.fastq"], id=IDS),
  reverse_reads = expand(["results/{id}/{id}_2.fastq"], id=IDS),
  contigs = expand(["results/{id}/{id}_iva"], id=IDS),
  initialization_directory = expand(["results/{id}/{id}_shiver_init_dir"], id=IDS),
  blast_hits = expand(["results/{id}/{id}.blast"], id=IDS),
  aligned_contigs_raw = expand(["results/{id}/{id}_raw_wRefs.fasta"], id=IDS),
  aligned_contigs_cut = expand(["results/{id}/{id}_cut_wRefs.fasta"], id=IDS),
  bam_file = expand(["results/{id}/{id}.bam"], id=IDS),
  ref_seqs = expand(["results/{id}/{id}_ref.fasta"], id=IDS),
  base_freqs = expand(["results/{id}/{id}_BaseFreqs.csv"], id=IDS),
  base_freqs_global_aln = expand(["results/{id}/{id}_BaseFreqs_ForGlobalAln.csv"], id=IDS),
  coords = expand(["results/{id}/{id}_coords.csv"], id=IDS),
  insert_size_dist = expand(["results/{id}/{id}_InsertSizeCounts.csv"], id=IDS),
  consensus_genome = expand(["results/{id}/{id}_remap_consensus_MinCov_10_30.fasta"], id=IDS),                       
  quast_results = expand(["results/{id}/{id}_quast_results"], id=IDS)

rule indexing:
 message: "Indexing the human reference genome"
 input:
  human_ref_genome = config['human_reference_genome']
 output:
    multiext(config['human_reference_genome'], ".amb", ".ann", ".bwt", ".pac", ".sa"),
 shell:
  "bwa index {input}"
  
rule map_to_human_genome:
 message: "Mapping reads to the human genome to remove human contaminants"
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
 input:
  mapped_reads = rules.map_to_human_genome.output.mapped_bam
 output:
  unmapped_bam = "results/{id}/{id}.bam"
 shell:
  "samtools view -b -f12 {input} > {output}"
  
rule bamtoFastq:
 message: "Converting BAM file into fastq files of forward and reverse reads"
 input:
  unmapped_bam_file = rules.extract_unmapped_reads.output.unmapped_bam
 output:
  forward_read = "results/{id}/{id}_1.fastq",
  reverse_read = "results/{id}/{id}_2.fastq"
 shell:
  "samtools fastq -N -1 {output[0]} -2 {output[1]} {input}"

rule assembly:
 message: "Assembly of forward and reverse reads"
 input:
  forward_read = rules.bamtoFastq.output.forward_read,
  reverse_read = rules.bamtoFastq.output.reverse_read
 output:
  contigs = directory("results/{id}/{id}_iva")
 threads: 8
 shell:
  "iva --reads_fwd {input.forward_read} -reads_rev {input.reverse_read} --threads {threads} {output}"



rule shiver_init:
 message: "Shiver initialization"
 input:
  Reference_alignment = config['viral_reference_genome'],
  Adapters = config['viral_sequencing_adapters'],
  Primers = config['viral_sequencing_primers']
 output:
  initialization_directory = directory("results/{id}/{id}_shiver_init_dir")
 shell:
   "shiver_init.sh {output} {wfbasedir}/config.sh {input[0]} {input[1]} {input[2]}"

rule align_contigs:
 message: "Aligning contigs"
 input:
  initialization_directory = rules.shiver_init.output.initialization_directory,
  contigs_file = rules.assembly.output.contigs
 output:
  blast_hits = "results/{id}/{id}.blast",
  aligned_contigs_raw = "results/{id}/{id}_raw_wRefs.fasta",
  aligned_contigs_cut = "results/{id}/{id}_cut_wRefs.fasta"
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
  ref_seqs = "results/{id}/{id}_ref.fasta",
  base_freqs = "results/{id}/{id}_BaseFreqs.csv",
  base_freqs_global_aln = "results/{id}/{id}_BaseFreqs_ForGlobalAln.csv",
  coords = "results/{id}/{id}_coords.csv",
  insert_size_dist = "results/{id}/{id}_InsertSizeCounts.csv",
  consensus_genome = "results/{id}/{id}_remap_consensus_MinCov_10_30.fasta"
 shell:
  "shiver_map_reads.sh {input} {wfbasedir}/config.sh {input[1]}/contigs.fasta 934 \
 {input[2]} {input[3]} {input[4]} {input[5]}"
  
 #934 in the shell of rule map should be changed to the sample ID

rule assembly_assessment:
  message: "Evaluate the quality of genome assembly"
  input:
   consensus_genome = "results/{id}/{id}_remap_consensus_MinCov_10_30.fasta",
   quast_ref_genome = config['viral_reference_genome'],
   gene_features = config['viral_reference_gene_features']
  output:
   quast_results = directory("results/{id}/{id}_quast_results")
  shell:
   "python quast.py -r {input[1]} -g {input[2]} -o {output.quast_results} {input[0]}"
