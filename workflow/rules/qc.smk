rule assembly_assessment:
    message:
        "Evaluate the quality of genome assembly"
    conda:
        "../envs/qc.yaml"
    # container: "docker://quay.io/biocontainers/quast:5.0.2--py36pl5262h30a8e3e_4"
    log:
        "results/logs/{id}/{id}_quast.log"
    input:
        consensus_genome="results/{id}/shiver/{id}_remap_consensus_MinCov_10_30.fasta",
        quast_ref_genome=config["viral_reference_genome"],
        gene_features=config["viral_reference_gene_features"],
    output:
        quast_results=directory("results/{id}/{id}_quast_results"),
    shell:
        "python quast.py -r {input.quast_ref_genome} -g {input.gene_features} -o {output.quast_results} {input.consensus_genome}"
