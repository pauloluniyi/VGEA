rule assembly_assessment:
    message:
        "Evaluate the quality of genome assembly: {wildcards.id}"
    conda:
        "../envs/qc.yaml"
    # container: "docker://quay.io/biocontainers/quast:5.0.2--py36pl5262h30a8e3e_4"
    log:
        "results/logs/{id}/{id}_quast.log"
    benchmark:
        "results/benchmarks/quast/{id}.tsv"
    input:
        consensus_genome="results/{id}/{id}.shiver_assembly.fasta",
        quast_ref_genome=config["viral_reference_genome"],
        gene_features=config["viral_reference_gene_features"]
    output:
        quast_results="results/{id}/{id}.quast_results/report.tsv",
    params:
        quast_dir="results/{id}/{id}.quast_results"
    shell:
        "quast -r {input.quast_ref_genome} -g {input.gene_features} -o {params.quast_dir} {input.consensus_genome} > {log} 2>&1"

rule multiqc:
    message:
        "Combine results into multiqc report"
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/multiqc.log"
    benchmark:
        "results/logs/multiqc.tsv"
    input:
        quast = expand(["results/{id}/{id}.quast_results/report.tsv"], id=IDS),
        flagstat = expand(["results/{id}/{id}.flagstat"], id=IDS),
        fastp = expand(["results/{id}/{id}.fastp.json"], id=IDS)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        """
        multiqc results -o results/multiqc > {log} 2>&1 
        """


