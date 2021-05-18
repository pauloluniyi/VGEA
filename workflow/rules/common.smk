from snakemake.utils import validate
import pandas as pd
from pathlib import Path

containerized: "docker://finlaymaguire/vgea:latest"

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##### complete species resource paths if following files defined with them #####
for config_resource in config.keys():
    # excluding non-path configs
    if config_resource not in ["viral_species"]:
        config[config_resource] = config[config_resource].format(
            viral_species=config["viral_species"]
        )
        config[config_resource] = str(Path(config[config_resource]).resolve())

##### load sample sheets #####
samples = pd.read_csv(config["sample_table"], sep="\t").set_index("id", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
IDS = samples.index.values

def get_fastq(wildcards):
    """Get fastq files of a given sample"""
    fastqs = samples.loc[wildcards.id, ["r1", "r2"]].dropna()
    return {"r1": fastqs.r1, "r2": fastqs.r2}
