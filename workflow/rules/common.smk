from snakemake.utils import validate
import pandas as pd
from pathlib import Path


container: "docker://continuumio/miniconda3:4.4.10"


configfile: "config/config.yaml"


##### complete species resource paths if following files defined with them #####
for config_resource in config.keys():
    # excluding non-path configs
    if config_resource not in ["viral_species"]:
        config[config_resource] = config[config_resource].format(
            viral_species=config["viral_species"]
        )
        config[config_resource] = str(Path(config[config_resource]).resolve())

##### load sample sheets #####
sample_table = pd.read_csv(config["sample_table"], sep="\t").set_index("id", drop=False)
sample_table.index.names = ["id"]

IDS = sorted(sample_table["id"].drop_duplicates().values)


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = sample_table.loc[wildcards.id, ["r1", "r2"]].dropna()
    return {"r1": fastqs.r1, "r2": fastqs.r2}
