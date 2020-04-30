__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell


threads = 1 if snakemake.threads <= 1 else snakemake.threads

memory = "2GB" if len(snakemake.resources) == 0 else snakemake.resources[0]

shell(
    "negative_training_sampler -i {snakemake.input[0]} "
    "-r {snakemake.input[1]} "
    "-g {snakemake.input[2]} "
    "--cores {threads} "
    "--memory {memory} | "
    "awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5}}' | "
    "bgzip -c > {snakemake.output[0]}"
)
