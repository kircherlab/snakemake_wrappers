__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell


threads = 1 if snakemake.threads <= 1 else snakemake.threads

memory = "2GB"
seed = ''


if len(snakemake.params) != 0:
    if "memory" in snakemake.params.keys():
        memory = snakemake.params["memory"]
    if "seed" in snakemake.params.keys():
        seed = "--seed %d " % int(snakemake.params["seed"])

log = ''
if snakemake.log != 0:
    log = "--log %s " % snakemake.log[0]

shell(
    "negative_training_sampler -i {snakemake.input[0]} "
    "-r {snakemake.input[1]} "
    "-g {snakemake.input[2]} "
    "--cores {threads} "
    "{seed}"
    "{log}"
    "--memory {memory} | "
    "bgzip -c > {snakemake.output[0]}"
)
