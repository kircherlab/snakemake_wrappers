__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell


shell(
    "bedtools makewindows "
    "-g {snakemake.input[0]} "
    "-w {snakemake.params[0]} -s {snakemake.params[1]} | "
    "awk -v 'OFS=\\t' '{{print $0,0}}' | "
    "bgzip -c > {snakemake.output[0]}"
)
