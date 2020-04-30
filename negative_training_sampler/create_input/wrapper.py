__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell


shell(
    "(zcat {snakemake.input[0]} | "
    "awk -v 'OFS=\\t' '{{print $1,$2,$3,1}}'; "
    "bedtools intersect -a {snakemake.input[1]} -b {snakemake.input[0]} -sorted -v) | "
    "sort -k1,1 -k2,2n | awk -v 'OFS=\\t' '{{print $1,$2,$3,$4}}' | "
    "bgzip -c > {snakemake.output}"
)
