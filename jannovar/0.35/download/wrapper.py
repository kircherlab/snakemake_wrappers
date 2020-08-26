__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell


shell(
    "jannovar annotate-vcf -d {snakemake.input.dd} -i {snakemake.input.i} -o {snakemake.output.o}"
)
