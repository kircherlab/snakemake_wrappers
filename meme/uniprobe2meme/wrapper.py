__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell

shell(
    "cat {snakemake.input[0]} | uniprobe2meme > {snakemake.output[0]}"
)
