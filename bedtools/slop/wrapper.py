__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell


if "b" in snakemake.params.keys():
    extend = "-b %d" % int(snakemake.params["b"])
elif "l" in snakemake.params.keys() and "r" in snakemake.params.keys():
    extend = "-s %d -l %d" % (int(snakemake.params["r"]),int(snakemake.params["r"]))
else:
    exit(1)


shell(
    "bedtools slop -g {snakemake.input.g} -s {extend} -i {snakemake.input.i} | sort -k1,1 -k2,2n | uniq | bgzip -c > {snakemake.output}"
)
