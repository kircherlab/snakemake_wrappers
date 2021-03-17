__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
__email__ = "max.schubach@bihealth.de"
__license__ = "MIT license"

from snakemake.shell import shell

input = " ".join(["<( zcat %s )" % i for i in snakemake.input])


# running the shell
shell(
    """
    paste {input} | gzip -c >  {snakemake.output}
    """
)
