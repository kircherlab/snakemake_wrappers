__author__ = "Max Schubach"
__copyright__ = "Copyright 2022"
__email__ = "max.schubach@bih-charite.de"
__license__ = "MIT license"

import os
from snakemake.shell import shell

scriptFolder = os.path.dirname(os.path.abspath(__file__))


class MissingParameterException(Exception):
    """Exception raised for errors in the input.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        parameter (string): Parameter which is missing.
    """

    def __init__(self, parameter):
        self.parameter = parameter

    def __str__(self):
        return("Parameter %s is missing!" % (self.parameter))


# Checking parameters, remove else when parameter is not necessary and add a default value

if "xname" in snakemake.params.keys():
    param_xname = "--xname %s" % snakemake.params["xname"]
else:
    param_xname = ""

if "yname" in snakemake.params.keys():
    param_yname = "--yname %s" % snakemake.params["yname"]
else:
    param_yname = ""



# running the shell
shell(
    """
    Rscript {scriptFolder}/pre_re_f1_f2.R \
    {param_xname} {param_yname} \
    --input {snakemake.input} \
    --output {snakemake.output}
    """
)
