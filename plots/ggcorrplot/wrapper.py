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

if "columns" in snakemake.params.keys():
    param_columns = "--column %s" % ",".join(snakemake.params["columns"])
else:
    param_columns = ""

if "bind" in snakemake.params.keys():
    param_bind = True
else:
    param_bind = False

if "arrange" in snakemake.params.keys():
    param_arrange = "--arrange %s" % snakemake.params["arrange"]
else:
    param_arrange = ""

if "method" in snakemake.params.keys():
    param_method = "--method %s " % snakemake.params["method"]
else:
    param_method = "spearman"


input = ",".join(snakemake.input)

# running the shell
shell(
    """
    Rscript {scriptFolder}/ggcorrplot.R \
    {param_columns} {param_bind} {param_arrange} {param_method} \
    --input {input} \
    --output {snakemake.output} 
    """
)
