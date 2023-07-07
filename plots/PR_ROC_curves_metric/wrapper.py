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

if "labelcolumns" in snakemake.params.keys():
    param_labelcolumns = "--labelcolumns %s" % snakemake.params["labelcolumns"]
else:
    param_labelcolumns = ""

if "type" in snakemake.params.keys():
    param_type = snakemake.params["type"]
else:
    raise MissingParameterException("type")

if "names" in snakemake.params.keys():
    param_name = "--name '%s'" % ",".join(snakemake.params["names"])
else:
    raise MissingParameterException("names")

inputs = ",".join(snakemake.input)


if len(snakemake.log) > 0:
    log = "&>  %s" % snakemake.log[0]
else:
    log = ""

# running the shell
shell(
    """
    Rscript {scriptFolder}/pr_roc_curve.R \
    {param_xname} {param_yname} {param_name} {param_labelcolumns}\
    --type {param_type} \
    --input {inputs} \
    --output {snakemake.output} {log}
    """
)
