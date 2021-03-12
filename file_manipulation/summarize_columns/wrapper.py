__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
__email__ = "max.schubach@bihealth.de"
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
    param_columns = " ".join(
        ["--column %s" % i for i in snakemake.params["columns"]])
else:
    raise MissingParameterException("columns")

if "new_columns" in snakemake.params.keys():
    param_new_columns = " ".join(
        ["--new-column-name %s" % i for i in snakemake.params["new_columns"]])
else:
    raise MissingParameterException("new_columns")

if "operations" in snakemake.params.keys():
    param_operations = " ".join(
        ["--operation %s" % i for i in snakemake.params["operations"]])
else:
    raise MissingParameterException("operation")


# running the shell
shell(
    """
    python  {scriptFolder}/summarize_columns.py \
    {param_operations} \
    {param_columns}  \
    {param_new_columns}  \
    --input {snakemake.input}  --output  {snakemake.output} 
    """
)
