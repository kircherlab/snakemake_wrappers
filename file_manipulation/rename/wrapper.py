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


# Checking parameters, remove else when parameter
# is not necessary and add a default value

if "columns" in snakemake.params.keys():
    param_columns = " ".join(
        ["--column %s %s" % (
            key, value
        ) for key, value in snakemake.params["columns"].items()]
    )
else:
    param_columns = ""

if "rows" in snakemake.params.keys():
    param_columns = " ".join(
        ["--row %s %s" % (
            key, value
        ) for key, value in snakemake.params["rows"].items()]
    )
else:
    param_rows = ""


# running the shell
shell(
    """
    python  {scriptFolder}/rename.py \
    {param_rows} \
    {param_columns} \
    --input {snakemake.input}  --output  {snakemake.output}
    """
)
