__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
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


# Checking parameters, remove else when parameter
# is not necessary and add a default value

if "column" in snakemake.params.keys():
    param_column = "--column %s" % str(snakemake.params["column"])
else:
    raise MissingParameterException("column")

if "header" in snakemake.params.keys():
    param_header = "--header" if snakemake.params["header"] else "--no-header"
else:
    param_header = "--no-header"

if "chunksize" in snakemake.params.keys():
    param_chunksize = "--chunksize %d" % snakemake.params["chunksize"]
else:
    param_chunksize = "--chunksize 10000"


# running the shell
shell(
    """
    python  {scriptFolder}/nucleotideCountPerPosition.py \
    {param_header} \
    {param_column} \
    {param_chunksize} \
    --input {snakemake.input}  --output  {snakemake.output}
    """
)
