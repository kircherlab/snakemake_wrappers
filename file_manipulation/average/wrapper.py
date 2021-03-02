__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
__email__ = "max.schubach@bihealth.de"
__license__ = "MIT license"

import os
from snakemake.shell import shell

scriptFolder = os.path.dirname(os.path.abspath(__file__))

class MissingInputException(Exception):
    """Exception raised for errors in the input.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        inp (string): Input which is missing.
    """

    def __init__(self, inp):
        self.inp = inp

    def __str__(self):
        return("Input %s is missing!" % (self.inp))

class MissingOutputException(Exception):
    """Exception raised for errors in the input.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        output (string): Output which is missing.
    """

    def __init__(self, output):
        self.output = output

    def __str__(self):
        return("Output %s is missing!" % (self.output))

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


# Checking inputs


# Checking outputs

# Checking parameters, remove else when parameter is not necessary and add a default value

if "columns" in snakemake.params.keys():
    param_columns = " ".join(["--column %s" % i for i in snakemake.params["columns"]])
else:
    raise MissingParameterException("columns")
 
if "new_column" in snakemake.params.keys():
    param_new_column = snakemake.params["new_column"]
else:
    raise MissingParameterException("new_column")

if "absolute" in snakemake.params.keys():
    param_absolute = "--abs" if snakemake.params["absolute"] else "--no-abs"
else:
    param_absolute = "--no-abs"


# running the shell
shell(
    """
    python  {scriptFolder}/average.py \
    {param_absolute} \
    {param_columns}  \
    --new-column-name {param_new_column}  \
    --input {snakemake.input}  --output  {snakemake.output} 
    """
)
