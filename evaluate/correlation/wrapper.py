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

if "a" in snakemake.input.keys():
    input_a = snakemake.input["a"]
else:
    raise MissingInputException("a")
 
if "b" in snakemake.input.keys():
    input_b = "--fileB %s" % snakemake.input["b"]
else:
    input_b = ""
 
 

# Checking parameters, remove else when parameter is not necessary and add a default value

if "value_a" in snakemake.params.keys():
    param_valueA = snakemake.params["value_a"]
else:
    raise MissingParameterException("value_a")
 
if "value_b" in snakemake.params.keys():
    param_valueB = snakemake.params["value_b"]
else:
    raise MissingParameterException("value_b")
 


# running the shell
shell(
    """
    python  {scriptFolder}/correlate.py \
    --values {param_valueA} {param_valueB} \
    --fileA {input_a}  {input_b}  \
    --output {snakemake.output} 
    """
)
