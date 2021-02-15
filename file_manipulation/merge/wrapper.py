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

if "left" in snakemake.input.keys():
    input_left = snakemake.input["left"]
else:
    raise MissingInputException("left")

if "right" in snakemake.input.keys():
    input_right = snakemake.input["right"]
else:
    raise MissingInputException("right")
# Checking parameters, remove else when parameter is not necessary and add a default value

if "how" in snakemake.params.keys():
    param_how = snakemake.params["how"]
else:
    raise MissingParameterException("how")

if "right_on" in snakemake.params.keys():
    param_right_on = " ".join(["--right-on %s" % i for i in snakemake.params["right_on"].split(" ")])
else:
    param_right_on = ""

if "left_on" in snakemake.params.keys():
    param_left_on = " ".join(["--left-on %s" % i for i in snakemake.params["left_on"].split(" ")])
else:
    raise MissingParameterException("left_on")

if "suffixes" in snakemake.params.keys():
    param_suffixes = "--suffixes %s" % snakemake.params["suffixes"]
else:
    param_suffixes = ""


# running the shell
shell(
    """
    python  {scriptFolder}/merge.py \
    --how {param_how} \
    {param_suffixes} \
    {param_left_on} {param_right_on} \
    --left {input_left} --right {input_right} --output {snakemake.output} 
    """
)
