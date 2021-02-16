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
if "regions" in snakemake.input.keys():
    input_regions = snakemake.input["regions"]
else:
    raise MissingInputException("regions")
 
if "model" in snakemake.input.keys():
    input_model = snakemake.input["model"]
else:
    raise MissingInputException("model")
 
if "weights" in snakemake.input.keys():
    input_weights = snakemake.input["weights"]
else:
    raise MissingInputException("weights")
 
if "reference" in snakemake.input.keys():
    input_reference = snakemake.input["reference"]
else:
    raise MissingInputException("reference")

if "reference_index" in snakemake.input.keys():
    input_reference = snakemake.input["reference_index"]
else:
    raise MissingInputException("reference_index")

if "genome" in snakemake.input.keys():
    input_genome = snakemake.input["genome"]
else:
    raise MissingInputException("genome")
 
 

# Checking parameters, remove else when parameter is not necessary and add a default value
if "leftEdge" in snakemake.params.keys():
    param_leftEdge = snakemake.params["leftEdge"]
else:
    param_leftEdge = 0

if "rightEdge" in snakemake.params.keys():
    param_rightEdge = snakemake.params["rightEdge"]
else:
    param_rightEdge = 0

if "refMinusAlt" in snakemake.params.keys():
    altMinusRef = "--refminusalt"
else:
    altMinusRef = "--altminusref"


# running the shell
shell(
    """
    python  {scriptFolder}/predictVariantsWithInsilicoSaturationMutagenesis.py \
    --edges {param_leftEdge} {param_rightEdge} {altMinusRef} \
    --regions {input_regions} --model {input_model} --weights {input_weights} \
    --reference {input_reference} --genome {input_genome} \
    --output {snakemake.output} 
    """
)
