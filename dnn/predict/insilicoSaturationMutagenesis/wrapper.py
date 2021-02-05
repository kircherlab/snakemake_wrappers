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
if "regions" in snakemake.inputs.keys():
    input_regions = snakemake.inputs["regions"]
else:
    raise MissingInputException("regions")
 
if "model" in snakemake.inputs.keys():
    input_model = snakemake.inputs["model"]
else:
    raise MissingInputException("model")
 
if "weights" in snakemake.inputs.keys():
    input_weights = snakemake.inputs["weights"]
else:
    raise MissingInputException("weights")
 
if "reference" in snakemake.inputs.keys():
    input_reference = snakemake.inputs["reference"]
else:
    raise MissingInputException("reference")

if "genome" in snakemake.inputs.keys():
    input_genome = snakemake.inputs["genome"]
else:
    raise MissingInputException("genome")
 

 # Checking outputs
if "output" in snakemake.outputs.keys():
    output_output = snakemake.outputs["output"]
else:
    raise MissingOuputException("output")
 

# Checking parameters, remove else when parameter is not necessary and add a default value
if "notPredictedEdge" in snakemake.params.keys():
    param_notPredictedEdge = snakemake.params["notPredictedEdge"]
else:
    raise MissingParameterException("notPredictedEdge")

if "refMinusAlt" in snakemake.params.keys():
    altMinusRef = "--refminusalt"
else:
    altMinusRef = "--altminusref"


# running the shell
shell(
    """
    python  {scriptFolder}/predictVariantsWithInsilicoSaturationMutagenesis.py \
    --edge {param_notPredictedEdge}  {altMinusRef} \
    --regions {input_regions} --model {input_model} --weights {input_weights} \
    --reference {input_reference} --genome {input_genome} \
    --output {output_output} 
    """
)
