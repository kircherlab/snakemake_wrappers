__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os
from snakemake.shell import shell

scriptFolder = os.path.dirname(os.path.abspath(__file__))

class MissingParameterException(Exception):
    """Exception raised for errors in the parameter input.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        parameter (string): Parameter which is missing.
    """

    def __init__(self, parameter):
        self.parameter = parameter

    def __str__(self):
        return("Parameter %s is missing!" % (self.parameter))

class MissingInputException(Exception):
    """Exception raised for errors in the input.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        input (string): Input which is missing.
    """

    def __init__(self, input):
        self.input = input

    def __str__(self):
        return("Input %s is missing!" % (self.input))

altMinusRef = "--refminusalt"


if "altMinusRef" in snakemake.params.keys():
    altMinusRef = "--altminusref" if snakemake.params["altMinusRef"] else "--refminusalt"
if "fileType" in snakemake.params.keys():
    fileType = "--file-type %s" snakemake.params["fileType"]
else:
    raise MissingParameterException("fileType")

if "variants" in snakemake.input.keys():
    raise MissingInputException("variants")
if "model" in snakemake.input.keys():
    raise MissingInputException("model")
if "weights" in snakemake.input.keys():
    raise MissingInputException("weights")
if "reference" in snakemake.input.keys():
    raise MissingInputException("reference")
if "genome" in snakemake.input.keys():
    raise MissingInputException("genome")


shell(
    """
    python {scriptFolder}/predictVariantsFromSequence.py \
    --variants {snakemake.variants} \
    --model {snakemake.model} --weights {snakemake.weights} \
    --reference {snakemake.reference} --genome {snakemake.genome} \
    {altMinusRef} {fileType} \
    --output {snakemake.output}"
    """
)