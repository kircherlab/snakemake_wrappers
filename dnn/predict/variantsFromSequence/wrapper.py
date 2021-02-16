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
fileType = "--file-type TSV"

if "altMinusRef" in snakemake.params.keys():
    altMinusRef = "--altminusref" if snakemake.params["altMinusRef"] else "--refminusalt"
if "fileType" in snakemake.params.keys():
    fileType = "--file-type %s" % snakemake.params["fileType"]

if "variants" not in snakemake.input.keys():
    raise MissingInputException("variants")
if "model" not in snakemake.input.keys():
    raise MissingInputException("model")
if "weights" not in snakemake.input.keys():
    raise MissingInputException("weights")
if "reference" not in snakemake.input.keys():
    raise MissingInputException("reference")
if "reference_index" not in snakemake.input.keys():
    raise MissingInputException("reference_index")
if "genome" not in snakemake.input.keys():
    raise MissingInputException("genome")


shell(
    """
    python {scriptFolder}/predictVariantsFromSequence.py \
    --variants {snakemake.input.variants} \
    --model {snakemake.input.model} --weights {snakemake.input.weights} \
    --reference {snakemake.input.reference} --genome {snakemake.input.genome} \
    {altMinusRef} {fileType} \
    --output {snakemake.output}
    """
)