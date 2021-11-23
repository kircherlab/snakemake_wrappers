__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
__email__ = "max.schubach@bih-charite.de"
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

if "input" in snakemake.inputs.keys():
    input_input = snakemake.inputs["predictions"]
else:
    raise MissingInputException("predictions")


# Checking outputs

if "output" in snakemake.outputs.keys():
    output_output = snakemake.outputs["output"]
else:
    raise MissingOuputException("output")


# Checking parameters, remove else when parameter is not necessary and add a default value

if "score_column" in snakemake.params.keys():
    param_score_column = snakemake.params["score_column"]
else:
    raise MissingParameterException("score_column")


if "label_column" in snakemake.params.keys():
    param_label_column = snakemake.params["label_column"]
else:
    raise MissingParameterException("label_column")


if "positive_label" in snakemake.params.keys():
    param_positive_label = snakemake.params["positive_label"]
else:
    param_positive_label = "1"

# running the shell
shell(
    """
    python  {scriptFolder}/metric_plots.py \
    --score-column {param_score_column}  --label-column {param_label_column} \
    --positive-label {param_positive_label} \
    --input {input_input}  --output {output_output}
    """
)
