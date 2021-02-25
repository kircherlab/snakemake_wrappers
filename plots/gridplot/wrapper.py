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


# Checking parameters, remove else when parameter is not necessary and add a default value

if "x" in snakemake.params.keys():
    param_x = snakemake.params["x"]
else:
    raise MissingParameterException("x")
 
if "y" in snakemake.params.keys():
    param_y = snakemake.params["y"]
else:
    raise MissingParameterException("y")
 
if "fill" in snakemake.params.keys():
    param_fill = "--fill %s" % snakemake.params["fill"]
else:
    param_fill = ""
 
if "group" in snakemake.params.keys():
    param_group = "--group %s " % snakemake.params["group"]
else:
    param_group = ""

if "colour" in snakemake.params.keys():
    param_colour = "--colour %s " % snakemake.params["colour"]
else:
    param_colour = ""
 
if "cols" in snakemake.params.keys():
    param_cols = snakemake.params["cols"]
else:
    raise MissingParameterException("cols")
 
if "rows" in snakemake.params.keys():
    param_rows = snakemake.params["rows"]
else:
    raise MissingParameterException("rows")
 
if "plot" in snakemake.params.keys():
    param_plot = snakemake.params["plot"]
else:
    raise MissingParameterException("plot")
 


# running the shell
shell(
    """
    Rscript {scriptFolder}/plotGrid.R \
    -x {param_x} -y {param_y} {param_fill} {param_group} {param_colour} --cols {param_cols} --rows {param_rows} --plot {param_plot}  \
    --input {snakemake.input} \
    --output {snakemake.output} 
    """
)
