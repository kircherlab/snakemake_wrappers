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

if "id_vars" in snakemake.params.keys():
    param_id_vars = " ".join(["--id %s" %i for i in snakemake.params["id_vars"]])
else:
    raise MissingParameterException("id_vars")

if "value_vars" in snakemake.params.keys():
    param_value_vars = " ".join(["--value %s" %i for i in snakemake.params["value_vars"]])
else:
    param_value_vars= ""

if "value_name" in snakemake.params.keys():
    param_value_name = "--value-name %s" % snakemake.params["value_name"]
else:
    param_value_name= ""

if "var_name" in snakemake.params.keys():
    param_var_name = "--var-name %s" % snakemake.params["var_name"]
else:
    param_var_name= ""

# running the shell
shell(
    """
    python  {scriptFolder}/melt.py \
    --input {snakemake.input} \
    {param_id_vars}  {param_value_vars}  \
    {param_value_name} {param_var_name} \
    --output {snakemake.output} 
    """
)
