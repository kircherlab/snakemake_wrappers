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

if "group_vars" in snakemake.params.keys():
    param_group_vars = " ".join(["--group %s" %i for i in snakemake.params["group_vars"]])
else:
    param_group_vars= ""

# running the shell
shell(
    """
    python  {scriptFolder}/scale.py \
    --input {snakemake.input} \
    {param_id_vars}  {param_group_vars}  \
    --output {snakemake.output} 
    """
)
