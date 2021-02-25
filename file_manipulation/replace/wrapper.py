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

if "columns" in snakemake.params.keys():
    param_columns = snakemake.params["columns"]
else:
    raise MissingParameterException("columns")

if "pat" in snakemake.params.keys():
    param_pat = snakemake.params["pat"]
else:
    raise MissingParameterException("pat")

if "replace" in snakemake.params.keys():
    param_replace = snakemake.params["replace"]
else:
    raise MissingParameterException("replace")

param_list=[]

for i in range(len(param_columns)):
    param_list.append('--replace %s "%s" %s' % (param_columns[i],param_pat[i],param_replace[i]))

params = " ".join(param_list)
# running the shell
shell(
    """
    python  {scriptFolder}/replace.py \
    {params} \
    --input {snakemake.input}  \
    --output {snakemake.output}
    """
)
