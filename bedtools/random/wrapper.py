__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
__email__ = "max.schubach@bihealth.de"
__license__ = "MIT license"

from snakemake.shell import shell

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

if "g" in snakemake.input.keys():
    input_g = snakemake.input["g"]
else:
    raise MissingInputException("g")
 
 
# Checking parameters, remove else when parameter is not necessary and add a default value

if "l" in snakemake.params.keys():
    param_l = snakemake.params["l"]
else:
    raise MissingParameterException("l")

if "n" in snakemake.params.keys():
    param_n = snakemake.params["n"]
else:
    raise MissingParameterException("n")

if "seed" in snakemake.params.keys():
    param_seed = "-seed %s" % str(snakemake.params["seed"])
else:
    param_seed = ""

# running the shell
shell(
    """
    bedtools random -g {input_g} \
    -l {param_l} -n {param_n} {param_seed} | \
    sort -k1,1 -k2,2n | \
    bgzip -c >  {snakemake.output} 
    """
)
