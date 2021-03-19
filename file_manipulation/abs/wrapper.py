__author__ = "Max Schubach"
__copyright__ = "Copyright 2021"
__email__ = "max.schubach@bihealth.de"
__license__ = "MIT license"

import os
from snakemake.shell import shell

scriptFolder = os.path.dirname(os.path.abspath(__file__))

# running the shell
shell(
    """
    python  {scriptFolder}/abs.py\
     --input {snakemake.input}  \
     --output {snakemake.output}
    """
)
