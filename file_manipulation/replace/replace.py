# Author Max Schubach, 2021

"""
Replace entries in table
"""

import pandas as pd
from pandas.api.types import is_numeric_dtype

import click

# options


@click.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    multiple=False,
    type=click.Path(exists=True, readable=True),
    help="Input TSV file with headers",
)
@click.option(
    "--replace",
    "replacements",
    multiple=True,
    required=True,
    type=(str, str, str),
    help="Column, patter, repleacement",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output TSV file with headers.",
)
def cli(input_file, replacements, output_file):

    df = pd.read_csv(input_file, sep="\t")
    for replacement in replacements:
        if is_numeric_dtype(df[replacement[0]]):
            replacement[1] = float(replacement[1])
    for replacement in replacements:
        df[replacement[0]] = df[replacement[0]].replace(replacement[1], replacement[2])

    df.to_csv(output_file, index=False, sep="\t")


if __name__ == "__main__":
    cli()
