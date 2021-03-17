# Author Max Schubach, 2021

"""
Concat tables with header using pandas
"""

import pandas as pd

import click

# options


@click.command()
@click.option('--input',
              'input_files',
              required=True,
              multiple=True,
              type=click.Path(exists=True, readable=True),
              help='Input TSV file with headers')
@click.option('--column',
              'columns',
              multiple=True,
              required=False,
              type=(str, str),
              help='Adding collumns with value. can be given for each file seperately or one for all')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output TSV file with headers.')
def cli(input_files, columns, output_file):
    if columns:
        if len(columns) != 1 and len(input_files) != len(columns):
            raise Exception("column option must be 1 or the length of the files")

    output = pd.DataFrame()

    for i, input_file in enumerate(input_files):
        df = pd.read_csv(input_file, sep="\t")
        if columns and len(columns) > 1:
            df[columns[i][0]] = columns[i][1]
        output = pd.concat([output, df])

    if columns and len(columns) == 1:
        output[columns[0][0]] = columns[0][1]
    output.to_csv(output_file, index=False, sep="\t")


if __name__ == '__main__':
    cli()
