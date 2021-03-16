# Author Max Schubach, 2021

"""
Extract columns of table
"""

import pandas as pd
import click

# options


@click.command()
@click.option('--input',
              'input_file',
              required=True,
              multiple=False,
              type=click.Path(exists=True, readable=True),
              help='Input TSV file with headers')
@click.option('--column',
              'columns',
              multiple=True,
              required=True,
              type=str,
              help='Column to extract')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output TSV file with headers.')
def cli(input_file, columns, output_file):
    
    pd.read_csv(input_file, sep="\t")[list(columns)].to_csv(output_file, index=False, sep="\t")


if __name__ == '__main__':
    cli()
