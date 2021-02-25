# Author Max Schubach, 2021

"""
Melt table with header using pandas
"""

import sys
import pandas as pd

import click

# options


@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input TSV file with headers')
@click.option('--id',
              'id_vars',
              multiple=True,
              required=True,
              type=str,
              help='Column to use as identifier variables.')
@click.option('--value',
              'value_vars',
              multiple=True,
              required=False,
              type=str,
              help='Column to unpivot. If not specified, uses all columns that are not set as id_vars.')
@click.option('--value-name',
              'value_name',
              required=False,
              default='value',
              type=str,
              help='Name to use for the ‘value’ column.')
@click.option('--var-name',
              'var_name',
              required=False,
              default='variable',
              type=str,
              help='Name to use for the ‘variable’ column.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output TSV file with headers.')

def cli(input_file, id_vars, value_vars,value_name,var_name, output_file):
    
    df = pd.read_csv(input_file, sep="\t")
    
    if value_vars:
        output = pd.melt(df, id_vars=list(id_vars),value_vars=list(value_vars),var_name=var_name, value_name=value_name)
    else:
        output = pd.melt(df, id_vars=list(id_vars),var_name=var_name, value_name=value_name)
    
    output.to_csv(output_file, index=False,sep="\t")
    
if __name__ == '__main__':
    cli()
