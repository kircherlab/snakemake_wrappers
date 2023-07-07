# Author Max Schubach, 2023

"""
Preprocess data by scaling each feature to a given range.
"""

from sklearn import preprocessing
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
@click.option('--group',
              'group_vars',
              multiple=True,
              required=False,
              type=str,
              help='Grouping on value.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output TSV file with headers.')

def cli(input_file, id_vars, group_vars, output_file):
    
    df = pd.read_csv(input_file, sep="\t")

    for id_var in id_vars:
        if group_vars:
            df[id_var] = df.groupby(list(group_vars))[id_var].transform(lambda x: preprocessing.scale(x))
        else:
            df[id_var] = preprocessing.scale(df[id_var])

    df.to_csv(output_file, index=False, sep="\t")
    
if __name__ == '__main__':
    cli()
