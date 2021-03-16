import click
import pandas as pd


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input tsv file to rename rows and cols.')
@click.option('--column',
              'columns',
              required=False,
              multiple=True,
              type=(str, str),
              help='Column and replacement.')
@click.option('--row',
              'rows',
              required=False,
              multiple=True,
              type=(str, str),
              help='row and replacement.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file (tsv file with headers).')
def cli(input_file, columns, rows, output_file):

    columns_dict = {key: value for key, value in columns}
    rows_dict = {key: value for key, value in rows}
    pd.read_csv(input_file, dtype=object, delimiter="\t").rename(
        columns=columns_dict
    ).rename(
        index=rows_dict
    ).to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == '__main__':
    cli()
