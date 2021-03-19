import click
import pandas as pd


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input tsv file containing headers with multple columns for absolute.')
@click.option('--column',
              'columns',
              required=True,
              multiple=True,
              type=str,
              help='Column used for abs')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with values replaced by absolute (tsv file with headers).')
def cli(input_file, columns, output_file):

    df = pd.read_csv(input_file, delimiter="\t")

    cols = [col for col in df.columns if col not in list(columns)]

    df = pd.concat([df[cols], df[list(columns)].abs()], axis=1)[df.columns]

    df.to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == '__main__':
    cli()
