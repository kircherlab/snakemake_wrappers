import click
import pandas as pd


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input tsv fiel with multple columns for average.')
@click.option('--column',
              'columns',
              required=False,
              multiple=True,
              type=str,
              help='Column used for generating the average')
@click.option('--new-column-name',
              'new_column_name',
              required=True,
              multiple=False,
              type=str,
              help='Name of the new column with averages of columns')
@click.option('--abs/--no-abs',
              'absolute',
              default=False,
              required=False,
              help='Use abs before average')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with average columb (tsv file with headers).')
def cli(input_file, columns, new_column_name, absolute, output_file):

    df = pd.read_csv(input_file, dtype=object, delimiter="\t")

    if absolute:
        df[new_column_name] = df[list(columns)].astype(float).abs().mean(axis=1)
    else:
        df[new_column_name] = df[list(columns)].astype(float).mean(axis=1)

    df.to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == '__main__':
    cli()
