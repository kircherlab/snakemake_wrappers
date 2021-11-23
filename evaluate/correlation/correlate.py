import click
import pandas as pd


# options
@click.command()
@click.option('--fileA',
              'file_a',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='File A with values')
@click.option('--fileB',
              'file_b',
              required=False,
              type=click.Path(exists=True, readable=True),
              help='File B with the other values. if not set file A will be used')
@click.option('--values',
              'values',
              required=True,
              type=(str, str),
              help='Two headers to correlate starting with file A')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with Pearson and Spearman correlations.')
def cli(file_a, file_b, values, output_file):

    df = pd.read_csv(file_a, delimiter="\t")

    value_a = df[[values[0]]].iloc[:, 0].rename("A")

    if file_b:
        value_b = pd.read_csv(file_a, delimiter="\t")[
            [values[1]]].iloc[:, 0].rename("B")
    else:
        value_b = df[[values[1]]].iloc[:, 0].rename("B")

    df = pd.concat([value_a, value_b], axis=1)

    output = []
    methods = pd.Series(['pearson', 'spearman'], name="metric")
    for method in methods:
        output.append(df.corr(method).loc["A", "B"])

    pd.concat([methods, pd.Series(output, name='value')],
              axis=1).to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == '__main__':
    cli()
