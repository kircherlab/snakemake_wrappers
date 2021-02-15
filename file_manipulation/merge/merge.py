import click
import pandas as pd


# options
@click.command()
@click.option('--left',
              'file_left',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Left file for join')
@click.option('--right',
              'file_right',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Right fiel for join')
@click.option('--left-on',
              'left_on',
              required=True,
              multiple=True,
              type=str,
              help='Column used for file left (or for both if right-on is not set')
@click.option('--right-on',
              'right_on',
              required=False,
              multiple=True,
              type=str,
              help='Column used for right file')
@click.option('--how',
              'how',
              required=True,
              multiple=False,
              type=click.Choice(
                  ['left', 'right', 'outer', 'inner', 'cross'], case_sensitive=False),
              help='How the join should be used')
@click.option('--suffixes',
              'suffixes',
              required=False,
              multiple=False,
              type=(str, str),
              default=('_X', '_Y'),
              help='A length-2 sequence where each element is optionally a string indicating the suffix to add to overlapping column names in left and right respectively')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with join.')
def cli(file_left, file_right, left_on, right_on, how, suffixes, output_file):

    df_left = pd.read_csv(file_left, dtype=object, delimiter="\t")

    df_right = pd.read_csv(file_right, dtype=object, delimiter="\t")

    if not right_on:
        right_on = left_on

    if (len(left_on) != len(right_on)):
        raise Exception("left-on and right-on has different length")

    df = df_left.merge(df_right, how=how, left_on=left_on,
                       right_on=right_on, suffixes=suffixes)

    df.to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == '__main__':
    cli()
