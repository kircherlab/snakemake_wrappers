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
              required=True,
              multiple=True,
              type=str,
              help='Column used for operation')
@click.option('--new-column-name',
              'new_column_names',
              required=True,
              multiple=True,
              type=str,
              help="Name of the new column with the operation. Must have the same length than --operation")
@click.option('--operation',
              'operations',
              required=True,
              multiple=True,
              type=click.Choice(['abs', 'max', 'min', 'mean', 'abs_max',
                                 'abs_mean', 'abs_min'], case_sensitive=False),
              help="Operation of the new column. Must have the same length than --new-column-name")
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with average columb (tsv file with headers).')
def cli(input_file, columns, new_column_names, operations, output_file):
    if len(new_column_names) != len(operations):
        raise CLIException("--new-column-name and --operation must have the same length to match operation to column.")

    df = pd.read_csv(input_file, dtype=object, delimiter="\t")

    switcher = {
        'abs': lambda df: df[list(columns)].astype(float).abs(axis=1),
        'max': lambda df: df[list(columns)].astype(float).max(axis=1),
        'min': lambda df: df[list(columns)].astype(float).min(axis=1),
        'mean': lambda df: df[list(columns)].astype(float).mean(axis=1),
        'abs_max': lambda df: df[list(columns)].astype(
            float).abs().max(axis=1),
        'abs_mean': lambda df:  df[list(columns)].astype(
            float).abs().mean(axis=1),
        'abs_min': lambda df:  df[list(columns)].astype(
            float).abs().min(axis=1),
    }

    for i, operation in enumerate(operations):
        if operation in ['abs'] and len(columns) != 1:
            raise CLIException(
                "Operation %s can only performed on one column" % operation)

        new_column_name = new_column_names[i]

        if operation not in switcher.keys():
            raise CLIException("Operation %s is not implemented" % operation)

        df[new_column_name] = switcher.get(operation)(df)

    df.to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == '__main__':
    cli()


class CLIException(Exception):
    """Exception raised for errors in the input.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        message (string): The message.
    """

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return(self.message)
