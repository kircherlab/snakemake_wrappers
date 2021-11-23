import click
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
import matplotlib.pyplot as plt


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='TSV file with scores and label')
@click.option('--score-column',
              'score_column',
              required=True,
              type=str,
              help='column name of the score')
@click.option('--label-column',
              'label_column',
              required=True,
              type=str,
              help='column name of the real label')
@click.option('--positive-label',
              'positive_label',
              required=True,
              default=1.0,
              type=float,
              help='column name of the positive label')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Final plot.')
def cli(input_file, score_column, label_column, positive_label, output_file):

    df = pd.read_csv(input_file, delimiter="\t")

    scores = df[[score_column]].iloc[:, 0]
    labels = df[[label_column]].iloc[:, 0]

    fpr, tpr, _ = roc_curve(labels, scores, pos_label=positive_label)
    roc_display = RocCurveDisplay(fpr=fpr, tpr=tpr).plot()

    prec, recall, _ = precision_recall_curve(labels, scores, pos_label=positive_label)
    pr_display = PrecisionRecallDisplay(precision=prec, recall=recall).plot()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))

    roc_display.plot(ax=ax1)
    pr_display.plot(ax=ax2)

    plt.savefig(output_file)


if __name__ == '__main__':
    cli()
