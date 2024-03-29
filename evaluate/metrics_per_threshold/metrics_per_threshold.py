import click
from click import progressbar
import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='TSV file with scores and label')
@click.option('--prediction-column',
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
              default=1,
              type=int,
              help='column name of the positive label')
@click.option('--steps',
              'steps',
              required=False,
              type=int,
              help='Maximum threshold steps done (if possible)')
@click.option('--decimals',
              'decimals',
              required=False,
              type=int,
              help='Threshold decimals for computation')
@click.option('--only-positive-thresholds/--all-thresholds',
              'useOnlyPositiveTrehshold',
              default=False,
              help='Use only positives thresholds (default: False)')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Final plot.')
def cli(input_file, score_column, label_column, positive_label, output_file, steps, decimals, useOnlyPositiveTrehshold):

    df = pd.read_csv(input_file, delimiter="\t")

    scores = df[[score_column]].iloc[:, 0]

    if (useOnlyPositiveTrehshold):
        thresholds = df[df[label_column] == positive_label][[score_column]].iloc[:, 0]
        thresholds = pd.concat([thresholds, df[[score_column]].min(), df[[score_column]].max()], ignore_index=False)
    else:
        thresholds = df[[score_column]].iloc[:, 0]
    
    if decimals:
        thresholds = thresholds.round(decimals).unique()
    else:
        thresholds = thresholds.unique()
    thresholds.sort()
    
    if steps and len(thresholds) // steps > 0:
        thresholds = np.unique(np.append(thresholds[0:len(thresholds):(len(thresholds) // steps)], thresholds[-1]))
        thresholds = np.sort(np.random.choice(thresholds, steps, replace=False))
    
    labels = df[[label_column]].iloc[:, 0]
    classes = labels.unique().tolist()
    classes.remove(positive_label)
    classes = classes + [positive_label]


    output = pd.DataFrame()

    with progressbar(thresholds, show_pos=True) as bar:
        for thresh in bar:
            scores_thresh = [1 if i >= thresh else 0 for i in scores]
            # prec, rec f1 score
            result_prec_rec = precision_recall_fscore_support(
                labels, scores_thresh, labels=[positive_label], pos_label=positive_label)
            result_prec_rec = [i[0] for i in result_prec_rec]

            # accuracy
            result_acc = [accuracy_score(labels, scores_thresh)] + [balanced_accuracy_score(labels, scores_thresh)]

            # TP TN FP FN
            result_conf = confusion_matrix(labels, scores_thresh, labels=classes).ravel().tolist()
 
            result = pd.DataFrame([result_conf + result_prec_rec + result_acc], columns=[
                "True-negatives", "False-positives", "False-negatives", "True-positives", "Precision",
                "Recall", "F1-score", "Support", "Accuracy", "Balanced-accuracy"
                ], index=[thresh])
            output = pd.concat([output,result])

    output["False-positive-rate"] = output["False-positives"]/(output["False-positives"]+output["True-negatives"])
    output["Specificity"] = output["True-negatives"]/(output["True-negatives"]+output["False-positives"])
    output["F2-score"] = (1+2 ** 2)*((output["Precision"]*output["Recall"])/((2 ** 2*output["Precision"])+output["Recall"]))

    output.index.name = 'Threshold'
    output = output[["True-negatives", "False-positives", "False-negatives", "True-positives",
                    "False-positive-rate", "Specificity", "Precision", "Recall", "F1-score", "F2-score",
                     "Accuracy", "Balanced-accuracy"]]

    output.to_csv(output_file, header=True, index=True, sep="\t", na_rep="NaN")


if __name__ == '__main__':
    cli()
