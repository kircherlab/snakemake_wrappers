import click
import numpy as np
import json
import csv
import gzip
import math

import copy

import tensorflow as tf

from seqiolib import Sequence, Interval, Variant, Encoder, VariantType
from seqiolib import utils


from pyfaidx import Fasta
import pybedtools


# options
@click.command()
@click.option('--regions',
              'regions_file',
              required=True,
              multiple=True,
              type=click.Path(exists=True, readable=True),
              help='Region file in BED format.')
@click.option('--model',
              'model_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Tensorflow model in json format.')
@click.option('--weights',
              'weights_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Model weights in hdf5 format.')
@click.option('--reference',
              'reference_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Reference sequence in FASTA format (indexed).')
@click.option('--genome',
              'genome_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Genome file of the reference with lengths of contigs.')
@click.option('--altminusref/--refminusalt',
              'altMinusRef',
              default=False,
              show_default=True,
              help='Creating delta by alt minus ref or ref minus alt. default: altminusref')
@click.option('--edge',
              'edge',
              type=int,
              default=0,
              show_default=True,
              help='Edge of the input sequence not used for in-silico mutagenesis')
@click.option('--output',
              'output_file',
              type=click.Path(writable=True),
              help='Output file with predictions in tsv.gz format.')
def cli(regions_file, model_file, weights_file, reference_file, genome_file, altMinusRef, edge, output_file):

    strategy = tf.distribute.MirroredStrategy()

    def loadAndPredict(sequences, model, variants=None):
        X = []
        i = 0
        for sequence in sequences:
            if (variants is not None):
                sequence.replace(variants[i])
            X.append(Encoder.one_hot_encode_along_channel_axis(
                sequence.getSequence()))
            i += 1
        prediction = model.predict(np.array(X))
        return(prediction)

    def extendIntervals(intervals, region_length, edge, genome_file):
        output = []
        for interval in intervals:
            extend = (region_length + edge*2) % region_length
            left = math.ceil((extend-1)/2)
            right = math.floor((extend-1)/2)
            output = output + list(map(pybedtoolsIntervalToInterval,
                                       pybedtools.BedTool([interval]).slop(r=right, l=left, g=str(genome_file))))
        return(output)

    def tilingIntervals(intervals, regions, region_length, edge):
        output = []
        for i, interval in enumerate(intervals):
            if regions[i].isReverse():
                interval = Interval(interval.contig, interval.end(), interval.start())
            output = output + interval.tiling(length=region_length, shift=region_length-edge)
        return(output)

    def regionToPybedtoolsInterval(region):
        if (region.isReverse()):
            return(pybedtools.Interval(region.contig, region.end()-1, region.start(), strand="-"))
        else:
            return(pybedtools.Interval(region.contig, region.start()-1, region.end(), strand="+"))

    def pybedtoolsIntervalToInterval(interval_pybed):
        return(Interval(interval_pybed.chrom, interval_pybed.start+1, interval_pybed.stop))

    # load regions
    click.echo("Loading regions...")
    regions = []
    for region_file in regions_file:
        regions += utils.io.IntervalIO.getIntervals(region_file)
    for i in regions:
        print(i)
    click.echo("Found %d regions" % len(regions))
    if len(regions) == 0:
        click.echo(
            "No regions found. Writing file with header only and exiting...")
        with gzip.open(output_file, 'wt') as score_file:
            names = ["#Chr", "Pos", "Ref", "Alt"]
            score_writer = csv.DictWriter(
                score_file, fieldnames=names, delimiter='\t')
            score_writer.writeheader()
        exit(0)
    # convert to intervals (pybedtools)
    click.echo("Convert to bed tools intervals...")
    intervals = pybedtools.BedTool(list(map(regionToPybedtoolsInterval, regions)))
    for i in intervals:
        print(i)

    with strategy.scope():
        click.echo("Load model...")
        model = utils.io.ModelIO.loadModel(model_file, weights_file)

        input_length = model.input_shape[1]

        click.echo("Extend intervals to fit tiling...")
        intervals = extendIntervals(intervals, input_length, edge, genome_file)

        for i in intervals:
            print(i)

        click.echo("Tiling the interval of length %d ans shift %d" % (input_length, input_length-edge))
        intervals = tilingIntervals(intervals, regions, input_length, edge)

        for i in intervals:
            print(i)

    #     # load sequence for variants
    #     reference = Fasta(reference_file)
    #     sequences_ref = []
    #     sequences_alt = []

    #     click.echo("Load reference and try to get ref and alt.")
    #     for i in range(len(variants)):
    #         variant = variants[i]
    #         interval = intervals[i]

    #         # can be problematic if we are on the edges of a chromose.
    #         # Workaround. It is possible to extend the intreval left or right to get the correct length
    #         if (interval.length != input_length):
    #             click.echo("Cannot use variant %s because of wrong size of interval %s " % (
    #                 str(variant), str(interval)))
    #             continue

    #         sequence_ref = utils.io.SequenceIO.readSequence(
    #             reference, interval)

    #         # INDEL
    #         if (variant.type == VariantType.DELETION or variant.type == VariantType.INSERTION):
    #             # DELETION
    #             if (variant.type == VariantType.DELETION):
    #                 extend = len(variant.ref) - len(variant.alt)
    #                 if interval.isReverse():
    #                     interval.position = interval.position + extend
    #                 else:
    #                     interval.position = interval.position - extend
    #                 interval.length = interval.length + extend
    #             # INSERTION
    #             elif (variant.type == VariantType.INSERTION):
    #                 extend = len(variant.alt) - len(variant.ref)
    #                 if interval.isReverse():
    #                     interval.position = interval.position - extend
    #                 else:
    #                     interval.position = interval.position + extend
    #                 interval.length = interval.length - extend
    #             if (interval.length > 0):
    #                 sequence_alt = utils.io.SequenceIO.readSequence(
    #                     reference, interval)
    #                 sequence_alt.replace(variant)
    #                 if (len(sequence_alt.sequence) == input_length):
    #                     # FIXME: This is a hack. it seems that for longer indels the replacement does not work
    #                     sequences_alt.append(sequence_alt)
    #                     sequences_ref.append(sequence_ref)
    #                 else:
    #                     print("Cannot use variant %s because of wrong interval %s has wrong size after InDel Correction" % (
    #                         str(variant), str(interval)))
    #             else:
    #                 print("Cannot use variant %s because interval %s has negative size" % (
    #                     str(variant), str(interval)))
    #         # SNV
    #         else:
    #             sequence_alt = copy.copy(sequence_ref)
    #             sequence_alt.replace(variant)
    #             sequences_alt.append(sequence_alt)
    #             sequences_ref.append(sequence_ref)
    #     click.echo("Predict reference...")
    #     results_ref = loadAndPredict(sequences_ref, model)
    #     click.echo("Predict alternative...")
    #     results_alt = loadAndPredict(sequences_alt, model)

    # with gzip.open(output_file, 'wt') as score_file:
    #     names = ["#Chr", "Pos", "Ref", "Alt"]
    #     for task in range(results_ref.shape[1]):
    #         names += ["Task_%d_PredictionDelta" % task,
    #                   "Task_%d_PredictionRef" % task, "Task_%d_PredictionAlt" % task]
    #     score_writer = csv.DictWriter(
    #         score_file, fieldnames=names, delimiter='\t')
    #     score_writer.writeheader()
    #     for i in range(results_ref.shape[0]):
    #         out = {"#Chr": variants[i].contig, "Pos": variants[i].position,
    #                "Ref": variants[i].ref, "Alt":  variants[i].alt}
    #         for task in range(results_ref.shape[1]):
    #             out["Task_%d_PredictionDelta" % task] = results_alt[i][task] - \
    #                 results_ref[i][task] if altMinusRef else results_ref[i][task] - \
    #                 results_alt[i][task]
    #             out["Task_%d_PredictionRef" % task] = results_ref[i][task]
    #             out["Task_%d_PredictionAlt" % task] = results_alt[i][task]
    #         score_writer.writerow(out)


if __name__ == '__main__':
    cli()
