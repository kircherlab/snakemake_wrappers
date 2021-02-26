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
@click.option('--verbose',
              'verbose',
              is_flag=True,
              help='Creating delta by alt minus ref or ref minus alt. default: altminusref')
@click.option('--edges',
              'edges',
              type=(int, int),
              default=(0, 0),
              show_default=True,
              help='Left and right of the input sequence not used for in-silico mutagenesis')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with predictions in tsv.gz format.')
def cli(regions_file, model_file, weights_file, reference_file, genome_file, altMinusRef, verbose, edges, output_file):

    strategy = tf.distribute.MirroredStrategy()
    
    def log(message):
        if verbose:
            click.echo(message)

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

    def extendIntervals(regions, region_length, edges, genome_file):
        # convert to intervals (pybedtools)
        click.echo("Convert to bed tools intervals...")
        intervals = pybedtools.BedTool(
            list(map(regionToPybedtoolsInterval, regions)))
        output = []
        shift = getShift(edges, region_length)
        for i, interval in enumerate(intervals):
            extended_interval = interval.length + edges[0] + edges[1]
            extend = (shift * (extended_interval//shift+1) +
                      (region_length % shift)) % extended_interval
            if regions[i].isReverse():
                right = math.ceil(extend/2)+edges[0]
                left = math.floor(extend/2)+edges[1]
            else:
                left = math.ceil(extend/2)+edges[0]
                right = math.floor(extend/2)+edges[1]
            output = output + list(map(pybedtoolsIntervalToInterval,
                                       pybedtools.BedTool([interval]).slop(r=right, l=left, g=str(genome_file))))
        return(output)

    def tilingInterval(interval, region, region_length, edges):
        if region.isReverse():
            interval = Interval(
                interval.contig, interval.end(), interval.start())
        output = interval.tiling(length=region_length,
                                 shift=getShift(edges, region_length))
        return(output)

    def regionToPybedtoolsInterval(region):
        if (region.isReverse()):
            return(pybedtools.Interval(region.contig, region.end()-1, region.start(), strand="-"))
        else:
            return(pybedtools.Interval(region.contig, region.start()-1, region.end(), strand="+"))

    def pybedtoolsIntervalToInterval(interval_pybed):
        return(Interval(interval_pybed.chrom, interval_pybed.start+1, interval_pybed.stop))

    def getShift(edge, length):
        return(length-(edge[0]+edge[1]))

    # load regions
    click.echo("Loading regions...")
    regions = []
    for region_file in regions_file:
        regions += utils.io.IntervalIO.getIntervals(region_file)
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

    reference = Fasta(reference_file)

    with strategy.scope():
        click.echo("Load model...")
        model = utils.io.ModelIO.loadModel(model_file, weights_file)

        input_length = model.input_shape[1]

        click.echo("Extend intervals to fit tiling...")
        intervals = extendIntervals(regions, input_length, edges, genome_file)

        click.echo("Tiling the interval of length %d and shift %d" %
                   (input_length, getShift(edges, input_length)))
        interval_i = 0

        with gzip.open(output_file, 'wt') as csvfile:
            writer = None

            for interval in intervals:
                log("Original: %s Extended: %s" %
                           (str(regions[interval_i]), str(interval)))

                tiled_intervals = tilingInterval(
                    interval, regions[interval_i], input_length, edges)

                log("Number of tiled intervals of interval %d: %d" %
                           (interval_i+1, len(tiled_intervals)))
                tiled_i = 0
                for tiled_interval in tiled_intervals:

                    if tiled_i % 1000 == 0:
                        log("Number of tiled intervals %d/%d of interval %d" %
                                   (tiled_i+1, len(tiled_intervals), interval_i+1))

                    log("Tiled interval %s" % tiled_interval)

                    sequence = utils.io.SequenceIO.readSequence(
                        reference, tiled_interval)
                    if tiled_interval.isReverse():
                        satMutSequences, variants = sequence.saturationMutagensis(
                            start=edges[1]+1, end=tiled_interval.length-edges[0])
                    else:
                        satMutSequences, variants = sequence.saturationMutagensis(
                            start=edges[0]+1, end=tiled_interval.length-edges[1])

                    X = []
                    for satMutSequence in satMutSequences:
                        X.append(Encoder.one_hot_encode_along_channel_axis(
                            satMutSequence.getSequence()))
                    # predict
                    prediction = model.predict(np.array(X))
                    n_tasks = np.shape(prediction)[1]

                    # initialize write once
                    if writer is None:
                        fieldnames = ["#Chr", "Pos", "Ref", "Alt"]
                        for task in range(n_tasks):
                            fieldnames += ["Task_%d_PredictionDelta" % task,
                                           "Task_%d_PredictionRef" % task, "Task_%d_PredictionAlt" % task]
                        writer = csv.DictWriter(
                            csvfile, fieldnames=fieldnames, delimiter='\t')
                        writer.writeheader()

                    for i in range(len(variants)):
                        variant = variants[i]
                        # only write out variants that ar ein original interval
                        if regions[interval_i].contains(variant):
                            toWrite = {"#Chr": variant.contig, "Pos": variant.position,
                                       "Ref": variant.ref, "Alt": variant.alt}
                            for task in range(n_tasks):
                                results = prediction[:, task].tolist()
                                toWrite["Task_%d_PredictionDelta" % task] = results[i+1] - \
                                    results[0] if altMinusRef else results[0] - \
                                    results[i+1]
                                toWrite["Task_%d_PredictionRef" %
                                        task] = results[0]
                                toWrite["Task_%d_PredictionAlt" %
                                        task] = results[i+1]

                            writer.writerow(toWrite)

                    tiled_i += 1
                interval_i += 1


if __name__ == '__main__':
    cli()
