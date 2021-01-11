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
@click.option('--variants',
              'variants_file',
              required=True,
              multiple=True,
              type=click.Path(exists=True, readable=True),
              help='Variant file(s) to predict in TSV format.')
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
              default = False,
              show_default=True,
              help='Creating delta by alt minus ref or ref minus alt. default: altminusref')
@click.option('--file-type',
              'fileType',
              type=click.Choice(['TSV', 'VCF'], case_sensitive=False),
              default = "TSV",
              show_default=True,
              help='Variant file type.')
@click.option('--output',
              'output_file',
              type=click.Path(writable=True),
              help='Output file with predictions in tsv.gz format.')




def cli(variants_file, model_file, weights_file, reference_file, genome_file, altMinusRef, fileType, output_file):

    strategy = tf.distribute.MirroredStrategy()
    def loadAndPredict(sequences, model, variants=None):
        X=[]
        i = 0
        for sequence in sequences:
            if (variants is not None):
                sequence.replace(variants[i])
            X.append(Encoder.one_hot_encode_along_channel_axis(sequence.getSequence()))
            i += 1
        prediction = model.predict(np.array(X))
        return(prediction)

    def extendIntervals(intervals, region_length, genome_file):
        left=math.ceil((region_length-1)/2)
        right=math.floor((region_length-1)/2)
        return(list(map(pybedtoolsIntervalToInterval,intervals.slop(r=right,l=left,g=str(genome_file)))))

    def variantToPybedtoolsInterval(variant):
        return(pybedtools.Interval(variant.contig, variant.position-1, variant.position))

    def pybedtoolsIntervalToInterval(interval_pybed):
        return(Interval(interval_pybed.chrom, interval_pybed.start+1, interval_pybed.stop))
    
    if fileType == "TSV":
        fileType = utils.FileType.TSV
    elif fileType == "VCF":
        fileType = utils.FileType.VCF

    # load variants
    variants = []
    for variant_input in variants_file:
        variants += utils.VariantIO.loadVariants(variant_input, fileType=fileType)
    if len(variants) == 0:
        with gzip.open(output_file, 'wt') as score_file:
            names=["#Chr","Pos","Ref","Alt"]
            score_writer = csv.DictWriter(score_file, fieldnames=names, delimiter='\t')
            score_writer.writeheader()
        exit(0)
    # convert to intervals (pybedtools)
    intervals = pybedtools.BedTool(list(map(variantToPybedtoolsInterval,variants)))

    with strategy.scope():
        model = utils.io.ModelIO.loadModel(model_file, weights_file)

        input_length = model.input_shape[1]
        intervals = extendIntervals(intervals, input_length, genome_file)

            # load sequence for variants
        reference = Fasta(reference_file)
        sequences_ref = []
        sequences_alt = []

        for i in range(len(variants)):
            variant = variants[i]
            interval = intervals[i]

            # can be problematic if we are on the edges of a chromose.
            # Workaround. It is possible to extend the intreval left or right to get the correct length
            if (interval.length != input_length):
                print("Cannot use variant %s because of wrong size of interval %s " % (str(variant), str(interval)))
                continue

            sequence_ref = utils.io.SequenceIO.readSequence(reference,interval)

            # INDEL
            if (variant.type == VariantType.DELETION or variant.type == VariantType.INSERTION):
                # DELETION
                if (variant.type == VariantType.DELETION):
                    extend = len(variant.ref) - len(variant.alt)
                    if interval.isReverse():
                        interval.position = interval.position + extend
                    else:
                        interval.position = interval.position - extend
                    interval.length = interval.length + extend
                # INSERTION
                elif (variant.type == VariantType.INSERTION):
                    extend = len(variant.alt) - len(variant.ref)
                    if interval.isReverse():
                        interval.position = interval.position - extend
                    else:
                        interval.position = interval.position + extend
                    interval.length = interval.length - extend
                if (interval.length > 0):
                    sequence_alt = utils.io.SequenceIO.readSequence(reference,interval)
                    sequence_alt.replace(variant)
                    if (len(sequence_alt.sequence) == input_length):
                        # FIXME: This is a hack. it seems that for longer indels the replacement does not work
                        sequences_alt.append(sequence_alt)
                        sequences_ref.append(sequence_ref)
                    else:
                        print("Cannot use variant %s because of wrong interval %s has wrong size after InDel Correction" % (str(variant), str(interval)))
                else:
                    print("Cannot use variant %s because interval %s has negative size" % (str(variant), str(interval)))
            # SNV
            else:
                sequence_alt = copy.copy(sequence_ref)
                sequence_alt.replace(variant)
                sequences_alt.append(sequence_alt)
                sequences_ref.append(sequence_ref)

        results_ref = loadAndPredict(sequences_ref,model)
        results_alt = loadAndPredict(sequences_alt,model)

    with gzip.open(output_file, 'wt') as score_file:
        names=["#Chr","Pos","Ref","Alt"]
        for task in range(results_ref.shape[1]):
            names += ["Task_%d_PredictionDelta" % task, "Task_%d_PredictionRef" % task,"Task_%d_PredictionAlt" % task]
        score_writer = csv.DictWriter(score_file, fieldnames=names, delimiter='\t')
        score_writer.writeheader()
        for i in range(results_ref.shape[0]):
            out =  {"#Chr": variants[i].contig, "Pos": variants[i].position, "Ref": variants[i].ref, "Alt":  variants[i].alt}
            for task in range(results_ref.shape[1]):
                out["Task_%d_PredictionDelta" % task] = results_alt[i][task]-results_ref[i][task] if altMinusRef else results_ref[i][task]-results_alt[i][task]
                out["Task_%d_PredictionRef" % task] = results_ref[i][task]
                out["Task_%d_PredictionAlt" % task] = results_alt[i][task]
            score_writer.writerow(out)
        
if __name__ == '__main__':
    cli()
