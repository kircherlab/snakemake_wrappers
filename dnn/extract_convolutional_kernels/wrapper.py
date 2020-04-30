__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os

import h5py
import numpy as np


f = h5py.File(snakemake.input[0], 'r')

# must be saved under conv1/conv1/kernel:0 I have no idea if this will generalize
kernels = np.array(f[snakemake.params[0]][snakemake.params[0]]['kernel:0'])

# The question is how to normalize?I woud suggest that the kernel integers range only from -1 to 1. Otherwise we can use the min max
#kernel_conv1_normMinMax = np.interp(kernel_conv1, (kernel_conv1.min(), kernel_conv1.max()), (0, 1))
kernels_norm = np.interp(kernels, (-1, 1), (0, 1))

# function to save a motif
def writeMotif(file, motif, name):
    file.write("%s\n" % (name))
    file.write("A: "+" ".join(list(map(str,motif.tolist()[0]))) + "\n")
    file.write("C: "+" ".join(list(map(str,motif.tolist()[1]))) + "\n")
    file.write("G: "+" ".join(list(map(str,motif.tolist()[2]))) + "\n")
    file.write("T: "+" ".join(list(map(str,motif.tolist()[3]))) + "\n")

# output file. Not compressed. plain text
with open(snakemake.output[0], 'w') as motifsOut:
    # iterate over the different kernels
    for i in range(kernels_norm.shape[2]):
        writeMotif(motifsOut, np.transpose(kernels_norm[:,:,i]), "kernel:%d" % i)
