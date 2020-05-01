__author__ = "Max Schubach"
__copyright__ = "Copyright 2020, Max Schubach"
__email__ = "Max Schubach"
__license__ = "MIT"


import os

import h5py
import numpy as np
from itertools import product


f = h5py.File(snakemake.input[0], 'r')

# must be saved under conv1/conv1/kernel:0 I have no idea if this will generalize
kernels = np.array(f[snakemake.params[0]][snakemake.params[0]]['kernel:0'])
print(kernels.shape)
# function to save a motif
def writeMotif(file, motif, name):
    motif = motif.T
    file.write("%s\n" % (name))
    file.write("A: "+" ".join(list(map(str,motif.tolist()[0]))) + "\n")
    file.write("C: "+" ".join(list(map(str,motif.tolist()[1]))) + "\n")
    file.write("G: "+" ".join(list(map(str,motif.tolist()[2]))) + "\n")
    file.write("T: "+" ".join(list(map(str,motif.tolist()[3]))) + "\n")

# default use the best 100. can be set as second parameter
if len(snakemake.params) > 1:
    max_n = snakemake.params[1]
else:
    max_n = 100

# product of alphabet using the kerne size
alphabet = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
prod = product(alphabet,repeat=kernels.shape[0])
prod = np.array(list(prod))
print(prod.shape)


# output file. Not compressed. plain text
with open(snakemake.output[0], 'w') as motifsOut:
    # iterate over the different kernels
    for i in range(kernels.shape[-1]):
        res = (prod* kernels[:,:,i]).sum(2).sum(1)
        highest_n_index = res.argsort()[-max_n:][::-1]
        pwm = prod[highest_n_index].sum(0)/max_n
        writeMotif(motifsOut, pwm, "kernel:%d" % i)
