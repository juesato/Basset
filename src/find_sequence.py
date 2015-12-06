from optparse import OptionParser
import subprocess

import h5py
import numpy as np
import random
import dna_io

def one_hot(seqs):
    bp_code = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }
    seq_vecs_list = []
    for seq in seqs:
        seq_len = len(seq)
        seq_code = np.zeros((4,1,seq_len), dtype='int8')
        for i in range(seq_len):
            if seq[i] in bp_code:
                seq_code[bp_code[seq[i]], 0, i] = 1
            else:
                seq_code[:,0, i] = .25

        # flatten and make a column vector 1 x len(seq)
        # seq_vec = seq_code.flatten()[None,:]
        print seq_code.shape
        seq_vecs_list.append(seq_code)

    seq_vecs = np.array(seq_vecs_list)
    return seq_vecs

def random_gene(length):
    bps = 'ACGT'
    seq = ''
    for i in range(length):
        a = random.randint(0,3)
        seq += bps[a]
    return seq


usage = 'usage: %prog [options] <model_file> <outputs_file>'
parser = OptionParser(usage)
(options,args) = parser.parse_args()
model_file = args[0]
outputs_file = args[1]

testSeqs = []
for i in range(30):
    testSeqs.append(random_gene(600))

testSeqs = one_hot(testSeqs)
print testSeqs.shape

sample_hdf5_file = 'randomSeqs_tmp.h5'
sample_hdf5_out = h5py.File(sample_hdf5_file, 'w')
sample_hdf5_out.create_dataset('test_in', data=testSeqs)
sample_hdf5_out.close()

torch_cmd = 'basset_predict.lua %s %s %s' % (model_file, sample_hdf5_file, outputs_file)
subprocess.call(torch_cmd, shell=True)