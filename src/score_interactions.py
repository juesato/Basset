from optparse import OptionParser
import subprocess

import h5py
import numpy as np
import random
import dna_io

# Iterate though motifs, and pick the sequence from test set which maximally activates each.
# Place each of these sequences in the middle of a sequence.
# We create a test set of sequences, the i^th sequence represents the motif which
# maximally activates the i^th filter.

# The remainder of the sequences represent interaction effects.
# For the i^th filter and j^th filter (j > i), we predict the output with the motif
# for filter j placed up to 10 base pairs before the motif for filter i, and up to
# 10 base pairs after.

####################
# Usage:
# Use score_interactions.py and analyze_interactions.py together
#
# juesato@juesato-asus:~/Basset/src$ python score_interactions.py -s 20 ../data/models/pretrained_model.th ../data/encode_roadmap.h5 interaction_scores.h5
# juesato@juesato-asus:~/Basset/src$ python analyze_interactions.py interaction_scores.h5 
#
# Right now, this script only looks at combinations of the first 30 filters (num_filters = 30 is hard-coded)
#
####################

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
        seq_code = np.zeros((4,1,seq_len), dtype='float')
        for i in range(seq_len):
            if seq[i] in bp_code:
                seq_code[bp_code[seq[i]], 0, i] = 1
            else:
                seq_code[:,0, i] = .25

        seq_vecs_list.append(seq_code)

    seq_vecs = np.array(seq_vecs_list)
    return seq_vecs

def maximal_activate_filter(model_file, sample_hdf5_file):
    torch_cmd = 'basset_motifs_predict.lua %s %s %s' % (model_file, sample_hdf5_file, 'tmp.hdf5')
    subprocess.call(torch_cmd, shell=True)

    model_hdf5_in = h5py.File('tmp.hdf5', 'r')
    filter_weights = np.array(model_hdf5_in['weights'])
    filter_outs = np.array(model_hdf5_in['outs']) # 
    seq_vecs = model_hdf5_in['sample_seqs']
    seqs = dna_io.vecs2dna(seq_vecs)
    model_hdf5_in.close()

    # store useful variables
    num_filters = filter_weights.shape[0]
    filter_size = filter_weights.shape[2]

    # For now, let's just look at what maximizes score for the first cell
    mxSequences = []
    for f in range(num_filters):
        mx = -999
        activations = filter_outs[:,f,:]
        for i in range(len(activations)):
            for j in range(len(activations[0])):
                if activations[i][j] > mx:
                    mx = activations[i][j]
                    mxI = i
                    mxJ = j
        kmer = seqs[mxI][mxJ:mxJ+filter_size]
        mxSequences.append(kmer)

    return mxSequences

def write_model_outputs_to_file(mxSequences, model_file, outputs_file):
    bp2vec = {
        'N': [.25, .25, .25, .25],
        'A': [1,0,0,0],
        'C': [0,1,0,0],
        'G': [0,0,1,0],
        'T': [0,0,0,1]
    }

    num_filters = len(mxSequences)
    len_filter = len(mxSequences[0])

    AP1_FILTERS = [54, 91, 264, 119, 174]
    CTCF_FILTERS = [9, 136, 147]
    UNK_FILTERS = [3,5,7,10]

    testSeqs = []
    BASE_SEQ = 'GCACAGTTTGTCCGCTGGGACTCTACTCGAACTTTTAACGAAACTCGTATGGCCAGGCAGGTGGACCCATAAACTCCTCGGCAAATCCCCACCGCGCCCATTTCTATGGATATTGTGCCGCATCGAACGAAAAATGATACTGCACTCGAGTCCTAAGCTCAGACCAGCGGGTACTTTCAGTCGCACTCTGGTACTACGGTATTTGGCTTCCGATAGTAAATTATGTTGCCCAGCTGCTTCAATTGGTTTGGAAGTTCGTTGGATCGGCGGAGTATTACATTTCCGACACCAGAATACCAGACGTTTTCACACATGAGAGCGGAGTTCCACCGGACGGCTGTCTCGTGGCGAAACTGTACCTAGTGTGTCGAGCTGTACTCGGTTTATAAATTAACGCATGCGAGTGAAAGAGAGATCTGTCCCGTGAACACCCACGTAGCCCTCGAGCTGCAAAGCTTCCGCCAGATAAGCACGACGGAGGCGTGCGGGGCCACACTACCGTGTTGGGAATCGAGGGATTACCAGATCCGGTGCATCAACTGATATTCCGTACCGCGACAGTGTAATGTATCGGTAGTACCCTGCCTGATTCTCTTGCCG'

    BASE_SEQS = []
    for i in range(10):
        bps = 'ACGT'
        seq = ''
        for j in range(600):
            seq += bps[random.randint(0,3)]
        BASE_SEQS.append(seq)

    # WINDOW_SIZE = 1

    FILTER_LIST = AP1_FILTERS  + CTCF_FILTERS + UNK_FILTERS
    # FILTER_LIST = AP1_FILTERS[:-2] + CTCF_FILTERS
    # FILTER_LIST = AP1_FILTERS[:-2] + CTCF_FILTERS + UNK_FILTERS[:-1]
    # FILTER_LIST = CTCF_FILTERS                                                                                       

    for z in range(10):
        testSeqs.append(BASE_SEQS[z])

    for i in FILTER_LIST:
        for z in range(10):
            base_seq = BASE_SEQS[z]
            seq = base_seq[:300] + mxSequences[i] + base_seq[300+len_filter:]
            testSeqs.append(seq)

    # for i in range(num_filters):
        # for j in range(i, num_filters):
    for i in FILTER_LIST:
        for j in FILTER_LIST:
            # for z in range(WINDOW_SIZE):
            #     k = z*10
            #     seq = BASE_SEQ[:300-len_filter-k] + mxSequences[j] + BASE_SEQ[300-k:300] + mxSequences[i] + BASE_SEQ[300+len_filter:]
            #     testSeqs.append(seq)
            # for z in range(WINDOW_SIZE):
            for z in range(10):
                base_seq = BASE_SEQS[z]
                k = 10
                seq = base_seq[:300] + mxSequences[i] + base_seq[300+len_filter:300+len_filter+k] + mxSequences[j] + base_seq[300+2*len_filter+k:]
                testSeqs.append(seq)

    # testSeqs = np.array(testSeqs, dtype='float')
    testSeqs = one_hot(testSeqs)
    print testSeqs.shape

    sample_hdf5_file = 'mxSequences_tmp.h5'
    sample_hdf5_out = h5py.File(sample_hdf5_file, 'w')
    sample_hdf5_out.create_dataset('test_in', data=testSeqs)
    sample_hdf5_out.close()

    torch_cmd = 'basset_predict.lua %s %s %s' % (model_file, sample_hdf5_file, outputs_file)
    subprocess.call(torch_cmd, shell=True)


def main():
    usage = 'usage: %prog [options] <model_file> <test_hdf5_file>'
    parser = OptionParser(usage)
    parser.add_option('-s', dest='sample', default=None, type='int', help='Sample sequences from the test set [Default:%default]')
    (options,args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Must provide Basset model file, test data in HDF5 format, and output file name.')
    else:
        model_file = args[0]
        test_hdf5_file = args[1]
        output_file = args[2]


    #################################################################
    # load data
    #################################################################
    # load sequences
    test_hdf5_in = h5py.File(test_hdf5_file, 'r')
    seq_vecs = np.array(test_hdf5_in['test_in'])
    seq_targets = np.array(test_hdf5_in['test_out'])
    try:
        target_names = list(test_hdf5_in['target_labels'])
    except KeyError:
        target_names = ['t%d'%ti for ti in range(seq_targets.shape[1])]
    test_hdf5_in.close()


    #################################################################
    # sample
    #################################################################
    if options.sample is not None:
        # choose sampled indexes
        sample_i = np.array(random.sample(xrange(seq_vecs.shape[0]), options.sample))

        # filter
        seq_vecs = seq_vecs[sample_i]
        seq_targets = seq_targets[sample_i]

        # create a new HDF5 file
        sample_hdf5_file = 'sample_tmp.h5'
        sample_hdf5_out = h5py.File(sample_hdf5_file, 'w')
        sample_hdf5_out.create_dataset('test_in', data=seq_vecs)
        sample_hdf5_out.close()

        # update test HDF5
        test_hdf5_file = sample_hdf5_file

    print "Finished creating sample file"


    mxSequences = maximal_activate_filter(model_file, sample_hdf5_file)
    print "Motifs which maximally activate each filter:", mxSequences
    write_model_outputs_to_file(mxSequences, model_file, output_file)

    #################################################################
    # Torch predict
    #################################################################

    # if options.model_hdf5_file is None:
    #     print "No model hdf5 file specified"
    #     options.model_hdf5_file = '%s/model_out.h5' % options.out_dir
    #     torch_cmd = 'basset_motifs_predict.lua %s %s %s' % (model_file, test_hdf5_file, options.model_hdf5_file)
    #     print torch_cmd
    #     subprocess.call(torch_cmd, shell=True)



if __name__ == '__main__':
    main()