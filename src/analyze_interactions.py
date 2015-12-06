from optparse import OptionParser
import h5py
import numpy as np

import dna_io

def parse_filter_scores_hdf5(scores_hdf5_file):

    ### single motif scores = [activation for each filter]
    ### paired motif scores = activation for [filter1][filter2][offset]
    ### seqs = motif that is used to represent each filter
    scores_hdf5_in = h5py.File(scores_hdf5_file, 'r')
    preds = np.array(scores_hdf5_in['preds'])
    seq_vecs = scores_hdf5_in['seqs']
    print preds.shape
    print seq_vecs.shape
    seqs = dna_io.vecs2dna(seq_vecs)
    scores_hdf5_in.close()

    num_seqs = len(seqs)
    window_size = (len(preds) - num_seqs) / (num_seqs * num_seqs)

    single_motif_scores = []
    paired_motif_scores = []
    for i in range(num_seqs):
        paired_motif_scores.append([])
        for j in range(i, num_seqs):
            paired_motif_scores[i].append([False] * window_size)

    for i in range(num_seqs):
        single_motif_scores.append(preds[i])

    z = num_seqs
    for i in range(num_seqs):
        for j in range(i, num_seqs):
            for k in range(window_size):
                paired_motif_scores[i][j][k] = preds[z]
                z += 1

    return (single_motif_scores, paired_motif_scores, seqs)

def main():
    usage = 'usage: %prog [options] <interaction_scores_hdf5_file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()
    scores_hdf5_file = args[0]

    (single_motif_scores, paired_motif_scores, seqs) = parse_filter_scores_hdf5(scores_hdf5_file)
    print len(single_motif_scores)
    print "Inserting one motif"
    print single_motif_scores
    # print len(paired_motif_scores)
    print "Inserting combinations of one other motif with first motif"
    print single_motif_scores[0]
    print len(paired_motif_scores[0])
    print seqs

if __name__ == '__main__':
    main()