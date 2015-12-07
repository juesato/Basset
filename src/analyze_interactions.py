from optparse import OptionParser
import h5py
import numpy as np

import dna_io

import matplotlib.pyplot as plt

##### CHECK THESE ARE RIGHT #####
NUM_SEQS = 6
WINDOWS = [-50, -40, -30, -20, -10, 10, 20, 30, 40, 50]
WINDOW_SIZE = len(WINDOWS)

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

    # num_seqs = len(seqs)
    assert(NUM_SEQS + WINDOW_SIZE * NUM_SEQS * NUM_SEQS == len(preds))
    num_seqs= NUM_SEQS
    window_size = WINDOW_SIZE
    # window_size = (len(preds) - num_seqs) / (num_seqs * num_seqs)

    single_motif_scores = []
    paired_motif_scores = []
    for i in range(num_seqs):
        paired_motif_scores.append([])
        for j in range(num_seqs):
            paired_motif_scores[i].append([False] * window_size)

    for i in range(num_seqs):
        single_motif_scores.append(preds[i])

    z = num_seqs
    for i in range(num_seqs):
        for j in range(num_seqs):
            for k in range(window_size):
                paired_motif_scores[i][j][k] = preds[z]
                z += 1

    return (single_motif_scores, paired_motif_scores, seqs)

def plot_scores(motif1_name, motif2_name, scores):
    x = WINDOWS
    y = scores
    plt.figure()
    plt.scatter(x,y)
    plt.savefig('interaction_plots/plot%s_%s.png' % (motif1_name[0], motif2_name[0]))
    plt.close()

def main():
    usage = 'usage: %prog [options] <interaction_scores_hdf5_file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()
    scores_hdf5_file = args[0]

    (single_motif_scores, paired_motif_scores, seqs) = parse_filter_scores_hdf5(scores_hdf5_file)
    for i in range(NUM_SEQS):
        for j in range(NUM_SEQS):
    # for i in range(1):
    #     for j in range(1):
            a = seqs[i].strip('N')
            b = seqs[j].strip('N')
            ij_scores = paired_motif_scores[i][j] # 0 at the end to take first cell type
            cell0_scores = [ij_scores[k][0] for k in range(WINDOW_SIZE)]
            plot_scores(str(i)+'_' +a,str(j)+'_'+b,cell0_scores)

    print len(single_motif_scores)
    # print "Inserting one motif"
    # print single_motif_scores
    # print len(paired_motif_scores)
    # print "Inserting combinations of one other motif with first motif"
    # print single_motif_scores[0]
    print len(paired_motif_scores[0])
    # print seqs

if __name__ == '__main__':
    main()