from optparse import OptionParser
import h5py
import numpy as np

import dna_io

import matplotlib.pyplot as plt

##### CHECK THESE ARE RIGHT #####
NUM_SEQS = 12
NUM_BASE_SEQS = 10

def parse_filter_scores_multiple_base_hdf5(scores_hdf5_file):
    ### single motif scores = [activation for each filter][base_seq]
    ### paired motif scores = activation for [filter1][filter2][base_seq]
    ### seqs = motif that is used to represent each filter
    scores_hdf5_in = h5py.File(scores_hdf5_file, 'r')
    preds = np.array(scores_hdf5_in['preds'])
    seq_vecs = scores_hdf5_in['seqs']
    # print preds.shape
    # print seq_vecs.shape
    seqs = dna_io.vecs2dna(seq_vecs)
    scores_hdf5_in.close()

    assert(NUM_BASE_SEQS + NUM_BASE_SEQS * NUM_SEQS + NUM_BASE_SEQS * NUM_SEQS * NUM_SEQS == len(preds))
    num_seqs= NUM_SEQS

    base_scores = []
    single_motif_scores = []
    paired_motif_scores = []

    for i in range(num_seqs):
        single_motif_scores.append([False] * NUM_BASE_SEQS)
        paired_motif_scores.append([])
        for j in range(num_seqs):
            paired_motif_scores[i].append([False] * NUM_BASE_SEQS)

    z = 0

    for i in range(NUM_BASE_SEQS):
        base_scores.append(preds[i])
        base_seqs.append(seqs[i])
        z += 1

    motif_seqs = []
    for i in range(NUM_SEQS):
        motif_seqs.append(seqs[z][300:319])
        for j in range(NUM_BASE_SEQS):
            single_motif_scores[i][j] = preds[z]
            z += 1

    for i in range(num_seqs):
        for j in range(num_seqs):
            for k in range(NUM_BASE_SEQS):
                paired_motif_scores[i][j][k] = preds[z]
                z += 1

    return (base_seqs, base_scores, single_motif_scores, paired_motif_scores, motif_seqs)

def plot_scores(motif1_name, motif2_name, scores):
    x = WINDOWS
    y = scores
    plt.figure()
    plt.scatter(x,y)
    plt.ylim(-1,1)
    plt.savefig('interaction_plots/plot%s_%s.png' % (motif1_name[0], motif2_name[0]))
    plt.close()

def main():
    usage = 'usage: %prog [options] <interaction_scores_hdf5_file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()
    scores_hdf5_file = args[0]

    (base_seqs, base_scores, single_motif_scores, paired_motif_scores, motif_seqs) = parse_filter_scores_multiple_base_hdf5(scores_hdf5_file)
    for i in range(NUM_SEQS):
        for j in range(NUM_SEQS):
            for k in range(NUM_BASE_SEQS)
                a = seqs[i].strip('N')
                b = seqs[j].strip('N')
                ij_scores = paired_motif_scores[i][j] # 0 at the end to take first cell type
                cell0_scores = [ij_scores[k][0] for k in range(WINDOW_SIZE)]
                expected_ij = base_scores[k][0] (single_motif_scores[i][0] + single_motif_scores[j][0]) / 2.0
                plot_scores(str(i)+'_'+a, str(j)+'_'+b, cell0_scores - avg_ij)

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