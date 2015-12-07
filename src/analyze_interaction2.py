from optparse import OptionParser
import h5py
import numpy as np

import dna_io

import matplotlib.pyplot as plt

##### IMPORTANT: CHECK THESE ARE RIGHT #####
NUM_SEQS = 12
NUM_BASE_SEQS = 10

def parse_interaction_scores_hdf5(scores_hdf5_file):
    ### single motif scores = [activation for each filter][base_seq]
    ### paired motif scores = activation for [filter1][filter2][base_seq]
    ### seqs = motif that is used to represent each filter

    def get_base_seqs(seqs):
        return seqs[:NUM_BASE_SEQS]

    def get_motifs(seqs):
        motif_seqs = []
        z = NUM_BASE_SEQS
        for i in range(NUM_SEQS):
            motif_seqs.append(seqs[z][300:319])
            for j in range(NUM_BASE_SEQS):
                z += 1
        return motif_seqs

    def get_base_scores(preds):
        base_scores = []
        for i in range(NUM_BASE_SEQS):
            base_scores.append(preds[i])
        return base_scores

    def get_single_motif_scores(preds):
        single_motif_scores = []
        z = NUM_BASE_SEQS
        for i in range(NUM_SEQS):
            for j in range(NUM_BASE_SEQS):
                single_motif_scores[i][j] = preds[z]
                z += 1

    def get_paired_motif_scores(preds):
        paired_motif_scores = []
        for i in range(num_seqs):
            paired_motif_scores.append([])
            for j in range(num_seqs):
                paired_motif_scores[i].append([False] * NUM_BASE_SEQS)

        z = NUM_BASE_SEQS + NUM_BASE_SEQS * NUM_SEQS
        for i in range(NUM_SEQS):
            for j in range(NUM_SEQS):
                for k in range(NUM_BASE_SEQS):
                    paired_motif_scores[i][j][k] = preds[z]
                    z += 1
        return paired_motif_scores

    ### Read in file
    scores_hdf5_in = h5py.File(scores_hdf5_file, 'r')
    preds = np.array(scores_hdf5_in['preds'])
    seq_vecs = scores_hdf5_in['seqs']
    seqs = dna_io.vecs2dna(seq_vecs)
    scores_hdf5_in.close()

    ### Make sure global variables are set properly
    assert(NUM_BASE_SEQS + NUM_BASE_SEQS * NUM_SEQS + NUM_BASE_SEQS * NUM_SEQS * NUM_SEQS == len(preds))

    base_seqs = get_base_seqs(seqs)
    base_scores = get_base_scores(preds)
    single_motif_scores = get_single_motif_scores(preds)
    paired_motif_scores = get_paired_motif_scores(preds)
    motif_seqs = get_motifs(seqs)

    return (base_seqs, base_scores, single_motif_scores, paired_motif_scores, motif_seqs)

def main():
    usage = 'usage: %prog [options] <interaction_scores_hdf5_file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()
    scores_hdf5_file = args[0]

    (base_seqs, base_scores, single_motif_scores, paired_motif_scores, motif_seqs) = parse_interaction_scores_hdf5(scores_hdf5_file)
    significant_interactions = []
    for i in range(NUM_SEQS):
        for j in range(NUM_SEQS):
            bool_interaction_positive = []
            for k in range(NUM_BASE_SEQS):
                ij_score = paired_motif_scores[i][j][k][0] # 0 at the end to take first cell type
                expected_ij = base_scores[k][0] + (single_motif_scores[i][k][0] - base_scores[k][0]) + (single_motif_scores[j][k][0] - base_scores[k][0])
                bool_interaction_positive.append(int(ij_score > expected_ij))
            if sum(bool_interaction_positive) >= 9 or sum(bool_interaction_positive) <= 1:
                significant_interactions.append((i,j))

    print significant_interactions

if __name__ == '__main__':
    main()