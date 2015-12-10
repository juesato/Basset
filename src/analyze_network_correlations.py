### Looks for interaction effects in weights between filters learned by the neural network

import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
import numpy as np

from scipy.stats.stats import pearsonr
from scipy.stats import ttest_rel, ttest_ind

def get_weight_correlations(weights_matrix, filters_to_analyze):
    """Regresses across all 200 * 11 data points for each pair in filters_to_analyze
    and finds correlation coefficients"""
    num_filters = len(weights_matrix)
    num_first_layer_filters = len(weights_matrix[0])
    len_filter = len(weights_matrix[0][0])
    num_filters_to_analyze = len(filters_to_analyze)

    def overall_correlation(weights_matrix, i, j):
        weights_i = []
        weights_j = []
        for f in range(num_filters):
            weights_i += list(weights_matrix[f][i])
            weights_j += list(weights_matrix[f][j])
        (r, p) = pearsonr(weights_i, weights_j)
        return r

    def get_correlation_overall_matrix():
        overall_correlations = []
        for i in filters_to_analyze:
            print i
            overall_correlations.append([] * num_filters_to_analyze)
            for j in filters_to_analyze:
                overall_correlations[-1].append(overall_correlation(weights_matrix, i, j))
        overall_correlations = np.array(overall_correlations, dtype='float')
        return overall_correlations

    overall_correlations = get_correlation_overall_matrix()
    return overall_correlations

def plot_weight_correlations(weights_matrix, filters_to_analyze):
    num_filters = len(weights_matrix)
    num_first_layer_filters = len(weights_matrix[0])
    len_filter = len(weights_matrix[0][0])
    num_filters_to_analyze = len(filters_to_analyze)

    def count_significant_correlations(weights_matrix, i, j):   
        significant_correlations = 0    
        for f in range(num_filters):
            weights_i = weights_matrix[f][i]
            weights_j = weights_matrix[f][j]
            (r, p) = pearsonr(weights_i, weights_j)
            if p < .05:
                significant_correlations += 1
        return significant_correlations

    def overall_correlation(weights_matrix, i, j):
        weights_i = []
        weights_j = []
        for f in range(num_filters):
            weights_i += list(weights_matrix[f][i])
            weights_j += list(weights_matrix[f][j])
        (r, p) = pearsonr(weights_i, weights_j)
        return r

    def get_correlation_ct_matrix():
        num_sig_correlations = []
        for i in filters_to_analyze:
            num_sig_correlations.append([] * num_filters_to_analyze)
            for j in filters_to_analyze:
                num_sig_correlations[-1].append(count_significant_correlations(weights_matrix, i, j))

        ### Ignore diagonals
        num_sig_correlations = np.array(num_sig_correlations, dtype='float')
        for i in range(len(num_sig_correlations)):
            num_sig_correlations[i][i] = get_matrix_average_ignore_diagonal(num_sig_correlations)
        ### Subtract mean so that plot is better
        avg = np.mean(num_sig_correlations)
        num_sig_correlations -= avg
        return num_sig_correlations 

    def get_matrix_average_ignore_diagonal(matrix):
        tot = 0.0
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if i == j:
                    continue
                tot += matrix[i][j]

        return (tot / (len(matrix) * len(matrix) - len(matrix)))

    def get_correlation_overall_matrix():
        overall_correlations = []
        for i in filters_to_analyze:
            overall_correlations.append([] * num_filters_to_analyze)
            for j in filters_to_analyze:
                overall_correlations[-1].append(overall_correlation(weights_matrix, i, j))

        overall_correlations = np.array(overall_correlations, dtype='float')
        print get_matrix_average_ignore_diagonal(overall_correlations)
        for i in range(len(overall_correlations)):
            overall_correlations[i][i] = get_matrix_average_ignore_diagonal(overall_correlations)

        avg = np.mean(overall_correlations)
        overall_correlations -= (avg + .04)     
        for i in range(len(overall_correlations)):
            overall_correlations[i][i] = 0  
        return overall_correlations

    num_sig_correlations = get_correlation_ct_matrix()
    overall_correlations = get_correlation_overall_matrix()

    # Insertion for our visualization
    return overall_correlations

    within_ap1 = []
    overall = []
    random = []
    within_ctcf = []
    for i in range(5):
        for j in range(5):
            if i == j:
                continue
            within_ap1.append(overall_correlations[i][j])
    for i in range(12):
        for j in range(12):
            if (i < 5 and j < 5) or (5 <= i <= 7 and 5 <= j <= 7) or i==j:
                continue
            overall.append(overall_correlations[i][j])
    for i in range(5,8):
        for j in range(5,8):
            if i ==j:
                continue
            within_ctcf.append(overall_correlations[i][j])
    for i in range(8,12):
        for j in range(8,12):
            if i ==j:
                continue
            random.append(overall_correlations[i][j])
    (t,p) = ttest_ind(within_ap1, overall)
    print (t,p)
    (t,p) = ttest_ind(within_ctcf, overall)
    print (t,p)
    (t,p) = ttest_ind(random, overall)
    print (t,p)

    plot_filter_heat(num_sig_correlations, 'interaction_plots/network_correlations_count.png')
    plot_filter_heat(overall_correlations, 'interaction_plots/network_correlations_overall.png')


################################################################################
# plot_filter_heat
#
# Plot a heatmap of the filter's parameters.
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def plot_filter_heat(param_matrix, out_pdf):
    param_range = abs(param_matrix).max()

    plt.figure(figsize=(param_matrix.shape[1], 4))
    sns.heatmap(param_matrix, cmap='PRGn', linewidths=0.2, vmin=-param_range, vmax=param_range)
    ax = plt.gca()
    ax.set_yticklabels(range(param_matrix.shape[1],0,-1))
    ax.set_xticklabels(range(1,param_matrix.shape[1]+1))
    # ax.set_yticklabels('TGCA', rotation='horizontal', size=10)
    plt.savefig(out_pdf)
    plt.close()


def main():
    usage = 'usage: %prog [options] <weights_file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide second layer weights as HDF5 file')
    else:
        weights_file = args[0]

    f = h5py.File(weights_file, 'r')
    weights_matrix = f['weights']

    AP1_FILTERS = [54, 91, 264, 119, 174]
    CTCF_FILTERS = [9, 136, 147]
    UNK_FILTERS = [3,5,7,10]

    FILTER_LIST = AP1_FILTERS  + CTCF_FILTERS + UNK_FILTERS

    corrs = get_weight_correlations(weights_matrix, range(300))
    # corrs = get_weight_correlations(weights_matrix, range(5))

    import sys
    import csv
    writer = csv.writer(sys.stdout)
    writer.writerows(corrs)
    # analyze_weights(matrix)


if __name__ == '__main__':
    main()