# Iterate though motifs, and pick the sequence from test set which maximally activates each.
# Place each of these sequences in the middle of a sequence.
# We create a test set of sequences, the i^th sequence represents the motif which
# maximally activates the i^th filter.

# The remainder of the sequences represent interaction effects.
# For the i^th filter and j^th filter (j > i), we predict the output with the motif
# for filter j placed up to 10 base pairs before the motif for filter i, and up to
# 10 base pairs after.

def maximal_activate_filter(model_file, sample_hdf5_file):
    torch_cmd = 'basset_motifs_predict.lua %s %s %s' % (model_file, sample_hdf5_file, 'tmp.hdf5')
    subprocess.call(torch_cmd, shell=True)

    model_hdf5_in = h5py.File(options.model_hdf5_file, 'r')
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

def get_filter_scores(mxSequences):
	num_filters = len(mxSequences)
	len_filter = len(mxSequences[0])
	testSeqs = []
	for seq in mxSequences:
		testSeqs.append('N'*300 + seq + 'N'*(300-len_filter))

	for i in range(num_filters):
		for j in range(i+1, num_filters):
			for k in range(10):
				x = 300 - len_filter - k
				seq = 'N'*x + mxSequences[j] + 'N'*k + mxSequences[i] + 'N'*(300-len_filter)
				testSeqs.append(seq)
			for k in range(10):
				x = 300 + len_filter + k
				seq = 'N'*300 + mxSequences[i] + 'N'*k + mxSequences[j] + 'N'*(600-x)
				testSeqs.append(seq)
	return testSeqs
	

def main():
    usage = 'usage: %prog [options] <model_file> <test_hdf5_file>'
    parser = OptionParser(usage)
    parser.add_option('-s', dest='sample', default=None, type='int', help='Sample sequences from the test set [Default:%default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide Basset model file and test data in HDF5 format.')
    else:
        model_file = args[0]
        test_hdf5_file = args[1]

	MIDPOINT = 300

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
        sample_hdf5_file = '%s/sample.h5' % options.out_dir
        sample_hdf5_out = h5py.File(sample_hdf5_file, 'w')
        sample_hdf5_out.create_dataset('test_in', data=seq_vecs)
        sample_hdf5_out.close()

        # update test HDF5
        test_hdf5_file = sample_hdf5_file

    print "Finished creating sample file"
    # convert to letters
    # seqs = dna_io.vecs2dna(seq_vecs)


    #################################################################
    # Torch predict
    #################################################################
    if options.model_hdf5_file is None:
        print "No model hdf5 file specified"
        options.model_hdf5_file = '%s/model_out.h5' % options.out_dir
        torch_cmd = 'basset_motifs_predict.lua %s %s %s' % (model_file, test_hdf5_file, options.model_hdf5_file)
        print torch_cmd
        subprocess.call(torch_cmd, shell=True)



if __name__ == '__main__':
	main()