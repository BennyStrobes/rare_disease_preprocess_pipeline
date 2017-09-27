import numpy as np 
import os
import sys
import pdb
import gzip

##################################################
# Functions

# Extract ordered list of sample names
def get_sample_names(sample_info_file):
    names = []
    f = open(sample_info_file)
    # Sample ids are in the first column
    for line in f:
        line = line.rstrip()
        data = line.split()
        names.append(data[0])
    return names

# main driver to merge outlier calls
def merge_outlier_calls(input_file_root, output_file_name, num_input_files, sample_names):
    # open output handle
    t = open(output_file_name, 'w')
    # Now looop through all nodes of input_file_root
    for input_file_number in range(num_input_files):
        # Get filename of the input_file_number'th file
        input_file_name = input_file_root + '_' + str(input_file_number) + '_' + str(num_input_files) + '_emperical_pvalue.txt'
        # Parse through input_file_name
        head_count = 0  # for header
        f = open(input_file_name)
        for line in f:
            line = line.rstrip()
            data = np.asarray(line.split())
            if head_count == 0:
                head_count = head_count + 1
                indices = []
                indices.append(0)  # print row label 
                for i,ele in enumerate(data):
                    if ele in sample_names:
                        indices.append(i)
                if len(indices) != len(sample_names) + 1:
                    print('FATAL ERROR')
                    pdb.set_trace()
                if input_file_number == 0:
                    t.write('\t'.join(data[indices]) + '\n')
                continue
            t.write('\t'.join(data[indices]) + '\n')
        f.close()
    t.close()








##################################################
# Load data

input_file_root = sys.argv[1]  # Prefix to files produced by outlier calling. There are $num_input_files of these. Go is to merge them 
output_file_name = sys.argv[2]  # Output file
num_input_files = int(sys.argv[3])  # Number of files to merge (how many parrallel nodes were ran)
sample_info_file = sys.argv[4]  # File containing ordered names of rare disease samples

# Extract ordered list of sample names
sample_names = get_sample_names(sample_info_file)

# main driver to merge outlier calls
merge_outlier_calls(input_file_root, output_file_name, num_input_files, sample_names)