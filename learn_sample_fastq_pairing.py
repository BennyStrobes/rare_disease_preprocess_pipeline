import numpy as np
import os
import sys
import pdb


# This script loops through each file in $fastq_directory to learn mapping from sample name to its two corresponding fastq files
# Also creates one manifest file (to be used for rail-rna) for each sample. 

fastq_directory = sys.argv[1]  # Directory containing all of the fastq files (2 per sample)
sample_fastq_pairing_output_file = sys.argv[2]  # Output file
sample_info_dir = sys.argv[3] # Directory to contain each of manifest files

# Initialize dictionary that will provide mapping from sample_name to fastq files
sample_names = {}
# Loop through files in directory
for file_name in os.listdir(fastq_directory):
    # Extract sample name from file
    sample_name = file_name.split('_')[0]

    # If sample_name has not been seen before, initialize it in dictionary
    if sample_name not in sample_names:
        sample_names[sample_name] = []
    # Add sample_name-file_name pair to dictionary
    sample_names[sample_name].append(fastq_directory + file_name)

# Open handle to write to output file
t = open(sample_fastq_pairing_output_file, 'w')
# Print one line for each sample
for sample_name in sorted(sample_names.keys()):
    # Open handle for sample's manifest file
    t_manifest = open(sample_info_dir + sample_name + '_rail_manifest.txt', 'w')

    t.write(sample_name)
    # Extract array of files
    file_array = sample_names[sample_name]
    # Simple error check
    if  len(file_array) != 2:
        print('ERROR MADE: Must check')
        quit()
    # Find fastq file corresponding to left (R1) reads
    for filer in file_array:
        if 'R1.trimm' in filer:
            t.write('\t' + filer)
            t_manifest.write(filer + '\t0\t')
    # Find fastq file corresponding to right (R1) reads
    for filer in file_array:
        if 'R2.trimm' in filer:
            t.write('\t' + filer + '\n')
            t_manifest.write(filer + '\t0\t')

    t_manifest.write(sample_name + '\n')
    t_manifest.close()

t.close()
print('Done with learning mapping from sample to fastq files')
