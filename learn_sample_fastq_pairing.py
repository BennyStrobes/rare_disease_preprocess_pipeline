import numpy as np
import os
import sys
import pdb
#RD083
#RD079

# This script loops through each file in $fastq_directory to learn mapping from sample name to its two corresponding fastq files
# Also creates one manifest file (to be used for rail-rna) for each sample. 

fastq_directory = sys.argv[1]  # Directory containing all of the fastq files (2 per sample)
sample_fastq_pairing_output_file = sys.argv[2]  # Output file
sample_info_dir = sys.argv[3] # Directory to contain each of manifest files
sample_names_conversion_file = sys.argv[4] #File that contains conversion from sample names (column 1) to instituition ids (column 2; used in vcf files) to case control status (column 3)
sample_info_pairing_file = sys.argv[5]  # output file for convresion from sample id to institution id
sample_info_pairing_gs_only_file = sys.argv[6] # output file for conversion from sample id to institution id. Excepted limited to samples that have WGS

sample_to_institution_id = {} # provides conversion from sample name to institution id
sample_to_case_control = {}  # provides conversion from sample name to case/control status
# Learn conversion from sample id to institution id & case/control status
f = open(sample_names_conversion_file)
for line in f:
    line = line.rstrip()
    data = line.split()
    sample_id = data[0]
    institution_id = data[1]
    case_control = data[2]
    sample_to_institution_id[sample_id] = institution_id
    sample_to_case_control[sample_id] = case_control

# Initialize dictionary that will provide mapping from sample_name to fastq files
sample_names = {}
# Loop through files in directory
for file_name in os.listdir(fastq_directory):
    # Extract sample name from file
    sample_name = file_name.split('_')[0]


    if sample_name == 'RD083' or sample_name == 'RD079':  # Skip these samples as they are outlier samples
        continue

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

# Print conversion from sample id to institution id & case/control status
f = open(sample_fastq_pairing_output_file)
t = open(sample_info_pairing_file,'w')
t2 = open(sample_info_pairing_gs_only_file, 'w')
for line in f:
    line = line.rstrip()
    data = line.split()
    sample_name = data[0]
    t.write(sample_name + '\t' + sample_to_institution_id[sample_name] + '\t' + sample_to_case_control[sample_name] + '\n')
    if sample_to_institution_id[sample_name].startswith('CGS'):  # Check to see if we have WGS for this sample
        t2.write(sample_name + '\t' + sample_to_institution_id[sample_name] + '\t' + sample_to_case_control[sample_name] + '\n')
t.close()


print('Done with learning mapping from sample to fastq files')



