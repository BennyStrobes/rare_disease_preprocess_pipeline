import numpy as np 
import os
import sys
import pdb
import gzip


##################################################
# Functions

# Get ordered list of the rare disease samples
def get_ordered_list_of_rare_disease_samples(sample_fastq_pairing_file):
    ordered_list = []
    # sample_fastq_pairing_file contains names of rare disease samples
    f = open(sample_fastq_pairing_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        # First column contains sample id
        ordered_list.append(data[0])
    return ordered_list

# Create a dictionary that maps from junction name (chromNum_startPos_endPos) to the line number at which that junction occurs in gtex_whole_blood_jxn_file
def create_dictionary_that_maps_from_junction_id_to_line_number(gtex_whole_blood_jxn_file):
    f = open(gtex_whole_blood_jxn_file)
    # Initialize the dictionary
    dicti = {} 
    # Loop through jxn file
    head_count = 0  # Used to skip header of the file because it does not contain junction names
    line_counter = 0  # Used to keep track of which line we are at.
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        full_jxn_id = data[0]  # contains not only junction position but also cluster name and mapped Ensamble ids
        # Extract jxn_id containing only chromNum_startPos_endPos
        jxn_id_info = full_jxn_id.split(':')  # Jxn id is sperated by ':'
        jxn_id = jxn_id_info[0] + ':' + jxn_id_info[1] + ':' + jxn_id_info[2]
        # Add jxn_id to dictionary
        dicti[jxn_id] = line_counter
        line_counter = line_counter + 1
    return dicti, line_counter

# Fill in rare_disease_jxn_mat for one sample using rail output for that sample (sample_rail_jxns_output)
def fill_in_rare_disease_jxn_mat_for_one_sample(rare_disease_jxn_mat, sample_num, sample_rail_jxns_output, jxn_id_to_line_number):
    # Create handle for rail output
    f = gzip.open(sample_rail_jxns_output)
    # loop through rail output
    for line in f:
        line = line.rstrip()
        data = line.split()
        # Extract jxn id for this line
        chrom_num = data[0].split('+')[0]
        jxn_start = str(int(data[1]) - 1) # Rail splice jxn positions and leafcutter splice jxn positions are off by 1
        jxn_end = str(int(data[2]) + 1)
        num_reads = data[4]
        jxn_id = chrom_num + ':' + jxn_start + ':' + jxn_end
        # Check to see if this junction was found in the gtex whole blood junction data
        if jxn_id in jxn_id_to_line_number:
            line_number = jxn_id_to_line_number[jxn_id]
            rare_disease_jxn_mat[line_number, sample_num] = int(num_reads)
    return rare_disease_jxn_mat

# Print appended junction matrix.
def print_gtex_rare_combined_jxn_matrix(gtex_whole_blood_jxn_file, rare_disease_jxn_mat, gtex_rare_combined_jxn_file, rare_disease_samples):
    rare_disease_jxn_mat = rare_disease_jxn_mat.astype(int)  # Convert junction mat to integers because we are working with count data
    t = open(gtex_rare_combined_jxn_file, 'w')  # Initialize handle to output file
    f = open(gtex_whole_blood_jxn_file)  # Initialize handle to input file
    head_count = 0  # Used to identify header
    line_counter = 0  # used to keep track what line we are at
    #  Loop through gtex whole blood junction matrix file
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Header
            head_count = head_count + 1
            t.write(line + '\t' + '\t'.join(rare_disease_samples) + '\n')  # Write new header
            continue
        rare_counts = rare_disease_jxn_mat[line_counter,:].astype(str)  # Convert counts for the rare disease samples to strings
        t.write(line + '\t' + '\t'.join(rare_counts) + '\n')  # Write new line
        line_counter = line_counter + 1
    f.close()
    t.close()

# Print appended junction matrix.
def print_rare_only_jxn_matrix(gtex_whole_blood_jxn_file, rare_disease_jxn_mat, rare_only_jxn_file, rare_disease_samples):
    rare_disease_jxn_mat = rare_disease_jxn_mat.astype(int)  # Convert junction mat to integers because we are working with count data
    t = open(rare_only_jxn_file, 'w')  # Initialize handle to output file
    f = open(gtex_whole_blood_jxn_file)  # Initialize handle to input file
    head_count = 0  # Used to identify header
    line_counter = 0  # used to keep track what line we are at
    #  Loop through gtex whole blood junction matrix file
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Header
            head_count = head_count + 1
            t.write(data[0] + '\t' + '\t'.join(rare_disease_samples) + '\n')  # Write new header
            continue
        rare_counts = rare_disease_jxn_mat[line_counter,:].astype(str)  # Convert counts for the rare disease samples to strings
        t.write(data[0] + '\t' + '\t'.join(rare_counts) + '\n')  # Write new line
        line_counter = line_counter + 1
    f.close()
    t.close()


# Main driver scripts for this analysis
def driver(gtex_whole_blood_jxn_file, sample_fastq_pairing_file, rail_rna_dir, gtex_rare_combined_jxn_file, rare_only_jxn_file):
    # Get ordered list of the rare disease samples
    rare_disease_samples = get_ordered_list_of_rare_disease_samples(sample_fastq_pairing_file)

    # Create a dictionary that maps from junction name (chromNum_startPos_endPos) to the line number at which that junction occurs in gtex_whole_blood_jxn_file
    # Also learn number of lines (junctions) in junction matrix
    jxn_id_to_line_number, num_lines = create_dictionary_that_maps_from_junction_id_to_line_number(gtex_whole_blood_jxn_file)

    # Initialize junction matrix for rare disease samples (numJxns X numSamples)
    rare_disease_jxn_mat = np.zeros((num_lines,len(rare_disease_samples)))

    # Loop through each of rare_disease_samples and fill in rare_disease_jxn_mat
    for sample_num, rare_disease_sample in enumerate(rare_disease_samples):
        # Path to rail jxn output file
        sample_rail_jxns_output = rail_rna_dir + rare_disease_sample + '_data/cross_sample_results/first_pass_junctions.tsv.gz'
        # Fill in rare_disease_jxn_mat for one sample using rail output for that sample (sample_rail_jxns_output)
        rare_disease_jxn_mat = fill_in_rare_disease_jxn_mat_for_one_sample(rare_disease_jxn_mat, sample_num, sample_rail_jxns_output, jxn_id_to_line_number)
    
    # Print appended junction matrices
    print_gtex_rare_combined_jxn_matrix(gtex_whole_blood_jxn_file, rare_disease_jxn_mat, gtex_rare_combined_jxn_file, rare_disease_samples)
    print_rare_only_jxn_matrix(gtex_whole_blood_jxn_file, rare_disease_jxn_mat, rare_only_jxn_file, rare_disease_samples)

##################################################
# Input data

gtex_whole_blood_jxn_file = sys.argv[1]  # This is an input file that contains junction matrix for gtex whole blood samples
sample_fastq_pairing_file = sys.argv[2]  # Input file containing all sample ids
rail_rna_dir = sys.argv[3]  # Input directory containing results of running rail-rna
gtex_rare_combined_jxn_file = sys.argv[4]  # Output file path. File is a junction matrix containing samples from both gtex whole blood and the rare disease samples
rare_only_jxn_file = sys.argv[5]  # Output file path. File is a junction matrix containing samples from both gtex whole blood and the rare disease samples


driver(gtex_whole_blood_jxn_file, sample_fastq_pairing_file, rail_rna_dir, gtex_rare_combined_jxn_file, rare_only_jxn_file)