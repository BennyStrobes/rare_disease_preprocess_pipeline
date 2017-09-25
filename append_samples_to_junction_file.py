import numpy as np 
import os
import sys
import pdb


##################################################
# Functions




##################################################
# Input data

gtex_whole_blood_jxn_file = sys.argv[1]  # This is an input file that contains junction matrix for gtex whole blood samples
sample_fastq_pairing_file = sys.argv[2]  # Input file containing all sample ids
rail_rna_dir = sys.argv[3]  # Input directory containing results of running rail-rna
gtex_rare_combined_jxn_file = sys.argv[4]  # Output file path. File is a junction matrix containing samples from both gtex whole blood and the rare disease samples
rare_only_jxn_file = sys.argv[4]  # Output file path. File is a junction matrix containing samples from both gtex whole blood and the rare disease samples


pdb.set_trace()