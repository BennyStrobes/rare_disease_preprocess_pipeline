import numpy as np
import os
import sys
import pdb

################## Scripts

def learn_dictionary_from_inst_to_sample(sample_info_pairing_file):
    f = open(sample_info_pairing_file)
    dicti = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        dicti[data[1]] = data[0]
    return dicti

def print_converter(input_vcf_file, output_vcf_file, institution_to_sample_ids):
    f = open(input_vcf_file)
    t = open(output_vcf_file, 'w')
    for line in f:
        line = line.rstrip()
        if line.startswith('#CHRO'): # header line
            data = line.split('\t')
            # Print first couple of columns (columns with no sample labels)
            t.write('\t'.join(data[0:9]))
            institution_ids = data[9:]
            # convert institution id one by one
            for institution_id in institution_ids: 
                t.write('\t' + institution_to_sample_ids[institution_id])
            t.write('\n')
        else: # just a normal line
            t.write(line + '\n')
    f.close()
    t.close()

################## Input data
input_vcf_file = sys.argv[1]  # Vcf file that has institution ids as column labels
output_vcf_file = sys.argv[2]  # Vcf file that has sample ids as column labels
sample_info_pairing_file = sys.argv[3] # File that has conversion from institution ids to sample ids

# Create dictionary with conversion from institution to sample ids
institution_to_sample_ids = learn_dictionary_from_inst_to_sample(sample_info_pairing_file)

# Now re-print vcf file with sample ids
print_converter(input_vcf_file, output_vcf_file, institution_to_sample_ids)

print('Conversion from institution ids to sample ids is complete')