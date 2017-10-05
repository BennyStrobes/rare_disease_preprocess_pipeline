import numpy as np
import os
import sys
import pdb
import gzip


############################################
# Functions
############################################


# Create mapping from site to European ancestry AF (in 1000 genomes)
def map_site_to_allele_frequency(genomes_1000_af_file):
    # Initialize mapping dictionary
    dicti = {}

    # Open input file handle
    f = gzip.open(genomes_1000_af_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        start_pos = int(data[1])
        end_pos = int(data[2])
        # Simple error checking
        if end_pos - start_pos != 1:
            print('ASSUMPTION ERROR')
            pdb.set_trace()

        site_id = data[0] + '_' + data[1]
        eur_af = data[6]
        # Simple error checking
        if site_id in dicti:
            print('ASSUMPTION ERROR')
            pdb.set_trace()
        if len(data) != 8:
            print('ASSUMPTION ERROR')
            pdb.set_trace()
        dicti[site_id] = eur_af
    return dicti



def extract_sample_ids_from_vcf_file(vcf_file):
    f = open(vcf_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):  # header
            return np.asarray(data[9:])
        if line.startswith('#'):
            continue


def print_bed_file(vcf_file, site_to_af, output_file,index):
    f = open(vcf_file)
    t = open(output_file, 'w')
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        sample_genotype = data[9+index]
        site_id = 'chr' + data[0] + '_' + str(int(data[1]) -1)
        if site_id not in site_to_af:
            maf = '0.0'
        else:
            maf = site_to_af[site_id]
        if sample_genotype == '.':  # Site homozygous or unobserved in this sample
            continue
        allele1 = int(sample_genotype.split('/')[0])
        allele2 = int(sample_genotype.split('/')[1])
        # Error checking
        if allele1 < 0 or allele2 < 0:
            print('ASSUMPTION ERROR')
            pdb.set_trace()

        # Case 1
        if allele1 == 0:  # only single variant
            # error check
            if allele2 == 0:
                print('ASSUMPTION ERROR')
                pdb.set_trace()
            ref_allele = data[3]
            if len(ref_allele) != 1:
                print('ASSUMPTION ERROR')
                pdb.set_trace()
            alt_allele = data[4].split(',')[allele2-1]
            if len(alt_allele) != 1:
                print('Assumptino ERROR')
                pdb.set_trace()
            t.write('chr' + data[0] + '\t' + str(int(data[1]) -1) + '\t' + data[1] + '\t' + str(maf) + '\t1\t' + ref_allele + '\t' + alt_allele + '\n')
        elif allele1 > 0:
            if allele2 == 0:
                print('Assumption ERROR')
                pdb.set_trace()
            # Presumably 1 of two things can occur. 1. have same variant on both alleles. Have two different variants
            if allele1 == allele2:
                ref_allele = data[3]
                if len(ref_allele) != 1:
                    print('ASSUMPTION ERROR')
                    pdb.set_trace()
                alt_allele = data[4].split(',')[allele2-1]
                if len(alt_allele) != 1:
                    print('Assumptino ERROR')
                    pdb.set_trace()
                t.write('chr' + data[0] + '\t' + str(int(data[1]) -1) + '\t' + data[1] + '\t' + str(maf) + '\t1\t' + ref_allele + '\t' + alt_allele + '\n')
            elif allele2 != allele1: 
                a1 = data[4].split(',')[allele1-1]
                a2 = data[4].split(',')[allele2-1]
                if len(a1) != 1 or len(a2) != 1:
                    print('ASSUMPTION ERRORR!!')
                    pdb.set_trace()
                t.write('chr' + data[0] + '\t' + str(int(data[1]) - 1) + '\t' + data[1] + '\t' + str(maf) + '\t2\t' + a1 + '\t' + a2 + '\n')
    t.close()







############################################
# Input Data
############################################


genomes_1000_af_file = sys.argv[1]  # File containing af of our variants in the 1000 genomes superpopulations
vcf_file = sys.argv[2]  # File containing sites of interest in our cohort
output_dir = sys.argv[3]


# Create mapping from site to European ancestry AF (in 1000 genomes)
site_to_af = map_site_to_allele_frequency(genomes_1000_af_file)
ordered_sample_ids = extract_sample_ids_from_vcf_file(vcf_file)
for index, sample_id in enumerate(ordered_sample_ids):
    print(index)
    output_file = output_dir + sample_id + '_SNPs.bed'
    print_bed_file(vcf_file, site_to_af, output_file, index)