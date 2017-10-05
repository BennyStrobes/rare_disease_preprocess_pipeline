import numpy as np 
import os
import sys
import pdb
import gzip
import time


def simple_debugging_checks(variant_bed_dir,input_individuals):
    f = open(input_individuals)
    indiz = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        indi = data[0]
        indiz[indi] = 1
    f.close()
    indiz2 = {}
    for file_name in os.listdir(variant_bed_dir):
        if file_name.endswith('SNPs.bed'):
            individual_id = file_name.split('_SNP')[0]
            indiz2[individual_id] = 1
    if indiz != indiz2:
        print('Error: Individuals in your individual input file do not match the individuals in your variant_bed_dir')

#Simple helper method to remove duplicated gene names at a given position in the chromosome data structure
def get_unique_genes(stringer):
    info = stringer.split(',')
    return ','.join(np.unique(info))

#Fill in chromosome object from chromosome[start-1000:end+1+1000] with the addition of this gene name
def method_1_add_gene_to_chromosome_object(chromosome,start,end,gene_name):
    #Consider gene to be all bp in gene, as well as all bp 1000 bases upstream from gene and all bp 1000 downstream from end of gene
    #Adjust start and end to do this exactly 
    new_start = start -1000
    if new_start < 0: #make sure new adjusted exon start isn't less than zero
        print('Gene Start went less than zero!!!!')
        new_start = 0
    new_end = end + 1001
    if new_end > 259250621: #make sure new adjusted end isn't greater than max chromosome size!
        print('Gene end went greater than max length of chromosome!!!')
        new_end = 259250621
    ###########################################################
    #Actually fill in the chromosome
    for pos in range(new_start,new_end+1):
        if chromosome[pos] == 'NULL': #No genes already mapped to this position
            chromosome[pos] = gene_name
        else: #at least one gene already mapped to this position
            chromosome[pos] = get_unique_genes(chromosome[pos] + ',' + gene_name)
    return chromosome


#Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
def method_1_make_chromosome_with_gene_names(chrom_num,gencode_hg19_gene_annotation_file):
    #initialize chromosome array
    chromosome = ['NULL']*259250621
    #loop through gencode file
    f = open(gencode_hg19_gene_annotation_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'): #ignore header lines
            continue
        gene_type = data[13].split('"')[1] #ie protein_coding,pseudo_gene,etc
        gene_name = data[9].split('"')[1] #ensamble id
        line_chrom_num = data[0]
        gene_part = data[2] #gene,UTR,exon,etc
        if gene_type != 'protein_coding': #limit to only protein coding genes
            continue
        if line_chrom_num != 'chr' + chrom_num: #limit to chromosome of interest
            continue
        if gene_part != 'gene': #method 1 is limited to gene only
            continue
        start = int(data[3]) #Start  of gene
        end = int(data[4]) #End (downstream) of gene
        if start > end:
            print('Error in underlying assumptions. need to reconsider everything')
        #Fill in chromosome object from chromosome[start-1000:end+1001] with the addition of this gene name
        chromosome = method_1_add_gene_to_chromosome_object(chromosome,start,end,gene_name)
    return chromosome

#extract individuals for which we have variant beds
def extract_individuals(input_individuals):
    individuals = []
    f = open(input_individuals)
    for line in f:
        line = line.rstrip()
        data = line.split()
        indi = data[0]
        individuals.append(indi)
    return sorted(individuals)


#stream individual's variant bed file. For each variant, use chromosome data structure to check which genes (if any) the variant overlaps
def stream_variant_bed(individual_bed_file,chrom_number,chromosome,out_handle,individual):
    f = open(individual_bed_file)
    chrom_number_string = 'chr' + chrom_number
    for line in f:
        line = line.rstrip()
        data = line.split()
        chrom_num = data[0]
        if chrom_num != chrom_number_string: #chromosome specific analysis
            continue
        variant_position = int(data[2])
        overlapping_genes = chromosome[variant_position]
        out_handle.write(line + '\t' + individual + '\t' + overlapping_genes + '\n')
        if variant_position - int(data[1]) != 1:
            print('DATA ASSUMPTION ERROR')
    return out_handle

#Perform analysis on each chromosome individually (for memory reasons)
def method_1_chromosome_specific_driver(variant_bed_dir,input_individuals,gencode_hg19_gene_annotation_file,out_handle,chrom_number):
    #Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
    chromosome = method_1_make_chromosome_with_gene_names(chrom_number,gencode_hg19_gene_annotation_file)
    #extract individuals for which we have variant beds
    individuals = extract_individuals(input_individuals)
    #loop through individuals
    for individual in individuals:
        individual_bed_file = variant_bed_dir + individual + '_SNPs.bed'
        #stream individual's variant bed file. For each variant, use chromosome data structure to check which genes (if any) the variant overlaps
        out_handle = stream_variant_bed(individual_bed_file,chrom_number,chromosome,out_handle,individual)
    return out_handle


#driver for method_1
def run_method_1(variant_bed_dir,input_individuals,all_variants_gene_mapping_file,gencode_hg19_gene_annotation_file):
    #Simple check to make sure the individuals in the input_individuals file are all present in the variant_bed_dir
    #Ie. Did all of the variant_extraction_nick_emily_parallel_part.sh jobs finish running
    simple_debugging_checks(variant_bed_dir, input_individuals)
    #output file handle (we are going to print while streaming variant files)
    out_handle = open(all_variants_gene_mapping_file,'w')
    #Perform analysis on each chromosome individually (for memory reasons)
    for chrom_number in range(1,23):
        print(chrom_number)
        out_handle = method_1_chromosome_specific_driver(variant_bed_dir,input_individuals,gencode_hg19_gene_annotation_file,out_handle,str(chrom_number))
        out_handle.flush()




variant_bed_dir = sys.argv[1] #Input directory that contains variant information for each individual
input_individuals = sys.argv[2] #Input file that contains list of all individuals in variant_bed_dir
all_variants_gene_mapping_file = sys.argv[3] #output File
gencode_hg19_gene_annotation_file = sys.argv[4] #gencode hg19 gene annotation file

#method_1: Limit to protein coding. A gene considered to be what gencode defines as 'gene', as well as 1 KB upstream and 1KB downstream
#if gene_mapping_method == 'method_1':
#driver for method_1
run_method_1(variant_bed_dir, input_individuals, all_variants_gene_mapping_file, gencode_hg19_gene_annotation_file)