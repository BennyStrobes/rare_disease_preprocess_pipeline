import numpy as np 
import os
import sys
import pdb

#Extract indiviudals for which we have both WGS and RNA-seq (in specified tissue)
#@param vcf_individuals : file_name that contains a list of individuals that we have WGS for
#@param rna_seq_sampls : file_name that contains a list of RNA-seq samples for the specified tissue
#@return union_indi : dictionary of individuals that we have both WGS and RNA-seq for in this tissue
def get_wgs_rna_seq_individuals(vcf_individuals,outlier_file):
    #extract vcf individuals
    vcf_indi = {}
    f = open(vcf_individuals)
    for line in f:
        line = line.rstrip()
        data = line.split()
        indi = data[0]
        vcf_indi[indi] = 1
    union_indi = {}
    f.close()
    #extract rna seq individuals
    #as well as compute rna seq individuals that are also in WGS
    f = open(outlier_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            indiz = data[1:]
            for indi in indiz:
                if indi in vcf_indi:
                    union_indi[indi] = 1
            continue
    return union_indi




def get_tissues(tissue_file):
    f = open(tissue_file)
    arr = []
    for line in f:
        line = line.rstrip()
        arr.append(line)
    return arr

#For given gene, extract top individual and corresponding pvalue
#@param pvalz : vector of p-values for all individuals
#@param position_to_indi : mapping from array index to individual id
#@return outlier_indi : Individual with smallest pvalue fror this gene
#@return outlier_pval : pvalue of outlier_indi
def extract_outliers_top_individual(pvalz,position_to_indi):
    outlier_pval = 1.1
    for position,pval in enumerate(pvalz):
        if pval < outlier_pval:
            outlier_pval = pval
            outlier_indi = position_to_indi[position]
    if outlier_pval == 1.1:
        print('erorororo')
        pdb.set_trace()
    return outlier_indi,outlier_pval


#extract outliers. Make object called gene_struct that stores enrichment information
#gene_struct has keys that are genes. Values are dictionaries with 2 keys: 'rv_indi' and 'outlier_indi'. Values of each of those keys is a dictionary list corresponding to those individuals
#@param outlier_file : file_name containing outlier calls
#@param union_indi : dictionary that contains list of all individuals with rna-seq and WGS
#@param thresh : nominal threshold (without bonferonni correction)
#@param outlier_calling_version : whether we take only top individual per gene ('top') or all significant outliers per gene ('all')
#@return gene_struct :object described above
def extract_outliers(outlier_file, union_indi, thresh, outlier_calling_version):
    #loop through outlier file to start filling in gene_struct object
    f = open(outlier_file)
    gene_struct = {}
    counter = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('CLUSTER_ID'): #ignore header
            #create mapping from index to individual id
            position_to_indi = {}
            for i,val in enumerate(data[1:]):
                position_to_indi[i] = val
            continue
        gene_id = data[0]
        #add key (gene) to gene_struct object
        gene_struct[gene_id] = {}
        gene_struct[gene_id]['outlier_indi'] = {}
        gene_struct[gene_id]['rv_indi'] = {}
        if outlier_calling_version == 'top':
            top_indi, line_pval = extract_outliers_top_individual(np.asarray(data[1:]).astype(float),position_to_indi)
            if top_indi in union_indi and line_pval < thresh:
                gene_struct[gene_id]['outlier_indi'][top_indi] = 1

        elif outlier_calling_version == 'all':
            pvalz = np.asarray(data[1:]).astype(float)
            for position,pval in enumerate(pvalz):
                indi = position_to_indi[position]
                if indi in union_indi and pval < thresh:
                    gene_struct[gene_id]['outlier_indi'][indi] = 1
        elif outlier_calling_version == 'null': 
            pvalz = np.asarray(data[1:]).astype(float)
            for position,pval in enumerate(pvalz):
                indi = position_to_indi[position]
                if indi in union_indi and np.random.rand() < .005:
                    gene_struct[gene_id]['outlier_indi'][indi] = 1

        '''
        #Simple test to make sure random is working
        randy = np.random.randint(len(union_indi))
        new_indi = union_indi.keys()[randy]
        gene_struct[gene_id]['outlier_indi'][new_indi] = 1
        if indi in union_indi and line_pval < thresh: #if indi has rna-seq/wgs and the pvalue is less than our threshold
            gene_struct[gene_id]['outlier_indi'][indi] = 1
        '''
    return gene_struct

#update gene_struct object to fill in 'rv_indi' dictionary
#@param variant_bed_file : file_name containing all variant calls across samples at all maf ranges
#@param gene_struct : gene_struct object
#@param union_indi : dictionary of individuals for which we have WGS and RNA-seq
#@param min_maf : lower bound on maf bin (inclusive)
#@param max_maf : upper bound on maf bin (exclusive)
#@param tissues: an array of tissues
#@return gene_struct
def extract_rv_indi(variant_bed_file,gene_structs,union_indiz,min_maf,max_maf):
    #loop through rv bed file
    f = open(variant_bed_file)
    countery = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        maf = float(data[3])
        indi = data[7]
        gene = data[8]
        ref = data[5]
        alt = data[6]
        if ref == '*' or alt == '*':
            continue
        if maf < min_maf or maf >= max_maf: #maf out of range
            continue
        if indi not in union_indiz: #don't have rnaseq for this indi
            continue
        if gene not in gene_structs: #gene not in our analysis
            continue
        gene_structs[gene]['rv_indi'][indi] = 1
    return gene_structs

#Run OR enrichment analysis and print output to output_handle
#@param tissue type
#@param gene_struct :gene_struct object
#@param output_handle : file_handle object to print to
#@param num_samp : total number of samples (individuals) in OR analysis. This is the number of individuals for which we have WGS and RNA-seq
def enrichment(namer,gene_struct,output_handle,num_samp):
    a = 0 #outliers with rare variant
    b = 0 #outliers
    c = 0 #non outliers with rare variant
    d = 0 #non outliers
    for gene in sorted(gene_struct.keys()):
        rv_indi = gene_struct[gene]['rv_indi'] #dictionary containing all individuals that have a RV for this gene
        outlier_indi = gene_struct[gene]['outlier_indi'] #dictionary containing all individuals that have a RV for this gene
        if len(outlier_indi) == 0:
            continue
        count = 0 #keep track of number of outlier individuals we have for this gene only
        count_rv = 0 #keep track of number of individuals that are outliers and have rv for this gene only
        for indi_o in outlier_indi.keys():
            count = count + 1
            if indi_o in rv_indi:
                count_rv = count_rv + 1
                a = a + 1
            b = b + 1
        num_non_outliers = num_samp - count
        temp_c = len(rv_indi) - count_rv
        temp_d =num_non_outliers
        c = c + temp_c
        d = d + temp_d
    num_frac = float(a)/float(b)
    den_frac = float(c)/float(d)
    odds_ratio = (num_frac/den_frac)
    output_handle.write(namer + '\t' + str(a) + '\t' + str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(odds_ratio) + '\n')
    return output_handle


def enrichment_analysis(union_indiz, outlier_file, variant_bed_file, output_file, thresh, min_maf, max_maf):
    gene_structs = extract_outliers(outlier_file, union_indiz, thresh, 'all')

    gene_structs = extract_rv_indi(variant_bed_file, gene_structs, union_indiz, min_maf, max_maf)
    output_handle = open(output_file, 'w')
    output_handle = enrichment(str(thresh), gene_structs, output_handle, len(union_indiz))



variant_bed_file = sys.argv[1]
output_file = sys.argv[2]
outlier_file = sys.argv[3]
vcf_individuals = sys.argv[4]
thresh = float(sys.argv[5])


union_indiz = get_wgs_rna_seq_individuals(vcf_individuals, outlier_file)

#pvalue threshold

#Rare variant only
min_maf = 0.0
max_maf = .01

enrichment_analysis(union_indiz, outlier_file, variant_bed_file, output_file, thresh, min_maf, max_maf)
