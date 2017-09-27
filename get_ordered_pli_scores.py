import numpy as np
import os
import sys
import pdb
import gzip




# Get mapping from clusters to gene names.
# Also get list of relevant genes
def get_mapping_from_clusters_to_gene_names(jxn_file):
    f = open(jxn_file)
    head_count = 0 # For header
    cluster_to_gene_name = {}
    gene_names = {}
    # Loop through jxn file
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        # row_label contains both cluster_id and genes the cluster overlaps
        row_label = data[0]
        row_label_info = row_label.split(':')  # relevant fields in row_label seperated by :
        cluster_id = row_label_info[3] 
        genes = row_label_info[4].split(',')
        # Add observed genes
        for gene in genes:
            gene_names[gene] = 1
        # If we've never seen this cluster before
        if cluster_id not in cluster_to_gene_name:
            cluster_to_gene_name[cluster_id] = genes
        # If we've seen this cluster before, then must append gene names onto existing array
        else:
            arr = cluster_to_gene_name[cluster_id]
            for gene in genes:
                arr.append(gene)
            cluster_to_gene_name[cluster_id] = arr
    # Remove duplicate genes from mapping
    for keys in cluster_to_gene_name.keys():
        arr = cluster_to_gene_name[keys]
        cluster_to_gene_name[keys] = np.unique(arr)
    return cluster_to_gene_name, gene_names

# Now get mapping from transcript_id to gene_name if gene is in gene_names
def get_mapping_from_transcript_id_to_gene_name(gencode_gene_annotation_file, gene_names):
    # loop through gencode file
    dicti = {}
    f = open(gencode_gene_annotation_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):  # ignore header lines
            continue
        gene_type = data[13].split('"')[1]  # ie protein_coding,pseudo_gene,etc
        gene_name = data[9].split('"')[1]  # ensamble id
        line_chrom_num = data[0]
        gene_part = data[2]  # gene,UTR,exon,etc
        if gene_name not in gene_names:  # Only consider genes we have (those in gene_names)
            continue
        if gene_part != 'transcript': # Only care about transcripts
            continue
        transcript_id = data[11].split('"')[1]
        if transcript_id not in dicti:
            dicti[transcript_id] = gene_name
        else: # We've seen this transcript before
            if gene_name != dicti[transcript_id]:
                print('ERROR!!') 
    return dicti

# get mapping from gene_name to pli score using transcript_id to gene_name mapping
# Need to consider what happens if multiple pli scores per gene (probably just take max)
def get_mapping_from_gene_name_to_pli_score(pli_score_file, transcript_id_to_gene_name):
    gene_name_to_pli_score = {}
    f = open(pli_score_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        transcript_id = data[0]
        pli_score = float(data[19])
        if transcript_id in transcript_id_to_gene_name:
            geney = transcript_id_to_gene_name[transcript_id]
            if geney not in gene_name_to_pli_score:
                gene_name_to_pli_score[geney] = pli_score
            else:
                print('ERRORORO')
                pdb.set_trace()
    return gene_name_to_pli_score

###Input data
jxn_file = sys.argv[1]  # Contains mapping from clusters to gene names
outlier_file = sys.argv[2]  # Contains ordered clusters (we will print pli scores in that order)
pli_score_file = sys.argv[3]  # Contains pli scores for transcript ids
gencode_gene_annotation_file = sys.argv[4]  # contains mapping from gene names to transcript ids
output_file = sys.argv[5]


# Get mapping from clusters to gene names.
# Also get list of relevant genes
# Jxn file row labels have the genes corresponding to each cluster. We are going to extract this mapping
cluster_to_gene_name, gene_names = get_mapping_from_clusters_to_gene_names(jxn_file)

# Now get mapping from transcript_id to gene_name if gene is in gene_names
transcript_id_to_gene_name = get_mapping_from_transcript_id_to_gene_name(gencode_gene_annotation_file, gene_names)

# get mapping from gene_name to pli score using transcript_id to gene_name mapping
# Need to consider what happens if multiple pli scores per gene (probably just take max)
gene_name_to_pli_score = get_mapping_from_gene_name_to_pli_score(pli_score_file, transcript_id_to_gene_name)

# open output handle
t = open(output_file, 'w')
t.write('clusterID\tgene_names\tmax_pli_score\n')
f = open(outlier_file)
head_count = 0
for line in f:
    line = line.rstrip()
    data = line.split()
    if head_count == 0:
        head_count = head_count + 1
        continue
    cluster_id = data[0]
    genes = cluster_to_gene_name[cluster_id]
    pli_scores = []
    for gene in genes:
        if gene in gene_name_to_pli_score:
            pli_scores.append(gene_name_to_pli_score[gene])
    if len(pli_scores) == 0:
        t.write(cluster_id + '\t' + ','.join(genes) + '\tNA\n')
    else:
        t.write(cluster_id + '\t' + ','.join(genes) + '\t' + str(np.max(pli_scores)) + '\n')

t.close()