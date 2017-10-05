import numpy as np
import os
import sys
import pdb
import time


def remove_duplicates(old_string, add_on):
    if old_string == 'NULL':
        return add_on
    else:
        arr = old_string.split(',')
        if add_on in arr:
            return old_string
        else:
            return old_string + ',' + add_on


def fill_in_chromosome(chrom, posi, cluster_id, distance):
    new_start = posi - distance
    if new_start < 0:
        new_start = 0
    new_end = posi + distance
    for pos in range(new_start, new_end + 1):
        chrom[pos] = remove_duplicates(chrom[pos], cluster_id)
    return chrom


def make_chromosome(cluster_file, chrom_num, distance):
    chrom = ['NULL']*259250621
    f = open(cluster_file)
    chrom_string = 'chr' + str(chrom_num)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        cluster_id = data[0]
        jxns = data[1].split(',')
        if len(jxns) < 2:
            print('ASSUMPTION ERROR')
            pdb.set_trace()
        chromer = jxns[0].split(':')[0]
        if chromer != chrom_string:  # Skip junctions not on this chromosome
            continue
        for jxn in jxns:
            jxn_info = jxn.split(':')
            start = int(jxn_info[1])
            end = int(jxn_info[2])
            chrom = fill_in_chromosome(chrom, start, cluster_id, distance)
            chrom = fill_in_chromosome(chrom, end, cluster_id, distance)
    return chrom


def stream_variant_bed_file(cluster_chromosome, variant_bed_file, t, chrom_num):
    chrom_string = 'chr' + str(chrom_num)
    f = open(variant_bed_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        line_chrom_num = data[0]
        if line_chrom_num != chrom_string:
            continue
        maf = float(data[3])
        if maf >= .01:
            continue
        var_pos = int(data[2])

        overlapping_clusters = cluster_chromosome[var_pos]
        clusterz = overlapping_clusters.split(',')
        for cluster in clusterz:
            t.write('\t'.join(data[0:8]) + '\t' + cluster + '\n')
    return t


variant_bed_file = sys.argv[1]
cluster_file = sys.argv[2]
output_file = sys.argv[3]
distance = int(sys.argv[4])


t = open(output_file, 'w')

for chrom_num in range(1, 23):
    start = time.time()
    print(chrom_num)
    cluster_chromosome = {}
    cluster_chromosome = make_chromosome(cluster_file, chrom_num, distance)
    t = stream_variant_bed_file(cluster_chromosome, variant_bed_file, t, chrom_num)
    print((time.time() - start) / 60.0)