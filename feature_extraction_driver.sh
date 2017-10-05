
######################################################################################
# Input Data
#######################################################################################
# Directory containing pre-processed data related to genomic annotations
genomic_annotation_data_dir="$1"
# Directory containing genomic annotations for all variants in each individual
individual_annotation_dir="$2"
# Genomic annotation input data
dann_input="$3"
vep_vcf="$4"
chrom_hmm_input="$5"
cadd_input="$6"
phylop_input_dir="$7"
phylop_suffix="$8"
# Gencode v19 gene annotation file
gencode_gene_annotation_file="$9"
# File created by variant_calling_part_1.sh
# Bed file formatted file that contains every variant, along with:
######1. Variant position
######2. MAF according to 1K genome EUR_ANCestry
######3. The the person's two alleles at the position
######4. If there is a nearby gene
variant_bed_file="${10}"
# File created by variant_calling_part_1.sh
# vcf file containing all variants 
# Indels were removed and so were variants on non-autosomal chromosomes
processed_vcf="${11}"
