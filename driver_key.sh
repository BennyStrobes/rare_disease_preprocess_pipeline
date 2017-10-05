
##########################################################################
# Input Files
##########################################################################
# Directory containing fastq file for each sample (now removed). Basically a place holder
fastq_directory="/srv/scratch/bstrobe1/rare_variant/input_data/fastq/"

# File that contains conversion from sample names (column 1) to instituition ids (column 2; used in vcf files) to case control status (column 3)
sample_names_conversion_file="/srv/scratch/bstrobe1/rare_variant/input_data/sample_names.txt"


# Absolute location of bowtie indices and fasta (hg38)
# This is required for Rail-RNA
bowtie_prefix="/srv/scratch/bstrobe1/rare_variant/genomes/bowtie_index/Homo_sapiens/UCSC/hg38/Sequence/BowtieIndex/genome"

# Absolute location of bowtie 2 indices and fasta (hg38)
# This is required for Rail-RNA
bowtie2_prefix="/srv/scratch/bstrobe1/rare_variant/genomes/bowtie_index/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"

# Absolute location of file created in gtex outlier splicing analysis on hhpc.
# Path to file on hhpc was: /scratch1/battle-fs1/bstrober/rare_variants/rare_splice/cluster_level/standard_outlier_calling/outlier_calling/clusters_filter/Whole_Blood_Analysis_hg19_filtered_xt_reclustered_gene_mapped.txt
# We are going to basically add columns (samples to this file). Without adding /deleting junctions
gtex_whole_blood_jxn_file="/srv/scratch/bstrobe1/rare_variant/input_data/Whole_Blood_Analysis_hg19_filtered_xt_reclustered_gene_mapped2.txt"

# Absolute location of file created in gtex outlier splicing analysis on hhpc.
# Path to file on hhpc was: /scratch1/battle-fs1/bstrober/rare_variants/rare_splice/cluster_level/standard_outlier_calling/outlier_calling/clusters_filter/Whole_Blood_Analysis_hg19_filtered_xt_reclustered_gene_mapped.txt
# Same as $gtex_whole_blood_jxn_file except contains all SRA samples that are predicted to be from Whole Blood
sra_whole_blood_jxn_file="/srv/scratch/bstrobe1/rare_variant/input_data/Whole_Blood_Analysis_sra_background_jxn_mat_hg19_reordered.txt.gz"

# Absolute location of directory containing vcf file for a subset of our CGS samples (missing later 18)
cgs_vcf_input_dir="/mnt/lab_data/montgomery/lfresard/CGS_vcf_temp/"

# pLI score file (provided by Laure)
pli_score_file="/mnt/lab_data/montgomery/shared/ExAC/release0.3/functional_gene_constraint/forweb_cleaned_exac_r03_march16_z_data_pLI.txt"

# Gencoge v19 gene annotation file
gencode_gene_annotation_file="/srv/scratch/bstrobe1/rare_variant/input_data/gencoge_v19_gene_annotation.gtf"

# Directory containing data on 1000 genomes variants
# Used to access AF in a specific sub-population 
# Downloaded from ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/
genomes_1000_dir="/srv/scratch/bstrobe1/rare_variant/1kg_data/"

# File containing mapping from cluster id to jxn location
# Created on hhpc. Can be found there at: /scratch1/battle-fs1/bstrober/rare_variants/rare_splice/cluster_level/standard_outlier_calling/outlier_calling/clusters_filter/cluster_info.txt
cluster_info_file="/srv/scratch/bstrobe1/rare_variant/input_data/cluster_info.txt"

#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####3. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory="/srv/scratch/bstrobe1/rare_variant/tools/liftOver_x86"

# Various genomic annotations downloaded from the interned to be used to annotate our variants
# DANN: dowloaded from https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz
genomic_anno_dann_input="/srv/scratch/bstrobe1/rare_variant/input_data/genomic_annotations/DANN_whole_genome_SNVs.tsv.bgz"
# VEP: Run software from http://www.ensembl.org/info/docs/tools/vep/index.html
genomic_anno_vep_vcf_input="/srv/scratch/bstrobe1/rare_variant/input_data/genomic_annotations/merge_sample_names_SNPs.recode.vep.vcf"
# CHROM-HMM: Downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz
genomic_anno_chrom_hmm_input="/srv/scratch/bstrobe1/rare_variant/input_data/genomic_annotations/wgEncodeBroadHmmGm12878HMM.bed.gz"
# CADD: Downloaded from http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs_inclAnno.tsv.gz
genomic_anno_cadd_input="/srv/scratch/bstrobe1/rare_variant/input_data/genomic_annotations/whole_genome_SNVs_inclAnno.tsv.gz"
# PhyloP: Downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/
genomic_anno_phylop_dir="/srv/scratch/bstrobe1/rare_variant/input_data/genomic_annotations/"
genomic_anno_phylop_suffix=".phyloP100way.wigFix.gz"

##########################################################################
# Output Directories (Assumes these directories exist before script starts)
##########################################################################

# Directory used to store information on each of our samples
# Also contains one manifest file (to be used for rail-rna) for each sample. 
# Manifest files named $sample_name"_rail_manifest.txt"
sample_info_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/sample_info/"

# Directory used to store output of rail-rna
# within rail_rna_dir, there will be a directory for each sample called $sample_id"_data/"
# Within $sample_id"_data/" we have have the output of rail-rna for that sample
rail_rna_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/rail_rna_output/"

# Directory to store log files and other intermediate outputs of rail_rna
rail_rna_log_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/rail_rna_log/"

# Directory containing processed junction matrices
junctions_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/junctions/"

# Directory containing outlier calling results
outlier_calling_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/outlier_calling/"

# Directory containing variant call information
variant_calling_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/variant_calling/"

# Directory containing visualization of outlier calling
visualize_outliers_call_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/visualize_outlier_calls/"

# Directory containing quantification and results of enrichment of rare variants within splicing outlier calls
variant_outlier_enrichment_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/variant_outlier_enrichment/"

# Directory containing pre-processed data related to genomic annotations
genomic_annotation_data_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/genomic_annotation_data/"

# Directory containing genomic annotations for all variants in each individual
individual_annotation_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/individual_genomic_annotations/"


##########################################################################
# Output Files (Do not exist before program starts). 
##########################################################################

# File that contains mapping from sample name (column 1) to complete path of fastq_left (column 2) and fastq_right (column 3)
# Created by `python learn_sample_fastq_pairing.py`
sample_fastq_pairing_file=$sample_info_dir"sample_fastq_pairing.txt"

# File that contains mapping from sample name (column 1) to institution id (column 2) and case-control status (column 3)
# Created by `python learn_sample_fastq_pairing.py`
sample_info_pairing_file=$sample_info_dir"sample_info_pairing.txt"

# File that contains mapping from sample name (column 1) to institution id (column 2) and case-control status (column 3)
# Created by `python learn_sample_fastq_pairing.py`
# different from sample_info_pairing_file in that this file is filtered to only samples that we have wgs for
sample_info_pairing_gs_only_file=$sample_info_dir"sample_info_pairing_gs_only.txt"

# File that contains the junction matrix (Junctions X Samples)
# Created by 'python append_samples_to_junction_file.py'
# There are three different versions of this file thatvary by which samples are included (1. CGS cohort only. 2. CGS cohort, and GTEx Whole Blood. 3. CGS cohort, GTEx Whole Blood and SRA-whole blood predicted)
gtex_rare_combined_jxn_file=$junctions_dir"gtex_rare_combined_junctions.txt"
rare_only_jxn_file=$junctions_dir"rare_only_junctions.txt"
sra_gtex_rare_combined_jxn_file=$junctions_dir"sra_gtex_rare_combined_junctions.txt"

# File that contains outlier calling p-values (Clusters X samples)
# Created by outlier_calling_dm.py.
# There are three different versions of this file that vary by which background distribution samples are included (1. CGS cohort only. 2. CGS cohort, and GTEx Whole Blood. 3. CGS cohort, GTEx Whole Blood and SRA-whole blood predicted)
gtex_rare_combined_outlier_file=$outlier_calling_dir"gtex_rare_combined__emperical_pvalue.txt"
rare_only_outlier_file=$outlier_calling_dir"rare_only__emperical_pvalue.txt"
sra_gtex_rare_combined_outlier_file=$outlier_calling_dir"sra_gtex_rare_combined__emperical_pvalue.txt"

# File that contains pli scores for each cluster
# Created by get_ordered_pli_scores.py
gtex_rare_combined_pli_score_file=$outlier_calling_dir"gtex_rare_combined_pli_score.txt"
sra_gtex_rare_combined_pli_score_file=$outlier_calling_dir"sra_gtex_rare_combined_pli_score.txt"
rare_only_pli_score_file=$outlier_calling_dir"rare_only_pli_score.txt"

# File created by variant_calling_part_1.sh
# vcf file containing all variants 
# Indels were removed and so were variants on non-autosomal chromosomes
processed_vcf_file=$variant_calling_dir"merge_sample_names_SNPs.recode.vcf"


# File created by variant_calling_part_1.sh
# Bed file formatted file that contains every variant, along with:
######1. Variant position
######2. MAF according to 1K genome EUR_ANCestry
######3. The the person's two alleles at the position
######4. If there is a nearby gene
variant_bed_file=$variant_calling_dir"all_variants_gene_mapping.txt"







##########################################################################
# Alignment/preprocess Scripts
##########################################################################

###############
# Create file that contains mapping from sample name (column 1) to complete path of fastq_left (column 2) and fastq_right (column 3)
# Also creates one manifest file (to be used for rail-rna) for each sample. 
# Manifest files named $sample_info_dir$sample_name"_rail_manifest.txt"
if false; then
python learn_sample_fastq_pairing.py $fastq_directory $sample_fastq_pairing_file $sample_info_dir $sample_names_conversion_file $sample_info_pairing_file $sample_info_pairing_gs_only_file
fi

manifest_file=$sample_info_dir$sample_id"_rail_manifest.txt"
sample_output_directory=$rail_rna_dir$sample_id"_data/"
# Only perform jxn 1st pass quantification (as done by Recount/Snaptron)
if false; then
###############
# Run rail-rna for each sample
# Takes a lot of compute power to run, so don't run all samples at the same time on DURGA
while read sample_id fq_file1 fq_file2
do
    echo $sample_id
    manifest_file=$sample_info_dir$sample_id"_rail_manifest.txt"
    sample_output_directory=$rail_rna_dir$sample_id"_data/"
    # Only perform jxn 1st pass quantification (as done by Recount/Snaptron)
    nohup sh rail_wrapper.sh $bowtie_prefix $bowtie2_prefix $manifest_file $sample_output_directory $rail_rna_log_dir $liftover_directory &>/dev/null &
done<$sample_fastq_pairing_file


###############
# Create junction matrices from both GTEx whole blood samples and rare disease samples
# Create three output files using same junctions as found in $gtex_whole_blood_junction_file:
#     1. Containing both gtex whole blood samples and rare disease samples
#     2. Only rare disease samples
#     3. gtex,sra, and rare
# Takes about 15 minutes to run
python append_samples_to_junction_file.py $gtex_whole_blood_jxn_file $sample_fastq_pairing_file $rail_rna_dir $sra_whole_blood_jxn_file $gtex_rare_combined_jxn_file $rare_only_jxn_file $sra_gtex_rare_combined_jxn_file
fi









##########################################################################
# Outlier Calling scripts
##########################################################################

# Run Dirichlet Multinomial Outlier calling
# We run outlier calling three different times. Each time we only calculate outlier status in the CGS samples, however we vary the background distribution.

############### VERSION 1: GTEx, CGS
# Run dirichlet multinomial outlier calling on our junction matrices
# Junction matrices from combined whole blood gtex and rare-cohort (cgs) 
max_dm_junctions="20"
jxn_filter_method="ignore_genes"
node_number="0"
total_nodes="15"
covariate_regression_method="none"
min_reads_per_individual="5"
min_individuals_per_gene="50"

if false; then
# Parrallelize outlier calling acroos $total_nodes
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    nohup python3 call_outlier_dm.py $gtex_rare_combined_jxn_file $outlier_calling_dir"gtex_rare_combined_" $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $min_reads_per_individual $min_individuals_per_gene &> log$node_number &
done
# Merge outlier calling results from the $total_nodes different output files
python merge_outlier_calling_dm_parrallel_runs.py $outlier_calling_dir"gtex_rare_combined_" $gtex_rare_combined_outlier_file $total_nodes $sample_info_pairing_file
# For each cluster, learn the pLI score
python get_ordered_pli_scores.py $gtex_rare_combined_jxn_file $gtex_rare_combined_outlier_file $pli_score_file $gencode_gene_annotation_file $gtex_rare_combined_pli_score_file
fi


############### VERSION 2: GTEx, CGS, SRA
# Run dirichlet multinomial outlier calling on our junction matrices
# Junction matrices from combined whole blood gtex, whole blood sra, and rare-cohort (cgs) 
max_dm_junctions="20"
jxn_filter_method="ignore_genes"
node_number="0"
total_nodes="10"
covariate_regression_method="none"
min_reads_per_individual="5"
min_individuals_per_gene="50"


if false; then
# Parrallelize outlier calling acroos $total_nodes
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    nohup python3 call_outlier_dm.py $sra_gtex_rare_combined_jxn_file $outlier_calling_dir"sra_gtex_rare_combined_" $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $min_reads_per_individual $min_individuals_per_gene &> log$node_number &
done

# Merge outlier calling results from the $total_nodes different output files
python merge_outlier_calling_dm_parrallel_runs.py $outlier_calling_dir"sra_gtex_rare_combined_" $sra_gtex_rare_combined_outlier_file $total_nodes $sample_info_pairing_file
# For each cluster, learn the pLI score
python get_ordered_pli_scores.py $sra_gtex_rare_combined_jxn_file $sra_gtex_rare_combined_outlier_file $pli_score_file $gencode_gene_annotation_file $sra_gtex_rare_combined_pli_score_file
fi


############### VERSION 3: CGS
# Run dirichlet multinomial outlier calling on our junction matrices
# Junction matrices from rare-cohort (cgs) only
max_dm_junctions="20"
jxn_filter_method="ignore_genes"
node_number="0"
total_nodes="15"
covariate_regression_method="none"
min_reads_per_individual="5"
min_individuals_per_gene="25"  # Had to decrease this number b/c there were only 36 samples. And  it was 50 (as previously)

if false; then
# Parrallelize outlier calling acroos $total_nodes
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    nohup python3 call_outlier_dm.py $rare_only_jxn_file $outlier_calling_dir"rare_only_" $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $min_reads_per_individual $min_individuals_per_gene &> log_rare_only_$node_number &
done
# Merge outlier calling results from the $total_nodes different output files
python merge_outlier_calling_dm_parrallel_runs.py $outlier_calling_dir"rare_only_" $rare_only_outlier_file $total_nodes $sample_info_pairing_file
# For each cluster, learn the pLI score
python get_ordered_pli_scores.py $rare_only_jxn_file $rare_only_outlier_file $pli_score_file $gencode_gene_annotation_file $rare_only_pli_score_file
fi



##########################################################################
# Outlier calling Visualization
##########################################################################
# Perform exploratory analysis of the outlier calls in the cohort (w/o genetic information)
# Compare outlier calling distributions in cases and controls
if false; then
Rscript visualize_outlier_calls.R $visualize_outliers_call_dir $gtex_rare_combined_outlier_file $rare_only_outlier_file $sra_gtex_rare_combined_outlier_file $sample_info_pairing_file $gtex_rare_combined_pli_score_file $rare_only_pli_score_file $sra_gtex_rare_combined_pli_score_file
fi




##########################################################################
# Variant extraction/processing scripts
##########################################################################

###############
# Make propper variant calls from $cgs_vcf_input_dir 
# The files in $cgs_vcf_input_dir have already gone through the MedGap 2.4 (stanford's) pipeline to remove low quality variants
# This script filters out indels and restricted to autosomal variants
# It also computes the AF of each of the variants according to 1K genomes European Ancestry individuals
if false; then
sh variant_calling_part_1.sh $cgs_vcf_input_dir $variant_calling_dir $sample_info_pairing_gs_only_file $genomes_1000_dir $gencode_gene_annotation_file
fi

##############
# Look for enrichment of rare variants within outlier calls (do this for each of our outlier calling versions (ie background distributions))
# We will only consider variants with in $distance bp from splice junctions
distance="8"
if false; then
version="gtex_rare_combined"
sh variant_enrichment_driver.sh $variant_bed_file $cluster_info_file $variant_outlier_enrichment_dir $distance $sample_info_pairing_gs_only_file $gtex_rare_combined_outlier_file $version &

version="sra_gtex_rare_combined"
sh variant_enrichment_driver.sh $variant_bed_file $cluster_info_file $variant_outlier_enrichment_dir $distance $sample_info_pairing_gs_only_file $sra_gtex_rare_combined_outlier_file $version &

version="rare_only"
sh variant_enrichment_driver.sh $variant_bed_file $cluster_info_file $variant_outlier_enrichment_dir $distance $sample_info_pairing_gs_only_file $rare_only_outlier_file $version &
fi





##########################################################################
# Feature generation
##########################################################################


sh feature_extraction_driver.sh $genomic_annotation_data_dir $individual_annotation_dir $genomic_anno_dann_input $genomic_anno_vep_vcf_input $genomic_anno_chrom_hmm_input $genomic_anno_cadd_input $genomic_anno_phylop_dir $genomic_anno_phylop_suffix $gencode_gene_annotation_file $variant_bed_file $processed_vcf_file



