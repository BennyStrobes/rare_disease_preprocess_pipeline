
##########################################################################
# Input Files
##########################################################################
# Directory containing fastq file for each sample
fastq_directory="/users/lfresard/CGS_fastq_temp/"

# Directory containing hg19 genome (fasta) 
genome_directory="/mnt/lab_data/montgomery/shared/genomes/"

# Absolute location of bowtie indices and fasta (hg19)
# This is required for Rail-RNA
bowtie_prefix="/mnt/lab_data/montgomery/shared/genomes/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome"

# Absolute location of bowtie 2 indices and fasta (hg19)
# This is required for Rail-RNA
bowtie2_prefix="/mnt/lab_data/montgomery/shared/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"

# Absolute location of file created in gtex outlier splicing analysis on hhpc.
# Path to file on hhpc was: /scratch1/battle-fs1/bstrober/rare_variants/rare_splice/cluster_level/standard_outlier_calling/outlier_calling/clusters_filter/Whole_Blood_Analysis_hg19_filtered_xt_reclustered_gene_mapped.txt
# We are going to basically add columns (samples to this file). Without adding /deleting junctions
gtex_whole_blood_jxn_file="/srv/scratch/bstrobe1/rare_variant/input_data/Whole_Blood_Analysis_hg19_filtered_xt_reclustered_gene_mapped2.txt"








##########################################################################
# Output Directories (Assumes these directories exist before script starts)
##########################################################################

# Directory used to store information on each of our samples
# Also contains one manifest file (to be used for rail-rna) for each sample. 
# Manifest files named $sample_name"_rail_manifest.txt"
sample_info_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/sample_info/"

# Directory used to store output of rail-rna
# within rail_rna_dir, there will be a directory for each sample called $sample_id"_data/"
# Within $sample_id"_data/" we have have the results of rail-rna for that sample
rail_rna_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/rail_rna_output/"

# Directory to store log files and other intermediate outputs of rail_rna
rail_rna_log_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/rail_rna_log/"

# Directory containing processed junction matrices
junctions_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/junctions/"

# Directory containing outlier calling results
outlier_calling_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/outlier_calling/"



##########################################################################
# Output Files (Do not exist before program starts). 
##########################################################################

# File that contains mapping from sample name (column 1) to complete path of fastq_left (column 2) and fastq_right (column 3)
# Created by `python learn_sample_fastq_pairing.py`
sample_fastq_pairing_file=$sample_info_dir"sample_fastq_pairing.txt"

# File that contains the junction matrix (Junctions X Samples) for both gtex whole blood samples and the rare disease samples
# Created by 'python append_samples_to_junction_file.py'
gtex_rare_combined_jxn_file=$junctions_dir"gtex_rare_combined_junctions.txt"

# File that contains the junction matrix (Junctions X Samples) for only the rare disease samples
# Created by 'python append_samples_to_junction_file.py'
rare_only_jxn_file=$junctions_dir"rare_only_junctions.txt"

# File that contains outlier calling p-values (Clusters X samples)
# Created by outlier_calling_dm.py.














##########################################################################
# Scripts/Functions/Analyses
##########################################################################

###############
# Create file that contains mapping from sample name (column 1) to complete path of fastq_left (column 2) and fastq_right (column 3)
# Also creates one manifest file (to be used for rail-rna) for each sample. 
# Manifest files named $sample_info_dir$sample_name"_rail_manifest.txt"
if false; then
python learn_sample_fastq_pairing.py $fastq_directory $sample_fastq_pairing_file $sample_info_dir
fi


###############
# Run rail-rna for each sample
# Takes a lot of compute power to run, so don't run all samples at the same time on DURGA
if false; then
while read sample_id fq_file1 fq_file2
do
    echo $sample_id
    manifest_file=$sample_info_dir$sample_id"_rail_manifest.txt"
    sample_output_directory=$rail_rna_dir$sample_id"_data/"
    if false; then
    # Only perform jxn 1st pass quantification (as done by Recount/Snaptron)
    sh rail_wrapper.sh $bowtie_prefix $bowtie2_prefix $manifest_file $sample_output_directory $rail_rna_log_dir &>/dev/null &
    fi
done<$sample_fastq_pairing_file
fi



###############
# Create junction matrices from both GTEx whole blood samples and rare disease samples
# Create two output files using same junctions as found in $gtex_whole_blood_junction_file:
#     1. Containing both gtex whole blood samples and rare disease samples
#     2. Only rare disease samples
if false; then
python append_samples_to_junction_file.py $gtex_whole_blood_jxn_file $sample_fastq_pairing_file $rail_rna_dir $gtex_rare_combined_jxn_file $rare_only_jxn_file
fi

###############
# Run dirichlet multinomial outlier calling on our junction matrices
# Junction matrices from combined whole blood gtex and rare-cohort (cgs) 
max_dm_junctions="20"
jxn_filter_method="ignore_genes"
node_number="0"
total_nodes="10"
covariate_regression_method="none"
min_reads_per_individual="5"
min_individuals_per_gene="50"

if false; then 
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    nohup python3 call_outlier_dm.py $gtex_rare_combined_jxn_file $outlier_calling_dir"gtex_rare_combined_" $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $min_reads_per_individual $min_individuals_per_gene &> log$node_number &
done
fi

###############
# Run dirichlet multinomial outlier calling on our junction matrices
# Junction matrices from rare-cohort (cgs) only
max_dm_junctions="20"
jxn_filter_method="ignore_genes"
node_number="0"
total_nodes="10"
covariate_regression_method="none"
min_reads_per_individual="5"
min_individuals_per_gene="30"

if false; then
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    nohup python3 call_outlier_dm.py $rare_only_jxn_file $outlier_calling_dir"rare_only_" $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $min_reads_per_individual $min_individuals_per_gene &> log_rare_only_$node_number &
done
fi
