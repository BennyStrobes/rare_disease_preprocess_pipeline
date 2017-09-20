
##########################################################################
# Input Files
##########################################################################
# Directory containing fastq file for each sample
fastq_directory="/users/lfresard/CGS_fastq_temp/"

# Directory containing hg19 genome (fasta) as well as bowtie indices
genome_directory="/mnt/lab_data/montgomery/shared/genomes/"









##########################################################################
# Output Directories (Assumes these directories exist before script starts)
##########################################################################

# Directory used to store information on each of our samples
# Also contains one manifest file (to be used for rail-rna) for each sample. 
# Manifest files named $sample_name"_rail_manifest.txt"
sample_info_dir="/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/sample_info/"






##########################################################################
# Output Files (Do not exist before program starts). 
##########################################################################

# File that contains mapping from sample name (column 1) to complete path of fastq_left (column 2) and fastq_right (column 3)
# Created by `python learn_sample_fastq_pairing.py`
sample_fastq_pairing_file=$sample_info_dir"sample_fastq_pairing.txt"



##########################################################################
# Scripts/Functions/Analyses
##########################################################################

# Create file that contains mapping from sample name (column 1) to complete path of fastq_left (column 2) and fastq_right (column 3)
# Also creates one manifest file (to be used for rail-rna) for each sample. 
# Manifest files named $sample_info_dir$sample_name"_rail_manifest.txt"
if false; then
python learn_sample_fastq_pairing.py $fastq_directory $sample_fastq_pairing_file $sample_info_dir
fi