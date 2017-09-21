
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



# Run rail-rna for each sample
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

# TEMPORARY
while read sample_id fq_file1 fq_file2
do
    echo $sample_id
    manifest_file=$sample_info_dir$sample_id"_rail_manifest.txt"
    sample_output_directory=$rail_rna_dir$sample_id"_data/"
    sample_log_dir=$rail_rna_log_dir$sample_id"_log/"
    mkdir $sample_log_dir
    # Only perform jxn 1st pass quantification (as done by Recount/Snaptron)
    nohup sh rail_wrapper.sh $bowtie_prefix $bowtie2_prefix $manifest_file $sample_output_directory $sample_log_dir &>/dev/null &
done<"/srv/scratch/bstrobe1/rare_variant/rare_disease/preprocess/sample_info/sample_fastq_pairing_part1.txt"
