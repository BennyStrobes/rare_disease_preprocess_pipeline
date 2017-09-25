

bowtie_prefix="$1"
bowtie2_prefix="$2"
manifest_file="$3"
sample_output_directory="$4"
rail_rna_dir="$5"

rail-rna go local -p 8 -d jx -x $bowtie_prefix $bowtie2_prefix -m $manifest_file -o $sample_output_directory --log $rail_rna_dir --force

