

bowtie_prefix="$1"
bowtie2_prefix="$2"
manifest_file="$3"
sample_output_directory="$4"
rail_rna_dir="$5"
liftover_directory="$6"

rail-rna go local -p 2 -d jx -x $bowtie_prefix $bowtie2_prefix -m $manifest_file -o $sample_output_directory --log $rail_rna_dir --force


if false; then
rail_rna_output_dir=$sample_output_directory"cross_sample_results/"
python liftover_jxns_hg38_to_hg19.py $rail_rna_output_dir"first_pass_junctions.tsv.gz" $rail_rna_output_dir $liftover_directory
fi