
cgs_vcf_input_dir="$1"
variant_calling_dir="$2"
sample_info_pairing_file="$3"
genomes_1000_dir="$4"
gencode_gene_annotation_file="$5"

if false; then
# Merge vcf files from each sample (1 core. About x min)
vcf-merge $(ls -1 $cgs_vcf_input_dir*.vcf.gz | perl -pe 's/\n/ /g') > $variant_calling_dir"merge.vcf"


# Change from institution id to sample id
# About a minute
python convert_merged_vcf_from_institution_id_to_sample_id.py $variant_calling_dir"merge.vcf" $variant_calling_dir"merge_sample_names.vcf" $sample_info_pairing_file


# Zip and tabix up our vcf files
# A few minutes
bgzip -c $variant_calling_dir"merge.vcf" > $variant_calling_dir"merge.vcf.gz"
tabix -p vcf $variant_calling_dir"merge.vcf.gz"

bgzip -c $variant_calling_dir"merge_sample_names.vcf" > $variant_calling_dir"merge_sample_names.vcf.gz"
tabix -p vcf $variant_calling_dir"merge_sample_names.vcf.gz"

# Filter out indels and remove flagged sites
# Filter out chr X, Y, and MT
# A few minutes
vcftools --gzvcf $variant_calling_dir"merge_sample_names.vcf.gz"  $variant_calling_dir"merge_sample_names_SNPs" --remove-filtered-all --remove-indels --not-chr X --not-chr Y --not-chr MT --recode &
vcftools --gzvcf $variant_calling_dir"merge_sample_names.vcf.gz" --out $variant_calling_dir"merge_sample_names_SNPs" --remove-filtered-all --remove-indels --not-chr X --not-chr Y --not-chr MT --freq &
vcftools --gzvcf $variant_calling_dir"merge_sample_names.vcf.gz" --out $variant_calling_dir"merge_sample_names_SNPs" --remove-filtered-all --remove-indels --not-chr X --not-chr Y --not-chr MT --extract-FORMAT-info GT &



## Extract AF of variants according to 1000 genomes super-populations
sh extract.1kg.AF.sh $genomes_1000_dir $variant_calling_dir


# output of extract.1kg.AF.sh
genomes_1000_af_file=$variant_calling_dir"SNPs.1kg.AF.bed.gz"

# output of vcftools filtering
filtered_vcf=$variant_calling_dir"merge_sample_names_SNPs.recode.vcf"

# Make a bed file for each sample.
# bed file contains all variants and their AF
python create_sample_specific_variant_bed_files.py $genomes_1000_af_file $filtered_vcf $variant_calling_dir

fi

# Map variants to nearby genes
python map_variants_to_genes.py $variant_calling_dir $sample_info_pairing_file $variant_calling_dir"all_variants_gene_mapping.txt" $gencode_gene_annotation_file



















################################################
# Old code (ignore for now)
#################################################


if false; then

vcftools --gzvcf $variant_calling_dir"merge_sample_names.vcf.gz" --out $variant_calling_dir"merge_sample_names_SNPs" --remove-filtered-all --remove-indels --freq &
vcftools --gzvcf $variant_calling_dir"merge_sample_names.vcf.gz" --out $variant_calling_dir"merge_sample_names_SNPs" --remove-filtered-all --remove-indels --extract-FORMAT-info GT &

#dir=${VCFPREFIX%/*} # keep same directory as vcf prefix
#outdir=${dir}/individuals # some future scripts (e.g. CADD one) assume this directory, so do not change this without changing downstream

freq=$variant_calling_dir"merge_sample_names_SNPs.frq"
gt=$variant_calling_dir"merge_sample_names_SNPs.GT.FORMAT"
# output
af=$variant_calling_dir/AF_SNPs.bed

##################
# SOME FUNCTIONS #
##################
# skip first line of file
alias skip_header="tail -n+2"

# to remove missing genotypes
alias filter_genotypes="grep -v '\./\.'"

# to turn genotypes into 0,1,2 (the number of non-reference alleles)
# 0/0 -> 0; 0/*,*/0 -> 1; */* -> 2 (where * is any number >0)
alias clean_genotypes="sed 's%0/0%0%' | sed 's%[1-9]/0%1%' | sed 's%0/[1-9]%1%' | sed 's%[1-9]/[1-9]%2%'"



# PROCESS ALLELE FREQUENCY FILE
# make bed file with minor allele frequency and whether the minor allele homozygote is 0 or 2
# for SNPs: columns 6 onwards are the alleles in order presented in the vcf with "." for columns up to 10 (ATCG*) if no extra alleles
# if there are multiple non reference alleles, the non-ref AF is set to 1-ref AF
# skips chrom X, Y, and MT

cat $freq | tr ":" "\t" | skip_header | \
    awk 'BEGIN{OFS="\t"}{
               af=$6;hz=0; 
               if(af=="-nan" || $1=="X" || $1=="Y" || $1=="MT"){next}; 
               if(af>0.5){af=1-$6;hz=2}; 
               printf "%i\t%i\t%i\t%f\t%i",$1,$2-1,$2,af,hz; 
               for(i=5;i<=NF-1;i=i+2){printf "\t%s",$i};
               for(j=NF+1;j<=13;j=j+2){printf "\t."};
               printf "\n"}' | \
    sort -k1,1 -k2,2n > $af

fi
