variant_bed_file="$1"
cluster_file="$2"
output_dir="$3"
distance="$4"
vcf_individuals="$5"
outlier_file="$6"
version="$7"


variant_cluster_bed_file=$output_dir"variant_cluster_bed_"$distance"_"$version".txt"
# First variants to only those near splice junctions
python map_rare_variants_to_clusters.py $variant_bed_file $cluster_file $variant_cluster_bed_file $distance


# Compute enrichments at various thresholds
thresh=".01"
python variant_enrichment_quantification.py $variant_cluster_bed_file $output_dir"variant_enrichment_"$distance"_"$version"_"$thresh".txt" $outlier_file $vcf_individuals $thresh

thresh=".001"
python variant_enrichment_quantification.py $variant_cluster_bed_file $output_dir"variant_enrichment_"$distance"_"$version"_"$thresh".txt" $outlier_file $vcf_individuals $thresh

thresh=".0001"
python variant_enrichment_quantification.py $variant_cluster_bed_file $output_dir"variant_enrichment_"$distance"_"$version"_"$thresh".txt" $outlier_file $vcf_individuals $thresh

thresh=".00001"
python variant_enrichment_quantification.py $variant_cluster_bed_file $output_dir"variant_enrichment_"$distance"_"$version"_"$thresh".txt" $outlier_file $vcf_individuals $thresh


cat $output_dir"variant_enrichment_"$distance"_"$version"_"*"txt" > $output_dir"variant_enrichment_"$distance"_"$version"_cat.txt"





# Visualize the enrichments
python or_enrichment_plotter.py $output_dir"variant_enrichment_"$distance"_"$version"_cat.txt" $output_dir"variant_enrichment_"$distance"_"$version"_boxplot.png"