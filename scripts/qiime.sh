qiime tools import --type 'FeatureData[Sequence]' \
    --input-path <(gunzip -c ../output/localtest/vsearch/all_reads_full_its/all_reads_full_its.centroids.fasta.gz) \
    --output-path centroids.qza

qiime tools import --type 'FeatureData[Sequence]' \
    --input-path <(gunzip -c ../data/db/unite-qiime-2020725/sh_refs_qiime_ver9_99_all_25.07.2023.fasta) \
    --output-path UNITE.qza

qiime tools import --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path ../data/db/unite-qiime-2020725/sh_taxonomy_qiime_ver9_99_all_25.07.2023.txt \
    --output-path UNITE_TAX.qza

qiime feature-classifier makeblastdb --i-sequences UNITE.qza --o-database UNITE_DB.qza

qiime feature-classifier classify-consensus-blast \
    --i-query centroids.qza \
    --i-blastdb UNITE_DB.qza \
    --i-reference-taxonomy UNITE_TAX.qza \
    --o-classification classifications.qza \
    --o-search-results search-results.qza \
    --verbose \
    --p-num-threads 4 \
    --parallel


#qiime tools extract --input-path classifications.qza  --output-path extracted


biom convert \
    -i ../output/localtest/vsearch-otu-map/all_reads_full_its/all_reads_full_its.otus.tsv \
    -o otus.biom \
    --to-hdf5
    
qiime tools import --type 'FeatureTable[Frequency]' \
  --input-path otus.biom \
  --output-path otus.gza


awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
   NR>1 {split($1, arr, "|"); print(arr[3]"_"arr[2]"_"arr[4]"\t"arr[5])}
  ' data/db/UNITE-full-all-10.15156-BIO-2938070-20230725/sh_general_release_dynamic_s_all_25.07.2023.fasta

awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
   NR>1 {split($1, arr, "|"); print(">"arr[3]"_"arr[2]"_"arr[4]"\n"$2)}
  ' data/db/UNITE-full-all-10.15156-BIO-2938070-20230725/sh_general_release_dynamic_s_all_25.07.2023.fasta  \
  > qiime/sh_general_release_dynamic_s_all_25.07.2023_qiime.fasta