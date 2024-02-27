qiime tools import --type 'FeatureData[Sequence]' \
    --input-path <(gunzip -c ../output/localtest/vsearch/all_reads_full_its/all_reads_full_its.centroids.fasta.gz) \
    --output-path centroids.qza

qiime tools import --type 'FeatureData[Sequence]' \
    --input-path <(gunzip -c ../data/db/unite-qiime-2020725/sh_refs_qiime_ver9_99_all_25.07.2023.fasta) \
    --output-path UNITE.qza

qiime tools import --type 'FeatureData[Taxonomy]' \
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
    --parallel


#qiime tools extract --input-path classifications.qza  --output-path extracted


biom convert \
    -i ../output/localtest/vsearch-otu-map/all_reads_full_its/all_reads_full_its.otus.tsv \
    -o otus.biom \
    --to-hdf5
    
qiime tools import --type 'FeatureTable[Frequency]' \                                                                                               ✔  4s  qiime2-amplicon-2024.2 Py
  --input-path otus.biom \
  --output-path otus.gza