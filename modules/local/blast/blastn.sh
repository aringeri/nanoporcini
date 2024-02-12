blastn -db db/blast/utax_reference_dataset_all_25.07.2023.blast \
    -query vsearch/its1/all-reads/all-reads.consensus.fasta \
    -outfmt 7 -out out.tsv