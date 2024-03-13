vsearch --fastx_filter scratch/Sample1.BC89.fastq \
    --eeout \
    --fastq_qmax 93 \
    --fastq_ascii 33 \
    --fastq_maxee 1 \
    --fastqout scratch/Sample1-BC89-filtered.fastq


vsearch \
    --fastq_eestats scratch/Sample1.BC89.fastq \
    --fastq_qmax 93 \
    --fastq_ascii 33 \
    --output scratch/eestats.tsv 

vsearch \
    --fastq_eestats2 scratch/Sample1.BC89.fastq \
    --fastq_qmax 93 \
    --fastq_ascii 33 \
    --output scratch/eestats2.tsv 