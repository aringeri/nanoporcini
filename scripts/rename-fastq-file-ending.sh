for a in $(ls -1 .); do mv $a "${a%.fq.gz}.fastq.gz"; done
