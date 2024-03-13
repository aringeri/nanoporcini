include { 
    ITSXPRESS as ITSXPRESS_ITS1; 
    ITSXPRESS as ITSXPRESS_ITS2;
    ITSXPRESS as ITSXPRESS_FULL } from '../modules/local/itsxpress'
include { SEQKIT_FQ2FA } from '../modules/local/seqkit/fq2fa'
include { ITSX } from '../modules/local/itsx'

workflow ExtractRegions {
    take:
        reads // reads to classify in fastq format

    main:
        its1 = ITSXPRESS_ITS1(reads, Regions.ITS1)
        its2 = ITSXPRESS_ITS2(reads, Regions.ITS2)
        full_its = ITSXPRESS_FULL(reads, Regions.FULL_ITS)
        itsx = SEQKIT_FQ2FA(reads) | ITSX
        
        lsu_mapped = RecoverRegionsFromFastq(reads.join(itsx.positions))
    
    emit:
        its1 = its1.reads
        its2 = its2.reads
        full_its = full_its.reads
        lsu = lsu_mapped.lsu
}

process RecoverRegionsFromFastq {
    container "biocontainers/biopython:1.81"

    input:
        tuple val(meta), path(fastq), path(positions)
    
    output:
        tuple val(meta), path("*LSU.fastq.gz"), emit: lsu

    script:
    """
    extract_region_by_pos.py --fastq $fastq \\
        --positions $positions \\
        -o ${fastq.baseName}.LSU.fastq.gz
    """

}