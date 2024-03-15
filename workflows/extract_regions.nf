include { 
    ITSXPRESS as ITSXPRESS_ITS1; 
    ITSXPRESS as ITSXPRESS_ITS2;
    ITSXPRESS as ITSXPRESS_FULL } from '../modules/local/itsxpress'
include { SEQKIT_FQ2FA } from '../modules/local/seqkit/fq2fa'
include { ITSX as ITSX_LSU } from '../modules/local/itsx'

workflow ExtractRegions {
    take:
        reads // reads to classify in fastq format

    main:
        its1 = ITSXPRESS_ITS1(reads, "ITS1")
        its2 = ITSXPRESS_ITS2(reads, "ITS2")
        full_its = ITSXPRESS_FULL(reads, "FULL_ITS")
        itsx = SEQKIT_FQ2FA(reads) | ITSX_LSU
        
        lsu_mapped = RecoverRegionsFromFastq(reads.join(itsx.positions))
    
    emit:
        its1 = its1.reads.map(addRegionToMetadata("ITS1"))
        its2 = its2.reads.map(addRegionToMetadata("ITS2"))
        full_its = full_its.reads.map(addRegionToMetadata("FULL_ITS"))
        lsu = lsu_mapped.lsu.map(addRegionToMetadata("LSU"))
}

def addRegionToMetadata(regionId) {
    addToMetadata([region: "$regionId"])
}

def addToMetadata(toAdd) {
    { meta, reads -> [meta + toAdd, reads] }
}

process RecoverRegionsFromFastq {
    tag "$meta.id - LSU"
    container "biocontainers/biopython:1.81"

    input:
        tuple val(meta), path(fastq), path(positions)
    
    output:
        tuple val(meta), path("*LSU.fastq.gz"), emit: lsu

    script:
    def outfile = "${fastq.baseName}".replaceAll(/.fastq/, '')
    """
    extract_region_by_pos.py --fastq $fastq \\
        --positions $positions \\
        -o ${outfile}.LSU.fastq.gz
    """

}