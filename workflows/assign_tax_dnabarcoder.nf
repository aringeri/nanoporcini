include { UnGzip } from './classify'
include { SEQKIT_FQ2FA } from '../modules/local/seqkit/fq2fa'


workflow assignTaxDnabarcoder {
    take:
        seqs_to_assign // Channel<Map<(Region, Scenario, Id)>, Fastq.GZ> - representative sequences from each cluster

    main:
        ref_db = Channel.fromPath(params.taxonomic_assignment.dnabarcoder.ref_db, checkIfExists: true)
        ref_classifications = Channel.fromPath(params.taxonomic_assignment.dnabarcoder.ref_classifications, checkIfExists: true)
        cutoffs = Channel.fromPath(params.taxonomic_assignment.dnabarcoder.cutoffs, checkIfExists: true)

        seqs_fasta = seqs_to_assign | SEQKIT_FQ2FA | UnGzip
        matches = dnabarcoderSearch(seqs_fasta.combine(ref_db)).best_matches

        dnabarcoderClassify(matches.combine(ref_classifications).combine(cutoffs))
}

process dnabarcoderSearch {
    tag "${meta.scenario.count}/${meta.scenario.rep}"
    container "ghcr.io/aringeri/dnabarcoder:3b21bc924486f3ccdba5d7ebe9b7cc5f7d069bc1"
    label "mega_mem"
    label "med_cpu"

    input:
        tuple val(meta), path(fasta), path(ref_db)
    output:
        tuple val(meta), path("*.bestmatch"), emit: best_matches
        tuple val(meta), path("*.log"),       emit: log

    script:
    """
    dnabarcoder.py search \\
        --ncpus ${task.cpus} \\
        -i $fasta \\
        -r $ref_db \\
        -o . > dnabarcoder_search.log
    """
}

process dnabarcoderClassify {
    tag "${meta.scenario.count}/${meta.scenario.rep}"
    container "ghcr.io/aringeri/dnabarcoder:3b21bc924486f3ccdba5d7ebe9b7cc5f7d069bc1"
    label "large_mem"
    label "large_cpu"

    input:
        tuple val(meta), path(best_matches), path(ref_classifications), path(cutoffs)
    output:
        tuple val(meta), path("*.classified"),     emit: classified
        tuple val(meta), path("*.classification"), emit: classification
        tuple val(meta), path("*.krona.report"),   emit: krona_report
        tuple val(meta), path("*.log"),            emit: log

    script:
    """
    dnabarcoder.py classify \\
        -i $best_matches \\
        -c $ref_classifications \\
        -cutoffs $cutoffs \\
        -o . > dnabarcoder_classify.log
    """
}