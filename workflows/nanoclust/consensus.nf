process readCorrection {
     memory { 7.GB }
     container "quay.io/biocontainers/canu:2.2--ha47f30e_0"

     input:
        tuple val(meta), path(reads)

     output:
        tuple val(meta), path('corrected_reads.correctedReads.fasta'), emit: corrected_reads

     script:
     genomeSize=1000
     minReadLength=400
     stopOnLowCoverage=1
     minInputCoverage=2
     """
     canu -correct -p corrected_reads -nanopore-raw \\
        $reads genomeSize=${genomeSize} \\
        stopOnLowCoverage=${stopOnLowCoverage} minInputCoverage=${minInputCoverage} \\
        minReadLength=${minReadLength} minOverlapLength=200

     gunzip corrected_reads.correctedReads.fasta.gz
     """
 }

//  process draft_selection {
//      publishDir "${params.outdir}/${barcode}/cluster${cluster_id}", mode: 'copy', pattern: 'draft_read.fasta'
//      errorStrategy 'retry'
//
//      input:
//      tuple val(barcode), val(cluster_id), file(cluster_log), file(reads) from corrected_reads
//
//      output:
//      tuple val(barcode), val(cluster_id), file('*_draft.log'), file('draft_read.fasta'), file(reads) into draft
//
//      script:
//      """
//      split -l 2 $reads split_reads
//      find split_reads* > read_list.txt
//
//      fastANI --ql read_list.txt --rl read_list.txt -o fastani_output.ani -t 48 -k 16 --fragLen 160
//
//      DRAFT=\$(awk 'NR>1{name[\$1] = \$1; arr[\$1] += \$3; count[\$1] += 1}  END{for (a in arr) {print arr[a] / count[a], name[a] }}' fastani_output.ani | sort -rg | cut -d " " -f2 | head -n1)
//      cat \$DRAFT > draft_read.fasta
//      ID=\$(head -n1 draft_read.fasta | sed 's/>//g')
//      cat $cluster_log > ${cluster_id}_draft.log
//      echo -n \$ID >> ${cluster_id}_draft.log
//     """
//  }
//
//  process racon_pass {
//      memory { 7.GB * task.attempt }
//      time { 1.hour * task.attempt }
//      errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
//      maxRetries 3
//
//      input:
//      tuple val(barcode), val(cluster_id), file(cluster_log), file(draft_read), file(corrected_reads) from draft
//
//      output:
//      tuple val(barcode), val(cluster_id), file(cluster_log), file('racon_consensus.fasta'), file(corrected_reads), env(success) into racon_output
//
//      script:
//      """
//      success=1
//      minimap2 -ax map-ont --no-long-join -r100 -a $draft_read $corrected_reads -o aligned.sam
//      if racon --quality-threshold=9 -w 250 $corrected_reads aligned.sam $draft_read > racon_consensus.fasta ; then
//         success=1
//      else
//         success=0
//         cat $draft_read > racon_consensus.fasta
//      fi
//
//      """
//  }
//
//  process medaka_pass {
//      memory { 7.GB * task.attempt }
//      time { 1.hour * task.attempt }
//      errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
//      maxRetries 3
//
//      publishDir "${params.outdir}/${barcode}/cluster${cluster_id}", mode: 'copy', pattern: 'consensus_medaka.fasta/consensus.fasta'
//
//      input:
//      tuple val(barcode), val(cluster_id), file(cluster_log), file(draft), file(corrected_reads), val(success) from racon_output
//
//      output:
//      tuple val(barcode), val(cluster_id), file(cluster_log), file('consensus_medaka.fasta/consensus.fasta') into final_consensus
//
//      script:
//      if(success == "0"){
//         log.warn """Sample $barcode : Racon correction for cluster $cluster_id failed due to not enough overlaps. Taking draft read as consensus"""
//         racon_warnings.add("""Sample $barcode : Racon correction for cluster $cluster_id failed due to not enough overlaps. Taking draft read as consensus""")
//      }
//      """
//      if medaka_consensus -i $corrected_reads -d $draft -o consensus_medaka.fasta -t 4 -m r941_min_high_g303 ; then
//         echo "Command succeeded"
//      else
//         cat $draft > consensus_medaka.fasta
//      fi
//      """
//
//  }