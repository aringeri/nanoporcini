import java.util.Random
import java.util.stream.Collectors
import java.util.Collections

include { seqtk_sample } from '../modules/local/seqtk/seqtk'

def generate_seed_for_each_rep(int initialSeed, int repetitions) {
    return new Random(initialSeed).ints()
        .limit(repetitions)
        .boxed()
        .collect(Collectors.toList())
}

workflow subsample {
    take:
        qcd_reads // Channel<(Map, Fastq.GZ)> - One tuple per sample

    main:
        if (params.subsample.enabled) {
            out = channel.fromList(params.subsample.scenarios)
                .flatMap { scenario ->
                    reps = params.subsample.repetitions
                    seeds = generate_seed_for_each_rep(scenario.seed, reps)
                    log.debug("Subsample seeds for scenario (count=${scenario.count}): $seeds")
                    // zip subsample count with each seed
                    return [ Collections.nCopies(reps, scenario.count), (1..reps).toList(), seeds ].transpose()
                }
                .combine(qcd_reads)
                .map { count, rep, seed, meta, reads ->  [meta + [scenario: [count: count, seed: seed, rep: rep]], reads] }
                //.view { "subsample scenario: $it"}
                | seqtk_sample
            subsampled = out.subsampled
        } else {
            log.info("subsample.enabled = false: Skipping subsampling")
            subsampled = qcd_reads
        }
    emit:
        subsampled = subsampled
}