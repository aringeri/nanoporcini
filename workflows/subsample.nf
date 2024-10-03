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
            if (params.subsample.type == 'uneven') {
                 sampling = channel.fromList(params.subsample.scenarios)
                    .flatMap { scenario ->
                        reps = params.subsample.repetitions
                        seeds = generate_seed_for_each_rep(scenario.seed, reps)
                        log.debug("Uneven subsample seeds for scenario (count=${scenario.count}): $seeds")
                        [  Collections.nCopies(reps, scenario), (1..reps).toList(), seeds ].transpose()
                    }
                    .combine(qcd_reads)
                    .map { scenario, rep, seed, meta, reads ->
                        sample_depth = scenario.count
                        if (meta.id in scenario.mixin) {
                            sample_depth = scenario.mixin[meta.id].prop * scenario.count
                        }
                        [meta + [scenario: [count: scenario.count, seed: seed, rep: rep]] + [ seqtk_sample_depth: sample_depth ], reads]
                    }
            } else if (params.subsample.type == 'even') {
                sampling = channel.fromList(params.subsample.scenarios)
                    .flatMap { scenario ->
                        reps = params.subsample.repetitions
                        seeds = generate_seed_for_each_rep(scenario.seed, reps)
                        log.debug("Even subsample seeds for scenario (count=${scenario.count}): $seeds")
                        // zip subsample count with each seed
                        return [ Collections.nCopies(reps, scenario.count), (1..reps).toList(), seeds ].transpose()
                    }
                    .combine(qcd_reads)
                    .map { count, rep, seed, meta, reads ->
                        [meta + [scenario: [count: count, seed: seed, rep: rep]] + [ seqtk_sample_depth: count ], reads]
                    }
            } else {
                error("Subsampling type (${params.subsample.type}) not recognized. Expecting one of: ('even', 'uneven') ")
            }

            subsampled = seqtk_sample(sampling).subsampled
        } else {
            log.info("subsample.enabled = false: Skipping subsampling")
            subsampled = qcd_reads.map{ meta, reads ->
                 [ meta + [scenario: [count: "ALL_READS", seed: "NA", rep: "1"]], reads ]
            }
        }
    emit:
        subsampled = subsampled
}