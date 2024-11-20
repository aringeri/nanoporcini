//file:noinspection GroovyAssignabilityCheck
class WfParamConfig {

    private static List<WfParam> workflowConfiguration() {
        return [
            new WfParam(path: "primers.fwd", required: true),
            new WfParam(path: "primers.rev_rc", required: true),
            new WfParam(path: "input", required: true),
            new WfParam(path: "outdir", required: true),
            new WfParam(path: "chimera_filtering.ref_db", required: true),

            new WfParam(path: "qualityFiltering.FULL_ITS.minQualityPhred", required: false, defaultValue: 20),
            new WfParam(path: "qualityFiltering.FULL_ITS.minLength", required: false, defaultValue: 300),
            new WfParam(path: "qualityFiltering.FULL_ITS.maxLength", required: false, defaultValue: 6000),

            new WfParam(path: "extract.ITS1", required: false, defaultValue: false),
            new WfParam(path: "extract.ITS2", required: false, defaultValue: false),
            new WfParam(path: "extract.FULL_ITS", required: false, defaultValue: true),
            new WfParam(path: "extract.LSU", required: false, defaultValue: false),

            new WfParam(path: "subsample.enabled", required: false, defaultValue: false),
            new WfParam(path: "subsample.repetitions", required: false, defaultValue: 1),

            new WfParam(path: "cluster.methods", required: false, defaultValue: ['vsearch']),
            new WfParam(path: "cluster.shuffle.enabled", required: false, defaultValue: false),
            new WfParam(path: "cluster.shuffle.seed", required: false, defaultValue: 14),
            new WfParam(path: "cluster.vsearch.min_cluster_size", required: false, defaultValue: 1),
            new WfParam(path: "cluster.hdbscan.min_cluster_sizes", required: false, defaultValue: [5]),
            new WfParam(path: "cluster.gather_min_cluster_size_stats", required: false, defaultValue: false),

            new WfParam(path: "taxonomic_assignment.enabled", required: false, defaultValue: false), // TODO to remove
            new WfParam(path: "taxonomic_assignment.dnabarcoder.ref_db", required: true),
            new WfParam(path: "taxonomic_assignment.dnabarcoder.ref_classifications", required: true),
            new WfParam(path: "taxonomic_assignment.dnabarcoder.cutoffs", required: true),

            new WfParam(path: "consensus.methods", required: false, defaultValue: []),
            new WfParam(path: "consensus.num_polishing_reads", required: false, defaultValue: 200),

            new WfParam(path: "qc_quality_profile", required: false, defaultValue: false),
            new WfParam(path: "qc_plot_sample_level", required: false, defaultValue: false),
            new WfParam(path: "sample_barcode_in_file_name", required: false, defaultValue: false),
            new WfParam(path: "exclude", required: false, defaultValue: []),
        ]
    }

    static List<ParamValidationError> validateParams(params) {
        def config = workflowConfiguration()

        return config.collectMany{ param ->
            validateParam(param, params)
        }
    }

    private static List<ParamValidationError> validateParam(WfParam param, params) {
        def actual = lookupParam(params, param)
        if (param.required && actual == null) {
            return [new ParamValidationError(message: "Required parameter (${param.path}) is missing.")]
        }

        if (actual == null) {
            replaceWithDefault(params, param)
        }

        return []
    }

    private static Object lookupParam(params, param) {
        def p = params
        for (pathSegment in param.path.tokenize('.')) {
            if (p.containsKey(pathSegment)) {
                p = p[pathSegment]
            } else {
                return null
            }
        }
        return p
    }

    private static replaceWithDefault(params, WfParam param) {
        def p = params
        def path = param.path.tokenize('.')
        def firstSegments = path.size() > 1 ? path[0..-2] : []
        for (pathSegment in firstSegments) {
            if (!p.containsKey(pathSegment)) {
                p[pathSegment] = [:]
            }
            p = p[pathSegment]
        }

        def lastSegment = path[-1]
        p[lastSegment] = param.defaultValue
    }
}

class ParamValidationError {
    String message
}