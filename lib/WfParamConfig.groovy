//file:noinspection GroovyAssignabilityCheck
class WfParamConfig {

    private static final List<WfParam> WORKFLOW_CONFIG = [
        new WfParam(path: "primers.fwd", required: true, description: 'The forward primer sequence in 5\' to 3\' orientation.'),
        new WfParam(path: "primers.rev_rc", required: true, description: 'The reverse primer sequence (which has been reverse complemented)'),
        new WfParam(path: "input", required: true, description: 'A set of fastq files. One file per sample. File names will be used to identify samples. Uses glob syntax (`*`) to select multiple files.'),
        new WfParam(path: "outdir", required: true, description: 'The output directory. Will be created if it does not exist.'),
        new WfParam(path: "chimera_filtering.ref_db", required: true, description: 'The reference database (fasta file) to be used for chimera detection.'),

        new WfParam(path: 'trim_adapters', required: false, defaultValue: false),
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
        new WfParam(path: "taxonomic_assignment.dnabarcoder.ref_db", required: true, description: 'The reference database (fasta file) to be used for taxonomic assignments.'),
        new WfParam(path: "taxonomic_assignment.dnabarcoder.ref_classifications", required: true, description: 'The reference database classifications file required by dnabarcoder.'),
        new WfParam(path: "taxonomic_assignment.dnabarcoder.cutoffs", required: true, description: 'The reference database cutoffs (json) required by dnabarcoder'),

        new WfParam(path: "consensus.methods", required: false, defaultValue: []),
        new WfParam(path: "consensus.num_polishing_reads", required: false, defaultValue: 200),

        new WfParam(path: "qc_quality_profile", required: false, defaultValue: false),
        new WfParam(path: "qc_plot_sample_level", required: false, defaultValue: false),
        new WfParam(path: "sample_barcode_in_file_name", required: false, defaultValue: false),
        new WfParam(path: "exclude", required: false, defaultValue: []),
        new WfParam(path: "help", required: false, description: 'Print this help message'),
    ]

    private static String SPACE = "   "

    static String helpString() {
        def config = WORKFLOW_CONFIG

        return """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
▗▖  ▗▖ ▗▄▖ ▗▖  ▗▖ ▗▄▖ ▗▄▄▖  ▗▄▖ ▗▄▄▖  ▗▄▄▖▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖
▐▛▚▖▐▌▐▌ ▐▌▐▛▚▖▐▌▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▌     █  ▐▛▚▖▐▌  █  
▐▌ ▝▜▌▐▛▀▜▌▐▌ ▝▜▌▐▌ ▐▌▐▛▀▘ ▐▌ ▐▌▐▛▀▚▖▐▌     █  ▐▌ ▝▜▌  █  
▐▌  ▐▌▐▌ ▐▌▐▌  ▐▌▝▚▄▞▘▐▌   ▝▚▄▞▘▐▌ ▐▌▝▚▄▄▖▗▄█▄▖▐▌  ▐▌▗▄█▄▖
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage:
    nextflow run https://github.com/aringeri/nanoporcini -params-file <parameter-file-yaml> [-c <env-config>]

    nextflow run https://github.com/aringeri/nanoporcini 
        --input <input> --outdir <outdir>
        --primers.fwd <fwd-primer> --primers.rev_rc <rev-primer>
        --chimera_filtering.ref_db <ref-db>
        --taxonomic_assignment.dnabarcoder.ref_db <ref-db>
        --taxonomic_assignment.dnabarcoder.ref_classifications <ref-classifications>
        --taxonomic_assignment.dnabarcoder.cutoffs <ref-cutoffs>
        [options]
        [-params-file <parameter-file-yaml>]
        [-c <env-config>]

Options:
    Required: 
${describeOptions(SPACE + SPACE, config.findAll { it.required })}
    
    Optional:
${describeOptions(SPACE + SPACE, config.findAll { !it.required })}

    Nextflow specific:
        -params-file
            Load script parameters from a JSON/YAML file
        -c, -config
            Add the specified file to configuration set
        -p, -profile
            Choose a configuration profile: docker (default), singularity
        -h, -help
            Print the nextflow command usage
"""
    }

    static String describeOptions(String tab, params) {
        return params
            .collect { "${tab}--${it.path}\n${SPACE}${tab}${it.description}"}
            .join("\n")
    }

    static List<ParamValidationError> validateParams(params) {
        def config = WORKFLOW_CONFIG

        return config.collectMany{ param ->
            validateParam(param, params)
        }
    }

    private static List<ParamValidationError> validateParam(WfParam param, params) {
        def actual = lookupParam(params, param)
        if (param.required && actual == null) {
            return [new ParamValidationError(message: "Required parameter (${param.path}) is missing.")]
        }

        if (actual == null && param.defaultValue != null) {
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