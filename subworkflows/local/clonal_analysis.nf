include { CHANGEO_DEFINECLONES       } from '../../modules/local/changeo/defineclones/main'
include { SCOPER_HIERARCHICALCLONES  } from '../../modules/local/scoper/hierarchicalclones/main'
include { SHAZAM_DISTTONEAREST       } from '../../modules/local/shazam/disttonearest/main'

workflow CLONAL_ANALYSIS {
    take:
    ch_airr   // [meta, *_parse-pass.tsv]

    main:
    // Group samples by cloneby (e.g. subject_id) so clones are defined across
    // all biological/technical replicates from the same donor
    ch_airr
        .map { meta, tab ->
            def group_key = meta[params.cloneby]
            def group_meta = [
                id         : group_key,
                subject_id : meta.subject_id,
                species    : meta.species,
                locus      : meta.locus
            ]
            [ group_meta, tab ]
        }
        .groupTuple(by: [0])
        .set { ch_grouped }

    // Method switch: 'exact' (Briney parity) or 'hierarchical' (SCOPer)
    ch_threshold = Channel.empty()
    if (params.cloning_method == 'exact') {
        CHANGEO_DEFINECLONES(ch_grouped)
        ch_cloned = CHANGEO_DEFINECLONES.out.tab
    } else if (params.cloning_method == 'hierarchical') {
        SCOPER_HIERARCHICALCLONES(ch_grouped)
        ch_cloned = SCOPER_HIERARCHICALCLONES.out.tab

        // Optional drift check: compute the data-driven valley per patient and
        // warn if it has drifted from the fixed clonal_threshold. Informational
        // only — clustering above always uses the fixed threshold.
        if (params.clone_threshold_mode == 'disttonearest') {
            SHAZAM_DISTTONEAREST(ch_grouped)
            ch_threshold = SHAZAM_DISTTONEAREST.out.threshold
        }
    } else {
        error "Invalid params.cloning_method: '${params.cloning_method}'. Must be 'exact' or 'hierarchical'."
    }

    emit:
    cloned_tab = ch_cloned
    threshold  = ch_threshold
}
