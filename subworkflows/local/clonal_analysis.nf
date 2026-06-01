include { SCOPER_HIERARCHICALCLONES } from '../../modules/local/scoper/hierarchicalclones/main'
include { SHAZAM_DISTTONEAREST      } from '../../modules/local/shazam/disttonearest/main'

workflow CLONAL_ANALYSIS {
    take:
    ch_airr   // [meta, *db-pass.tsv]

    main:
    // Group samples by cloneby field (e.g. subject_id) before running DefineClones
    // so clones are defined across all samples from the same donor
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

    SCOPER_HIERARCHICALCLONES(ch_grouped)

    // Optional drift check: compute the data-driven valley per patient and
    // warn if it has drifted from the fixed clonal_threshold. Informational
    // only — clustering above always uses the fixed threshold.
    ch_threshold = Channel.empty()
    if (params.clone_threshold_mode == 'disttonearest') {
        SHAZAM_DISTTONEAREST(ch_grouped)
        ch_threshold = SHAZAM_DISTTONEAREST.out.threshold
    }

    emit:
    cloned_tab = SCOPER_HIERARCHICALCLONES.out.tab
    threshold  = ch_threshold
}
