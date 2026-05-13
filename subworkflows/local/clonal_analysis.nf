include { CHANGEO_DEFINECLONES } from '../../modules/local/changeo/defineclones/main'

workflow CLONAL_ANALYSIS {
    take:
    ch_airr   // [meta, *db-pass.tsv]

    main:
    // Group samples by cloneby field (e.g. subject_id) before running DefineClones
    // so clones are defined across all samples from the same donor
    ch_airr
        .map { meta, tab ->
            def group_key = meta."${params.cloneby}"
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

    // DefineClones with fixed threshold — no Shazam auto-detection
    CHANGEO_DEFINECLONES(ch_grouped)

    emit:
    cloned_tab = CHANGEO_DEFINECLONES.out.tab
}
