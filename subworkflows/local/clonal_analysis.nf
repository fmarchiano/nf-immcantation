include { CHANGEO_DEFINECLONES       } from '../../modules/local/changeo/defineclones/main'
include { CHANGEO_CREATEGERMLINES    } from '../../modules/local/changeo/creategermlines/main'
include { SCOPER_HIERARCHICALCLONES  } from '../../modules/local/scoper/hierarchicalclones/main'
include { SCOPER_SPECTRALCLONES_FAST } from '../../modules/local/scoper/spectralclones_fast/main'
include { SHAZAM_DISTTONEAREST       } from '../../modules/local/shazam/disttonearest/main'
include { SHAZAM_CLONALITYMEASURES   } from '../../modules/local/shazam/clonalitymeasures/main'

workflow CLONAL_ANALYSIS {
    take:
    ch_airr       // [meta, *_parse-pass.tsv]
    ch_germlines  // path to germlines dir

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
    } else if (params.cloning_method == 'spectral_fast') {
        SCOPER_SPECTRALCLONES_FAST(ch_grouped)
        ch_cloned = SCOPER_SPECTRALCLONES_FAST.out.tab
    } else {
        error "Invalid params.cloning_method: '${params.cloning_method}'. Must be 'exact', 'hierarchical', or 'spectral_fast'."
    }

    // Reconstruct germline sequences using clonal consensus
    CHANGEO_CREATEGERMLINES(ch_cloned, ch_germlines)

    // Per-sequence SHM + clone-level clonality measures
    SHAZAM_CLONALITYMEASURES(CHANGEO_CREATEGERMLINES.out.tab)

    emit:
    cloned_tab = ch_cloned
    shm_tab    = SHAZAM_CLONALITYMEASURES.out.tab
    measures   = SHAZAM_CLONALITYMEASURES.out.measures
    threshold  = ch_threshold
}
