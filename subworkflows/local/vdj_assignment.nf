include { CHANGEO_ASSIGNGENES } from '../../modules/local/changeo/assigngenes/main'
include { CHANGEO_MAKEDB      } from '../../modules/local/changeo/makedb/main'

workflow VDJ_ASSIGNMENT {
    take:
    ch_fasta   // [meta, *_atleast-2.fasta]

    main:
    // 1. IgBLAST V(D)J alignment — uses $IGDATA inside container
    CHANGEO_ASSIGNGENES(ch_fasta)

    // 2. Parse IgBLAST output into AIRR-compliant TSV
    //    germlines from /usr/local/share/germlines/imgt/{species}/vdj/ inside container
    CHANGEO_MAKEDB(
        CHANGEO_ASSIGNGENES.out.fasta,
        CHANGEO_ASSIGNGENES.out.blast
    )

    emit:
    airr_tab = CHANGEO_MAKEDB.out.tab
}
