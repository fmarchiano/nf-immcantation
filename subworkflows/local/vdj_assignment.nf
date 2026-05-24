include { CHANGEO_ASSIGNGENES } from '../../modules/local/changeo/assigngenes/main'
include { CHANGEO_MAKEDB      } from '../../modules/local/changeo/makedb/main'
include { CHANGEO_PARSEDB     } from '../../modules/local/changeo/parsedb/main'

workflow VDJ_ASSIGNMENT {
    take:
    ch_fasta   // [meta, *_atleast-2.fasta]

    main:
    // 1. IgBLAST V(D)J alignment — uses $IGDATA inside container
    CHANGEO_ASSIGNGENES(ch_fasta)

    // 2. Parse IgBLAST output into AIRR-compliant TSV
    CHANGEO_MAKEDB(
        CHANGEO_ASSIGNGENES.out.fasta,
        CHANGEO_ASSIGNGENES.out.blast
    )

    // 3. Filter productive + strip allele into v_call_gene / j_call_gene
    //    (Briney clonotype definition uses gene-level V/J)
    CHANGEO_PARSEDB(CHANGEO_MAKEDB.out.tab)

    emit:
    airr_tab = CHANGEO_PARSEDB.out.tab
}
