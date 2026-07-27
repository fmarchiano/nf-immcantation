include { CHANGEO_ASSIGNGENES      } from '../../modules/local/changeo/assigngenes/main'
include { PRE_IGBLAST_ASSIGNGENES } from '../../modules/local/pre_igblast/assigngenes/main'
include { CHANGEO_MAKEDB          } from '../../modules/local/changeo/makedb/main'
include { CHANGEO_PARSEDB         } from '../../modules/local/changeo/parsedb/main'

workflow VDJ_ASSIGNMENT {
    take:
    ch_fasta      // [meta, *_atleast-2.fasta]
    ch_igblast    // value channel: igblast/ ref dir (staged per-task)
    ch_germlines  // value channel: germlines/ ref dir (staged per-task)

    main:
    // 1. IgBLAST V(D)J alignment
    if (params.pre_igblast) {
        PRE_IGBLAST_ASSIGNGENES(ch_fasta, ch_igblast, ch_germlines)
        ch_blast = PRE_IGBLAST_ASSIGNGENES.out.blast
        ch_fasta_out = PRE_IGBLAST_ASSIGNGENES.out.fasta
    } else {
        CHANGEO_ASSIGNGENES(ch_fasta, ch_igblast)
        ch_blast = CHANGEO_ASSIGNGENES.out.blast
        ch_fasta_out = CHANGEO_ASSIGNGENES.out.fasta
    }

    // 2. Parse IgBLAST output into AIRR-compliant TSV
    CHANGEO_MAKEDB(
        ch_fasta_out,
        ch_blast,
        ch_germlines
    )

    // 3. Filter productive + strip allele into v_call_gene / j_call_gene
    //    (Briney clonotype definition uses gene-level V/J)
    CHANGEO_PARSEDB(CHANGEO_MAKEDB.out.tab)

    emit:
    airr_tab = CHANGEO_PARSEDB.out.tab
}
