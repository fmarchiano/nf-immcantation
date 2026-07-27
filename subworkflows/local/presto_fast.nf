include { FASTP                                                           } from '../../modules/nf-core/fastp/main'
include { GUNZIP                                                          } from '../../modules/local/presto/gunzip/main'
include { PRESTO_MASKPRIMERS_FAST as PRESTO_MASKPRIMERS_C                 } from '../../../implementations/presto/maskprimers-fast/nextflow/maskprimers_fast/main'
include { PRESTO_MASKPRIMERS_FAST as PRESTO_MASKPRIMERS_V                 } from '../../../implementations/presto/maskprimers-fast/nextflow/maskprimers_fast/main'
include { PRESTO_PAIRSEQ_FAST     as PRESTO_PAIRSEQ                       } from '../../../implementations/presto/pairseq-fast/nextflow/pairseq_fast/main'
include { PRESTO_PAIRSEQ_FAST     as PAIRSEQ_BARCODE                      } from '../../../implementations/presto/pairseq-fast/nextflow/pairseq_fast/main'
include { PRESTO_PAIRSEQ_FAST     as PAIRSEQ_CONSENSUS                    } from '../../../implementations/presto/pairseq-fast/nextflow/pairseq_fast/main'
include { PRESTO_BUILDCONSENSUS_FAST as BUILDCONSENSUS_V                  } from '../../../implementations/presto/buildconsensus-fast/nextflow/buildconsensus_fast/main'
include { PRESTO_BUILDCONSENSUS_FAST as BUILDCONSENSUS_C                  } from '../../../implementations/presto/buildconsensus-fast/nextflow/buildconsensus_fast/main'
include { PEAR                                                            } from '../../modules/local/pear/main'
include { PRESTO_ASSEMBLEPAIRS_FAST as PRESTO_ASSEMBLEPAIRS               } from '../../../implementations/presto/assemblepairs-fast/nextflow/assemblepairs_fast/main'
include { PRESTO_COLLAPSESEQ_FAST as PRESTO_COLLAPSESEQ                  } from '../../../implementations/presto/collapseseq-fast/nextflow/collapseseq_fast/main'
include { PRESTO_SPLITSEQ                                                 } from '../../modules/local/presto/splitseq/main'

workflow PRESTO_FAST {
    take:
    ch_reads
    ch_cprimers
    ch_vprimers

    main:
    ch_reads_fastp = ch_reads.map { meta, reads -> [ meta, reads, [] ] }
    FASTP(ch_reads_fastp, false, false, false)

    GUNZIP(FASTP.out.reads)

    ch_cread = GUNZIP.out.reads.map { meta, reads -> [ meta, reads[0] ] }
    ch_vread = GUNZIP.out.reads.map { meta, reads -> [ meta, reads[1] ] }

    PRESTO_MASKPRIMERS_C(ch_cread, ch_cprimers.collect(), 'C')
    PRESTO_MASKPRIMERS_V(ch_vread, ch_vprimers.collect(), 'V')

    ch_paired = PRESTO_MASKPRIMERS_V.out.reads
        .join(PRESTO_MASKPRIMERS_C.out.reads, by: [0])

    if (params.umi) {
        PAIRSEQ_BARCODE(ch_paired)

        ch_v_bc = PAIRSEQ_BARCODE.out.reads.map { meta, r1, r2 -> [ meta, r1 ] }
        ch_c_bc = PAIRSEQ_BARCODE.out.reads.map { meta, r1, r2 -> [ meta, r2 ] }
        BUILDCONSENSUS_V(ch_v_bc, 'V')
        BUILDCONSENSUS_C(ch_c_bc, 'C')

        ch_cons = BUILDCONSENSUS_V.out.reads
            .join(BUILDCONSENSUS_C.out.reads, by: [0])
        PAIRSEQ_CONSENSUS(ch_cons)

        ch_assemble = PAIRSEQ_CONSENSUS.out.reads
    } else {
        PRESTO_PAIRSEQ(ch_paired)
        ch_assemble = PRESTO_PAIRSEQ.out.reads
    }

    if (params.assembler == 'assemblepairs') {
        PRESTO_ASSEMBLEPAIRS(ch_assemble)
        ch_assembled = PRESTO_ASSEMBLEPAIRS.out.reads
    } else {
        PEAR(ch_assemble)
        ch_assembled = PEAR.out.reads
    }

    PRESTO_COLLAPSESEQ(ch_assembled)

    PRESTO_SPLITSEQ(PRESTO_COLLAPSESEQ.out.reads)

    emit:
    fasta      = PRESTO_SPLITSEQ.out.fasta
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
}
