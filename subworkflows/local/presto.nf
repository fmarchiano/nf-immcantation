include { FASTP                                              } from '../../modules/nf-core/fastp/main'
include { GUNZIP                                             } from '../../modules/local/presto/gunzip/main'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_MASKPRIMERS_C   } from '../../modules/local/presto/maskprimers_align/main'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_MASKPRIMERS_V   } from '../../modules/local/presto/maskprimers_align/main'
include { PRESTO_PAIRSEQ                                     } from '../../modules/local/presto/pairseq/main'
include { PEAR                                               } from '../../modules/local/pear/main'
include { PRESTO_COLLAPSESEQ                                 } from '../../modules/local/presto/collapseseq/main'
include { PRESTO_SPLITSEQ                                    } from '../../modules/local/presto/splitseq/main'

workflow PRESTO {
    take:
    ch_reads     // [meta, [R1, R2]]
    ch_cprimers  // C-region primers FASTA (isotype-specific: IgM, IgG, ...)
    ch_vprimers  // V-region primers FASTA (VH1-VH6)

    main:
    // 1. Adapter trimming + 3' quality trim
    ch_reads_fastp = ch_reads.map { meta, reads -> [ meta, reads, [] ] }
    FASTP(ch_reads_fastp, false, false, false)

    // 2. Decompress fastp output — pRESTO tools need uncompressed FASTQ
    GUNZIP(FASTP.out.reads)

    // 3. Split into C-read (reads[0]) and V-read (reads[1])
    //    Briney 2019 SRA dump has inverted convention: R1=C-side, R2=V-side
    ch_cread = GUNZIP.out.reads.map { meta, reads -> [ meta, reads[0] ] }
    ch_vread = GUNZIP.out.reads.map { meta, reads -> [ meta, reads[1] ] }

    // 4a. MaskPrimers on C-read with isotype-specific primers
    //     Header annotated with C_CALL (isotype) for downstream propagation
    PRESTO_MASKPRIMERS_C(ch_cread, ch_cprimers.collect(), 'C')

    // 4b. MaskPrimers on V-read with V-region primers
    //     Briney removed these with cutadapt post-assembly; we remove pre-assembly
    PRESTO_MASKPRIMERS_V(ch_vread, ch_vprimers.collect(), 'V')

    // 5. PairSeq — re-synchronize C-read + V-read after primer masking
    //    Order matters: V-read first (PairSeq -1, PEAR -f), C-read second (PairSeq -2, PEAR -r)
    //    --2f C_CALL copies isotype from -2 (C-read) onto -1 (V-read) so it survives PEAR assembly
    ch_paired = PRESTO_MASKPRIMERS_V.out.reads
        .join(PRESTO_MASKPRIMERS_C.out.reads, by: [0])

    PRESTO_PAIRSEQ(ch_paired)

    // 6. PEAR — merge paired-end reads by overlap
    PEAR(PRESTO_PAIRSEQ.out.reads)

    // 7. CollapseSeq — exact-sequence PCR-duplicate collapse
    PRESTO_COLLAPSESEQ(PEAR.out.reads)

    // 8. SplitSeq — filter by DUPCOUNT >= splitseq_min_count; outputs FASTA
    PRESTO_SPLITSEQ(PRESTO_COLLAPSESEQ.out.reads)

    emit:
    fasta      = PRESTO_SPLITSEQ.out.fasta
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
}
