include { FASTP                                              } from '../../modules/nf-core/fastp/main'
include { GUNZIP                                             } from '../../modules/local/presto/gunzip/main'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_MASKPRIMERS_C   } from '../../modules/local/presto/maskprimers_align/main'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_MASKPRIMERS_V   } from '../../modules/local/presto/maskprimers_align/main'
include { PRESTO_PAIRSEQ                                     } from '../../modules/local/presto/pairseq/main'
include { PRESTO_PAIRSEQ as PAIRSEQ_BARCODE                  } from '../../modules/local/presto/pairseq/main'
include { PRESTO_PAIRSEQ as PAIRSEQ_CONSENSUS                } from '../../modules/local/presto/pairseq/main'
include { PRESTO_BUILDCONSENSUS as BUILDCONSENSUS_V          } from '../../modules/local/presto/buildconsensus/main'
include { PRESTO_BUILDCONSENSUS as BUILDCONSENSUS_C          } from '../../modules/local/presto/buildconsensus/main'
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

    // 4a. MaskPrimers on C-read with isotype-specific primers.
    //     Header annotated with C_CALL (isotype); in UMI mode (--umi) also
    //     extracts the UMI 5' of the primer into BARCODE (see modules.config).
    PRESTO_MASKPRIMERS_C(ch_cread, ch_cprimers.collect(), 'C')

    // 4b. MaskPrimers on V-read with V-region primers
    //     Briney removed these with cutadapt post-assembly; we remove pre-assembly
    PRESTO_MASKPRIMERS_V(ch_vread, ch_vprimers.collect(), 'V')

    // join: V-read first (PairSeq -1 / PEAR -f), C-read second (PairSeq -2 / PEAR -r)
    ch_paired = PRESTO_MASKPRIMERS_V.out.reads
        .join(PRESTO_MASKPRIMERS_C.out.reads, by: [0])

    if (params.umi) {
        // 5. Copy BARCODE (+ C_CALL isotype) from the C-read onto the V-read so
        //    both mates can be consensus-built per UMI.
        PAIRSEQ_BARCODE(ch_paired)

        // 6. Per-UMI consensus on each mate (grouped by the BARCODE field)
        ch_v_bc = PAIRSEQ_BARCODE.out.reads.map { meta, r1, r2 -> [ meta, r1 ] }
        ch_c_bc = PAIRSEQ_BARCODE.out.reads.map { meta, r1, r2 -> [ meta, r2 ] }
        BUILDCONSENSUS_V(ch_v_bc, 'V')
        BUILDCONSENSUS_C(ch_c_bc, 'C')

        // 7. Re-sync the two consensus files (BuildConsensus may drop a UMI from
        //    one mate but not the other) before assembly.
        ch_cons = BUILDCONSENSUS_V.out.reads
            .join(BUILDCONSENSUS_C.out.reads, by: [0])
        PAIRSEQ_CONSENSUS(ch_cons)

        ch_assemble = PAIRSEQ_CONSENSUS.out.reads
    } else {
        // 5. PairSeq — re-synchronize C-read + V-read after primer masking.
        //    --2f C_CALL copies isotype from -2 (C-read) onto -1 (V-read) so it
        //    survives PEAR assembly.
        PRESTO_PAIRSEQ(ch_paired)
        ch_assemble = PRESTO_PAIRSEQ.out.reads
    }

    // 8. PEAR — merge paired-end reads by overlap
    PEAR(ch_assemble)

    // 9. CollapseSeq — exact-sequence PCR-duplicate collapse
    PRESTO_COLLAPSESEQ(PEAR.out.reads)

    // 10. SplitSeq — filter by DUPCOUNT >= splitseq_min_count; outputs FASTA
    PRESTO_SPLITSEQ(PRESTO_COLLAPSESEQ.out.reads)

    emit:
    fasta      = PRESTO_SPLITSEQ.out.fasta
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
}
