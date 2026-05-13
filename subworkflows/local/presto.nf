include { FASTP                  } from '../../modules/nf-core/fastp/main'
include { GUNZIP                 } from '../../modules/local/presto/gunzip/main'
include { PRESTO_MASKPRIMERS_ALIGN } from '../../modules/local/presto/maskprimers_align/main'
include { PRESTO_PAIRSEQ         } from '../../modules/local/presto/pairseq/main'
include { ASSEMBLEPAIRS          } from '../../modules/local/assemblepairs/main'
include { PRESTO_COLLAPSESEQ     } from '../../modules/local/presto/collapseseq/main'
include { PRESTO_SPLITSEQ        } from '../../modules/local/presto/splitseq/main'

workflow PRESTO {
    take:
    ch_reads    // [meta, [R1, R2]]
    ch_cprimers // cprimers.fasta

    main:
    // 1. Adapter trimming + 3' quality trim
    // fastp module expects [meta, reads, adapter_fasta] — pass [] for adapter_fasta
    ch_reads_fastp = ch_reads.map { meta, reads -> [ meta, reads, [] ] }
    FASTP(ch_reads_fastp, false, false, false)

    // 2. Decompress fastp output — pRESTO tools need uncompressed FASTQ
    GUNZIP(FASTP.out.reads)

    // 3. Split R1 (V-read) and R2 (C-read) into separate channels
    ch_r1 = GUNZIP.out.reads.map { meta, reads -> [ meta, reads[0] ] }
    ch_r2 = GUNZIP.out.reads.map { meta, reads -> [ meta, reads[1] ] }

    // 4. MaskPrimers align on R2 only — finds CH1 primer by alignment,
    //    cuts it along with the variable 2/4/6bp stagger offset upstream
    //    R1 has no primer (V-gene primers enzymatically removed before sequencing)
    PRESTO_MASKPRIMERS_ALIGN(ch_r2, ch_cprimers.collect())

    // 5. PairSeq — synchronize R1 + R2, propagate C_CALL annotation
    ch_paired = ch_r1
        .join(PRESTO_MASKPRIMERS_ALIGN.out.reads, by: [0])

    PRESTO_PAIRSEQ(ch_paired)

    // 6. AssemblePairs sequential
    //    -1 = C-read (R2_pair-pass), -2 = V-read (R1_pair-pass)
    //    blastn reference: $IGDATA/fasta/imgt_{species}_ig_v.fasta inside container
    ASSEMBLEPAIRS(PRESTO_PAIRSEQ.out.reads)

    // 7. CollapseSeq — remove PCR duplicates
    PRESTO_COLLAPSESEQ(ASSEMBLEPAIRS.out.reads)

    // 8. SplitSeq — filter by DUPCOUNT >= splitseq_min_count; outputs FASTA
    PRESTO_SPLITSEQ(PRESTO_COLLAPSESEQ.out.reads)

    emit:
    fasta      = PRESTO_SPLITSEQ.out.fasta
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
}
