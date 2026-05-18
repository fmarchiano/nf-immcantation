include { FASTP                  } from '../../modules/nf-core/fastp/main'
include { GUNZIP                 } from '../../modules/local/presto/gunzip/main'
include { PRESTO_MASKPRIMERS_ALIGN } from '../../modules/local/presto/maskprimers_align/main'
include { PRESTO_PAIRSEQ         } from '../../modules/local/presto/pairseq/main'
include { PEAR                   } from '../../modules/local/pear/main'
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

    // 4. MaskPrimers align on R1 only — R1 is the C-read in this dataset
    //    (reads from CH1 end toward VH); finds CH1 primer by alignment and cuts
    //    the variable upstream region (CDR3/JH) along with the primer itself
    PRESTO_MASKPRIMERS_ALIGN(ch_r1, ch_cprimers.collect())

    // 5. PairSeq — synchronize R2 (V-read) + masked R1 (C-read), propagate C_CALL annotation
    ch_paired = ch_r2
        .join(PRESTO_MASKPRIMERS_ALIGN.out.reads, by: [0])

    PRESTO_PAIRSEQ(ch_paired)

    // 6. PEAR — merge paired-end reads by overlap
    //    -f = C-read (R2_pair-pass, original Illumina R1), -r = V-read (R1_pair-pass, original Illumina R2)
    PEAR(PRESTO_PAIRSEQ.out.reads)

    // 7. CollapseSeq — remove PCR duplicates
    PRESTO_COLLAPSESEQ(PEAR.out.reads)

    // 8. SplitSeq — filter by DUPCOUNT >= splitseq_min_count; outputs FASTA
    PRESTO_SPLITSEQ(PRESTO_COLLAPSESEQ.out.reads)

    emit:
    fasta      = PRESTO_SPLITSEQ.out.fasta
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
}
