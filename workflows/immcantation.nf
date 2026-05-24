include { INPUT_CHECK      } from '../subworkflows/local/input_check'
include { PRESTO           } from '../subworkflows/local/presto'
include { VDJ_ASSIGNMENT   } from '../subworkflows/local/vdj_assignment'
include { CLONAL_ANALYSIS  } from '../subworkflows/local/clonal_analysis'

workflow IMMCANTATION {

    // Validate required params
    if (!params.input)    error "Parameter --input is required"
    if (!params.outdir)   error "Parameter --outdir is required"
    if (!params.cprimers) error "Parameter --cprimers is required (C-region/isotype primers FASTA)"
    if (!params.vprimers) error "Parameter --vprimers is required (V-region primers FASTA)"

    ch_cprimers = Channel.fromPath(params.cprimers, checkIfExists: true)
    ch_vprimers = Channel.fromPath(params.vprimers, checkIfExists: true)

    // Parse samplesheet
    INPUT_CHECK()
    ch_reads = INPUT_CHECK.out.reads

    // QC: adapter trim + pRESTO masking + assembly
    PRESTO(ch_reads, ch_cprimers, ch_vprimers)

    // V(D)J gene assignment — uses $IGDATA + germlines inside container
    VDJ_ASSIGNMENT(PRESTO.out.fasta)

    // Clonal analysis — exact-match (Briney parity) or hierarchical (SCOPer)
    CLONAL_ANALYSIS(VDJ_ASSIGNMENT.out.airr_tab)
}
