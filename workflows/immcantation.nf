include { INPUT_CHECK      } from '../subworkflows/local/input_check'
include { PRESTO           } from '../subworkflows/local/presto'
include { VDJ_ASSIGNMENT   } from '../subworkflows/local/vdj_assignment'
include { CLONAL_ANALYSIS  } from '../subworkflows/local/clonal_analysis'

workflow IMMCANTATION {

    // Validate required params
    if (!params.input)    error "Parameter --input is required"
    if (!params.outdir)   error "Parameter --outdir is required"
    if (!params.cprimers) error "Parameter --cprimers is required"

    ch_cprimers = Channel.fromPath(params.cprimers, checkIfExists: true)

    // Parse samplesheet
    INPUT_CHECK()
    ch_reads = INPUT_CHECK.out.reads

    // QC: adapter trim + pRESTO assembly
    PRESTO(ch_reads, ch_cprimers)

    // V(D)J gene assignment — uses $IGDATA + germlines inside container
    VDJ_ASSIGNMENT(PRESTO.out.fasta)

    // Clonal analysis with fixed threshold
    CLONAL_ANALYSIS(VDJ_ASSIGNMENT.out.airr_tab)
}
