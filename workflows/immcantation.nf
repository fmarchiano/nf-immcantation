include { INPUT_CHECK      } from '../subworkflows/local/input_check'
include { PRESTO           } from '../subworkflows/local/presto'
include { PRESTO_FAST      } from '../subworkflows/local/presto_fast'
include { VDJ_ASSIGNMENT   } from '../subworkflows/local/vdj_assignment'
include { CLONAL_ANALYSIS  } from '../subworkflows/local/clonal_analysis'

workflow IMMCANTATION {

    // Validate required params
    if (!params.input)    error "Parameter --input is required"
    if (!params.outdir)   error "Parameter --outdir is required"
    if (!params.cprimers) error "Parameter --cprimers is required (C-region/isotype primers FASTA)"
    if (!params.vprimers) error "Parameter --vprimers is required (V-region primers FASTA)"

    ch_cprimers  = Channel.fromPath(params.cprimers, checkIfExists: true)
    ch_vprimers  = Channel.fromPath(params.vprimers, checkIfExists: true)
    ch_igblast   = Channel.fromPath("${params.ref_dir}/igblast",   type: 'dir', checkIfExists: true).first()
    ch_germlines = Channel.fromPath("${params.ref_dir}/germlines", type: 'dir', checkIfExists: true).first()

    // Parse samplesheet
    INPUT_CHECK()
    ch_reads = INPUT_CHECK.out.reads

    // QC: adapter trim + pRESTO masking + assembly
    if (params.fast) {
        PRESTO_FAST(ch_reads, ch_cprimers, ch_vprimers)
        ch_presto_fasta = PRESTO_FAST.out.fasta
    } else {
        PRESTO(ch_reads, ch_cprimers, ch_vprimers)
        ch_presto_fasta = PRESTO.out.fasta
    }

    // V(D)J gene assignment — uses $IGDATA + germlines inside container
    VDJ_ASSIGNMENT(ch_presto_fasta, ch_igblast, ch_germlines)

    // Clonal analysis — method set by params.cloning_method
    // ('hierarchical' = SCOPer on this branch; 'exact' = DefineClones)
    if (!params.skip_clonal) {
        CLONAL_ANALYSIS(VDJ_ASSIGNMENT.out.airr_tab)
    }
}
