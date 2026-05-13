#!/usr/bin/env nextflow

include { IMMCANTATION } from './workflows/immcantation'

workflow {
    IMMCANTATION()
}
