include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    main:
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { row ->
            def meta = [
                id         : row[0].id,
                subject_id : row[0].subject_id,
                species    : row[0].species,
                locus      : row[0].locus,
                single_end : false
            ]
            [ meta, [ file(row[1]), file(row[2]) ] ]
        }
        .set { ch_reads }

    emit:
    reads = ch_reads
}
