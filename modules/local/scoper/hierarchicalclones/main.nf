process SCOPER_HIERARCHICALCLONES {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper'

    input:
    tuple val(meta), path(tabs)   // one or more AIRR TSV files grouped by subject

    output:
    tuple val(meta), path("*_clone-pass.tsv"), emit: tab
    path "*_scoper_summary.txt", optional: true, emit: summary
    path "*_inter_intra.tsv", optional: true, emit: inter_intra
    path "*_eff_threshold.tsv", optional: true, emit: eff_threshold
    path "*_vjl_groups.tsv", optional: true, emit: vjl_groups
    path "*_spectral_density.pdf", optional: true, emit: spectral_plot
    path "*_clone_summary.pdf", optional: true, emit: clone_plot
    path "versions.yml",                       emit: versions

    script:
    def method    = task.ext.method    ?: 'nt'
    def linkage   = task.ext.linkage   ?: 'complete'
    def threshold = task.ext.threshold ?: params.clonal_threshold
    def tab_list  = tabs instanceof List ? tabs.join(' ') : tabs
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(alakazam)
        library(scoper)
    })

    outname   <- "${meta.id}"
    method    <- "${method}"
    linkage   <- "${linkage}"
    threshold <- as.numeric("${threshold}")
    nproc     <- as.integer("${task.cpus}")
    tsv_files <- strsplit("${tab_list}", " ")[[1]]

    db <- do.call(rbind, lapply(tsv_files, function(f) {
        read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
    }))

    message("Loaded ", nrow(db), " sequences from ", length(tsv_files), " file(s)")

    result <- if (method == "novj") {
        spectralClones(db, method = "novj", nproc = nproc, summarize_clones = TRUE)
    } else {
        hierarchicalClones(db,
                           threshold = threshold,
                           method    = method,
                           linkage   = linkage,
                           v_call    = "v_call",
                           j_call    = "j_call",
                           junction  = "junction",
                           clone     = "clone_id",
                           nproc     = nproc,
                           summarize_clones = TRUE)
    }

    is_s4  <- isS4(result)
    out_db <- if (is_s4) as.data.frame(result) else result
    n_clones <- length(unique(na.omit(out_db\$clone_id)))
    message("Assigned ", n_clones, " clones")

    write.table(out_db,
                file      = paste0(outname, "_clone-pass.tsv"),
                sep       = "\t",
                quote     = FALSE,
                row.names = FALSE)

    # Summary stats
    tryCatch({
        sink(paste0(outname, "_scoper_summary.txt"))
        cat("Sample: ", outname, "\\n")
        cat("Method: ", method, "\\n")
        cat("Sequences: ", nrow(out_db), "\\n")
        cat("Clones: ", n_clones, "\\n\\n")
        if (is_s4) {
            summary(result)
            et <- slot(result, "eff_threshold")
            cat("\\n--- Effective thresholds per VJL group ---\\n")
            if (length(et) > 0) {
                print(head(et, 50))
                cat("... (", length(et), " VJL groups total)\\n")
            }
        }
        sink()
        message("Wrote summary: ", outname, "_scoper_summary.txt")
    }, error = function(e) { try(sink(), silent=TRUE); message("Warning: summary failed: ", e\$message) })

    # Export raw data for downstream plotting
    if (is_s4) {
        tryCatch({
            write.table(slot(result, "inter_intra"),
                        file = paste0(outname, "_inter_intra.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            message("Wrote inter/intra distances: ", outname, "_inter_intra.tsv")
        }, error = function(e) message("Warning: inter_intra export failed: ", e\$message))

        tryCatch({
            et <- slot(result, "eff_threshold")
            write.table(data.frame(vjl_group = names(et), threshold = unname(et)),
                        file = paste0(outname, "_eff_threshold.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            message("Wrote effective thresholds: ", outname, "_eff_threshold.tsv")
        }, error = function(e) message("Warning: eff_threshold export failed: ", e\$message))

        tryCatch({
            write.table(as.data.frame(slot(result, "vjl_groups")),
                        file = paste0(outname, "_vjl_groups.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            message("Wrote VJL groups: ", outname, "_vjl_groups.tsv")
        }, error = function(e) message("Warning: vjl_groups export failed: ", e\$message))
    }

    # QC plots (PDF)
    if (is_s4) {
        tryCatch({
            pdf(paste0(outname, "_spectral_density.pdf"), width = 10, height = 8)
            print(plot(result))
            dev.off()
            message("Wrote spectral density plot: ", outname, "_spectral_density.pdf")
        }, error = function(e) message("Warning: spectral density plot failed: ", e\$message))
    }

    tryCatch({
        pdf(paste0(outname, "_clone_summary.pdf"), width = 10, height = 6)
        print(plotCloneSummary(out_db))
        dev.off()
        message("Wrote clone summary plot: ", outname, "_clone_summary.pdf")
    }, error = function(e) message("Warning: clone summary plot failed: ", e\$message))

    writeLines(paste0(
        '"SCOPER_HIERARCHICALCLONES":\\n',
        '    scoper: ',   as.character(packageVersion("scoper")),   '\\n',
        '    alakazam: ', as.character(packageVersion("alakazam")), '\\n',
        '    r-base: ',   paste(R.version\$major, R.version\$minor, sep = ".")
    ), con = "versions.yml")
    """

    stub:
    """
    touch ${meta.id}_clone-pass.tsv
    printf '"SCOPER_HIERARCHICALCLONES":\\n    scoper: 1.5.0\\n    alakazam: 1.4.3\\n    r-base: 4.4.2\\n' > versions.yml
    """
}
