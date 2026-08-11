process SHAZAM_CLONALITYMEASURES {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper'

    input:
    tuple val(meta), path(tab)   // clone-pass.tsv from SCOPer

    output:
    tuple val(meta), path("*_clonality_measures.tsv"), emit: measures
    path "versions.yml",                               emit: versions

    script:
    def tab_list = tab instanceof List ? tab.join(' ') : tab
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(alakazam)
        library(shazam)
    })

    outname   <- "${meta.id}"
    nproc     <- as.integer("${task.cpus}")
    tsv_files <- strsplit("${tab_list}", " ")[[1]]

    out_db <- do.call(rbind, lapply(tsv_files, function(f) {
        read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
    }))
    message("Loaded ", nrow(out_db), " sequences from ", length(tsv_files), " file(s)")

    tryCatch({
        prod_db <- out_db[out_db\$productive == TRUE, ]
        n_total <- nrow(out_db)
        n_prod  <- nrow(prod_db)

        if (n_prod > 0 && "clone_id" %in% colnames(prod_db)) {
            cc <- countClones(prod_db, clone = "clone_id", copy = "duplicate_count")
            cc <- cc[order(-cc\$seq_freq), ]
            n_clones_prod <- nrow(cc)

            p <- cc\$seq_freq
            richness <- n_clones_prod

            shannon_H <- if (richness > 1) -sum(p[p > 0] * log(p[p > 0])) else 0
            simpson_lambda <- sum(p^2)
            gini_simpson <- 1 - simpson_lambda
            evenness_J <- if (richness > 1) shannon_H / log(richness) else 1
            clonality <- 1 - evenness_J

            s1 <- sum(cc\$seq_count == 1)
            s2 <- sum(cc\$seq_count == 2)
            chao1 <- if (s2 > 0) richness + (s1^2) / (2 * s2) else richness + s1 * (s1 - 1) / 2

            sorted <- sort(cc\$seq_count)
            n_cc <- length(sorted)
            gini <- if (n_cc > 1) (2 * sum(seq_len(n_cc) * sorted) - (n_cc + 1) * sum(sorted)) / (n_cc * sum(sorted)) else 0

            top_id    <- as.character(cc\$clone_id[1])
            top_freq  <- cc\$seq_freq[1]
            top_count <- cc\$seq_count[1]

            cumfreq <- cumsum(cc\$seq_freq)
            d10     <- min(which(cumfreq >= 0.10))
            d50     <- min(which(cumfreq >= 0.50))
            d10_pct <- d10 / richness
            d50_pct <- d50 / richness

            mean_cs   <- mean(cc\$seq_count)
            median_cs <- median(cc\$seq_count)
            max_cs    <- max(cc\$seq_count)

            shm_v_mean   <- NA_real_
            shm_cdr_r    <- NA_real_
            shm_cdr_s    <- NA_real_
            shm_fwr_r    <- NA_real_
            shm_fwr_s    <- NA_real_

            if ("sequence_alignment" %in% colnames(prod_db) && "germline_alignment" %in% colnames(prod_db)) {
                tryCatch({
                    shm_db <- prod_db[nchar(prod_db\$sequence_alignment) == nchar(prod_db\$germline_alignment), ]
                    message("SHM input: ", nrow(shm_db), "/", nrow(prod_db), " sequences with matched alignment lengths")
                    mut_db <- observedMutations(shm_db,
                        sequenceColumn   = "sequence_alignment",
                        germlineColumn   = "germline_alignment",
                        regionDefinition = IMGT_V,
                        frequency = TRUE, combine = FALSE, nproc = nproc)
                    shm_cdr_r <- mean(mut_db\$mu_freq_cdr_r, na.rm = TRUE)
                    shm_cdr_s <- mean(mut_db\$mu_freq_cdr_s, na.rm = TRUE)
                    shm_fwr_r <- mean(mut_db\$mu_freq_fwr_r, na.rm = TRUE)
                    shm_fwr_s <- mean(mut_db\$mu_freq_fwr_s, na.rm = TRUE)
                    shm_v_mean <- mean(shm_cdr_r + shm_cdr_s + shm_fwr_r + shm_fwr_s, na.rm = TRUE)
                    message("Computed SHM frequencies (V-region, IMGT_V)")
                }, error = function(e) message("Warning: SHM computation failed: ", e\$message))
            }

            measures <- data.frame(
                sample_id            = outname,
                total_sequences      = n_total,
                productive_sequences = n_prod,
                clone_count          = n_clones_prod,
                richness_chao1       = round(chao1, 2),
                shannon_entropy      = round(shannon_H, 6),
                simpson_index        = round(simpson_lambda, 6),
                gini_simpson         = round(gini_simpson, 6),
                evenness_pielou      = round(evenness_J, 6),
                clonality_index      = round(clonality, 6),
                gini_coefficient     = round(gini, 6),
                top_clone_id         = top_id,
                top_clone_freq       = round(top_freq, 6),
                top_clone_count      = top_count,
                d10                  = d10,
                d10_pct              = round(d10_pct, 6),
                d50                  = d50,
                d50_pct              = round(d50_pct, 6),
                mean_clone_size      = round(mean_cs, 2),
                median_clone_size    = median_cs,
                max_clone_size       = max_cs,
                shm_freq_v_mean      = round(shm_v_mean, 6),
                shm_freq_cdr_r_mean  = round(shm_cdr_r, 6),
                shm_freq_cdr_s_mean  = round(shm_cdr_s, 6),
                shm_freq_fwr_r_mean  = round(shm_fwr_r, 6),
                shm_freq_fwr_s_mean  = round(shm_fwr_s, 6),
                stringsAsFactors     = FALSE
            )

            write.table(measures, file = paste0(outname, "_clonality_measures.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            message("Wrote clonality measures: ", outname, "_clonality_measures.tsv")
        } else {
            message("Warning: no productive sequences with clone_id — skipping clonality measures")
        }
    }, error = function(e) message("Warning: clonality measures failed: ", e\$message))

    writeLines(paste0(
        '"SHAZAM_CLONALITYMEASURES":\\n',
        '    shazam: ',   as.character(packageVersion("shazam")),   '\\n',
        '    alakazam: ', as.character(packageVersion("alakazam")), '\\n',
        '    r-base: ',   paste(R.version\$major, R.version\$minor, sep = ".")
    ), con = "versions.yml")
    """

    stub:
    """
    touch ${meta.id}_clonality_measures.tsv
    printf '"SHAZAM_CLONALITYMEASURES":\\n    shazam: 1.3.2\\n    alakazam: 1.4.3\\n    r-base: 4.4.2\\n' > versions.yml
    """
}
