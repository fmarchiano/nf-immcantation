process SHAZAM_CLONALITYMEASURES {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper'

    input:
    tuple val(meta), path(tab)   // germ-pass.tsv from CreateGermlines

    output:
    tuple val(meta), path("*_shm-pass.tsv"),          emit: tab
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
        read.table(f, header = TRUE, sep = "\\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
    }))
    message("Loaded ", nrow(out_db), " sequences from ", length(tsv_files), " file(s)")

    # ── Per-sequence SHM frequencies ─────────────────────────────────
    prod_db <- out_db[out_db\$productive == TRUE, ]

    if (nrow(prod_db) > 0 &&
        "sequence_alignment" %in% colnames(prod_db) &&
        "germline_alignment_d_mask" %in% colnames(prod_db)) {

        shm_db <- prod_db[nchar(prod_db\$sequence_alignment) == nchar(prod_db\$germline_alignment_d_mask), ]
        message("SHM input: ", nrow(shm_db), "/", nrow(prod_db),
                " sequences with matched alignment lengths")

        if (nrow(shm_db) > 0) {
            shm_db <- observedMutations(shm_db,
                sequenceColumn   = "sequence_alignment",
                germlineColumn   = "germline_alignment_d_mask",
                regionDefinition = IMGT_V,
                frequency = TRUE, combine = FALSE, nproc = nproc)
            message("Computed per-sequence SHM frequencies")

            # Merge mu_freq columns back into full dataframe
            mu_cols <- grep("^mu_freq_", colnames(shm_db), value = TRUE)
            out_db[rownames(shm_db), mu_cols] <- shm_db[, mu_cols]
        }
    } else if (nrow(prod_db) > 0 &&
               "sequence_alignment" %in% colnames(prod_db) &&
               "germline_alignment" %in% colnames(prod_db)) {
        # Fallback: use germline_alignment if d_mask not available
        shm_db <- prod_db[nchar(prod_db\$sequence_alignment) == nchar(prod_db\$germline_alignment), ]
        message("SHM input (fallback germline_alignment): ", nrow(shm_db), "/", nrow(prod_db),
                " sequences with matched alignment lengths")

        if (nrow(shm_db) > 0) {
            shm_db <- observedMutations(shm_db,
                sequenceColumn   = "sequence_alignment",
                germlineColumn   = "germline_alignment",
                regionDefinition = IMGT_V,
                frequency = TRUE, combine = FALSE, nproc = nproc)
            message("Computed per-sequence SHM frequencies (fallback)")

            mu_cols <- grep("^mu_freq_", colnames(shm_db), value = TRUE)
            out_db[rownames(shm_db), mu_cols] <- shm_db[, mu_cols]
        }
    }

    # Write per-sequence output with mu_freq columns
    write.table(out_db, file = paste0(outname, "_shm-pass.tsv"),
                sep = "\\t", quote = FALSE, row.names = FALSE)
    message("Wrote per-sequence SHM: ", outname, "_shm-pass.tsv")

    # ── Clone-level + sample-level summary ───────────────────────────
    tryCatch({
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

            # Clone-level SHM: mean mu_freq per clone, then patient-level summary
            # Three views: weighted mean (captures expansion), median (typical clone),
            # top clone (direct expansion flag)
            mu_cols <- grep("^mu_freq_", colnames(out_db), value = TRUE)
            if (length(mu_cols) > 0) {
                clone_shm <- aggregate(
                    out_db[out_db\$productive == TRUE, mu_cols],
                    by = list(clone_id = out_db[out_db\$productive == TRUE, "clone_id"]),
                    FUN = function(x) mean(x, na.rm = TRUE)
                )
                clone_shm <- merge(clone_shm, cc[, c("clone_id", "seq_freq")], by = "clone_id")

                # Weighted mean by clone frequency
                shm_cdr_r_w <- weighted.mean(clone_shm\$mu_freq_cdr_r, clone_shm\$seq_freq, na.rm = TRUE)
                shm_cdr_s_w <- weighted.mean(clone_shm\$mu_freq_cdr_s, clone_shm\$seq_freq, na.rm = TRUE)
                shm_fwr_r_w <- weighted.mean(clone_shm\$mu_freq_fwr_r, clone_shm\$seq_freq, na.rm = TRUE)
                shm_fwr_s_w <- weighted.mean(clone_shm\$mu_freq_fwr_s, clone_shm\$seq_freq, na.rm = TRUE)

                # Unweighted median across clones
                shm_cdr_r_med <- median(clone_shm\$mu_freq_cdr_r, na.rm = TRUE)
                shm_cdr_s_med <- median(clone_shm\$mu_freq_cdr_s, na.rm = TRUE)
                shm_fwr_r_med <- median(clone_shm\$mu_freq_fwr_r, na.rm = TRUE)
                shm_fwr_s_med <- median(clone_shm\$mu_freq_fwr_s, na.rm = TRUE)

                # Top clone SHM
                top_shm <- clone_shm[clone_shm\$clone_id == top_id, ]
                shm_cdr_r_top <- if (nrow(top_shm) > 0) top_shm\$mu_freq_cdr_r else NA_real_
                shm_cdr_s_top <- if (nrow(top_shm) > 0) top_shm\$mu_freq_cdr_s else NA_real_
                shm_fwr_r_top <- if (nrow(top_shm) > 0) top_shm\$mu_freq_fwr_r else NA_real_
                shm_fwr_s_top <- if (nrow(top_shm) > 0) top_shm\$mu_freq_fwr_s else NA_real_
            } else {
                shm_cdr_r_w <- shm_cdr_s_w <- shm_fwr_r_w <- shm_fwr_s_w <- NA_real_
                shm_cdr_r_med <- shm_cdr_s_med <- shm_fwr_r_med <- shm_fwr_s_med <- NA_real_
                shm_cdr_r_top <- shm_cdr_s_top <- shm_fwr_r_top <- shm_fwr_s_top <- NA_real_
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
                shm_freq_cdr_r_weighted  = round(shm_cdr_r_w, 6),
                shm_freq_cdr_r_median    = round(shm_cdr_r_med, 6),
                shm_freq_cdr_r_top_clone = round(shm_cdr_r_top, 6),
                shm_freq_cdr_s_weighted  = round(shm_cdr_s_w, 6),
                shm_freq_cdr_s_median    = round(shm_cdr_s_med, 6),
                shm_freq_cdr_s_top_clone = round(shm_cdr_s_top, 6),
                shm_freq_fwr_r_weighted  = round(shm_fwr_r_w, 6),
                shm_freq_fwr_r_median    = round(shm_fwr_r_med, 6),
                shm_freq_fwr_r_top_clone = round(shm_fwr_r_top, 6),
                shm_freq_fwr_s_weighted  = round(shm_fwr_s_w, 6),
                shm_freq_fwr_s_median    = round(shm_fwr_s_med, 6),
                shm_freq_fwr_s_top_clone = round(shm_fwr_s_top, 6),
                stringsAsFactors     = FALSE
            )

            write.table(measures, file = paste0(outname, "_clonality_measures.tsv"),
                        sep = "\\t", quote = FALSE, row.names = FALSE)
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
    touch ${meta.id}_shm-pass.tsv
    touch ${meta.id}_clonality_measures.tsv
    printf '"SHAZAM_CLONALITYMEASURES":\\n    shazam: 1.3.2\\n    alakazam: 1.4.3\\n    r-base: 4.4.2\\n' > versions.yml
    """
}
