#!/usr/bin/env Rscript
# Usage: run_scoper.R <outname> <method> <linkage> <threshold> <nproc> <tsv1> [tsv2 ...]
#   method:   'nt' (hierarchicalClones nucleotide hamming) | 'novj' (spectralClones)
#   linkage:  'single' | 'average' | 'complete'  (used only when method='nt')
#   threshold: numeric distance cutoff             (used only when method='nt')
#
# Outputs:
#   <outname>_clone-pass.tsv        — AIRR table with clone_id
#   <outname>_scoper_summary.txt    — text summary (seq/clone counts, thresholds)
#   <outname>_inter_intra.tsv       — pairwise inter/intra-clonal distances
#   <outname>_eff_threshold.tsv     — per-VJL-group effective thresholds
#   <outname>_vjl_groups.tsv        — VJL grouping table
#   <outname>_spectral_density.pdf  — density plot (inter vs intra + threshold)
#   <outname>_clone_summary.pdf     — clone size distribution
#   <outname>_clonality_measures.tsv — diversity/clonality/SHM indexes
#   versions.yml                    — tool versions

suppressPackageStartupMessages({
    library(alakazam)
    library(scoper)
    library(shazam)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    stop("Usage: run_scoper.R <outname> <method> <linkage> <threshold> <nproc> <tsv1> [tsv2 ...]")
}

outname   <- args[1]
method    <- args[2]
linkage   <- args[3]
threshold <- as.numeric(args[4])
nproc     <- as.integer(args[5])
tsv_files <- args[6:length(args)]

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
n_clones <- length(unique(na.omit(out_db$clone_id)))
message("Assigned ", n_clones, " clones")

write.table(out_db,
            file      = paste0(outname, "_clone-pass.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# Summary stats
tryCatch({
    sink(paste0(outname, "_scoper_summary.txt"))
    cat("Sample: ", outname, "\n")
    cat("Method: ", method, "\n")
    cat("Sequences: ", nrow(out_db), "\n")
    cat("Clones: ", n_clones, "\n\n")
    if (is_s4) {
        summary(result)
        et <- slot(result, "eff_threshold")
        cat("\n--- Effective thresholds per VJL group ---\n")
        if (length(et) > 0) {
            print(head(et, 50))
            cat("... (", length(et), " VJL groups total)\n")
        }
    }
    sink()
    message("Wrote summary: ", outname, "_scoper_summary.txt")
}, error = function(e) { try(sink(), silent=TRUE); message("Warning: summary failed: ", e$message) })

# Export raw data for downstream plotting
if (is_s4) {
    tryCatch({
        write.table(slot(result, "inter_intra"),
                    file = paste0(outname, "_inter_intra.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        message("Wrote inter/intra distances: ", outname, "_inter_intra.tsv")
    }, error = function(e) message("Warning: inter_intra export failed: ", e$message))

    tryCatch({
        et <- slot(result, "eff_threshold")
        write.table(data.frame(vjl_group = names(et), threshold = unname(et)),
                    file = paste0(outname, "_eff_threshold.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        message("Wrote effective thresholds: ", outname, "_eff_threshold.tsv")
    }, error = function(e) message("Warning: eff_threshold export failed: ", e$message))

    tryCatch({
        write.table(as.data.frame(slot(result, "vjl_groups")),
                    file = paste0(outname, "_vjl_groups.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        message("Wrote VJL groups: ", outname, "_vjl_groups.tsv")
    }, error = function(e) message("Warning: vjl_groups export failed: ", e$message))
}

# QC plots (PDF)
if (is_s4) {
    tryCatch({
        pdf(paste0(outname, "_spectral_density.pdf"), width = 10, height = 8)
        print(plot(result))
        dev.off()
        message("Wrote spectral density plot: ", outname, "_spectral_density.pdf")
    }, error = function(e) message("Warning: spectral density plot failed: ", e$message))
}

tryCatch({
    pdf(paste0(outname, "_clone_summary.pdf"), width = 10, height = 6)
    print(plotCloneSummary(out_db))
    dev.off()
    message("Wrote clone summary plot: ", outname, "_clone_summary.pdf")
}, error = function(e) message("Warning: clone summary plot failed: ", e$message))

# Clonality and diversity measures
tryCatch({
    prod_db <- out_db[out_db$productive == TRUE, ]
    n_total <- nrow(out_db)
    n_prod  <- nrow(prod_db)

    if (n_prod > 0 && "clone_id" %in% colnames(prod_db)) {
        cc <- countClones(prod_db, clone = "clone_id", copy = "duplicate_count")
        cc <- cc[order(-cc$seq_freq), ]
        n_clones_prod <- nrow(cc)

        p <- cc$seq_freq
        richness <- n_clones_prod

        shannon_H <- if (richness > 1) -sum(p[p > 0] * log(p[p > 0])) else 0
        simpson_lambda <- sum(p^2)
        gini_simpson <- 1 - simpson_lambda
        evenness_J <- if (richness > 1) shannon_H / log(richness) else 1
        clonality <- 1 - evenness_J

        s1 <- sum(cc$seq_count == 1)
        s2 <- sum(cc$seq_count == 2)
        chao1 <- if (s2 > 0) richness + (s1^2) / (2 * s2) else richness + s1 * (s1 - 1) / 2

        sorted <- sort(cc$seq_count)
        n_cc <- length(sorted)
        gini <- if (n_cc > 1) (2 * sum(seq_len(n_cc) * sorted) - (n_cc + 1) * sum(sorted)) / (n_cc * sum(sorted)) else 0

        top_id    <- as.character(cc$clone_id[1])
        top_freq  <- cc$seq_freq[1]
        top_count <- cc$seq_count[1]

        cumfreq <- cumsum(cc$seq_freq)
        d10     <- min(which(cumfreq >= 0.10))
        d50     <- min(which(cumfreq >= 0.50))
        d10_pct <- d10 / richness
        d50_pct <- d50 / richness

        mean_cs   <- mean(cc$seq_count)
        median_cs <- median(cc$seq_count)
        max_cs    <- max(cc$seq_count)

        shm_v_mean   <- NA_real_
        shm_cdr_r    <- NA_real_
        shm_cdr_s    <- NA_real_
        shm_fwr_r    <- NA_real_
        shm_fwr_s    <- NA_real_

        if ("sequence_alignment" %in% colnames(prod_db) && "germline_alignment" %in% colnames(prod_db)) {
            tryCatch({
                shm_db <- prod_db[nchar(prod_db$sequence_alignment) == nchar(prod_db$germline_alignment), ]
                message("SHM input: ", nrow(shm_db), "/", nrow(prod_db), " sequences with matched alignment lengths")
                mut_db <- observedMutations(shm_db,
                    sequenceColumn   = "sequence_alignment",
                    germlineColumn   = "germline_alignment",
                    regionDefinition = IMGT_V,
                    frequency = TRUE, combine = FALSE, nproc = nproc)
                shm_cdr_r <- mean(mut_db$mu_freq_cdr_r, na.rm = TRUE)
                shm_cdr_s <- mean(mut_db$mu_freq_cdr_s, na.rm = TRUE)
                shm_fwr_r <- mean(mut_db$mu_freq_fwr_r, na.rm = TRUE)
                shm_fwr_s <- mean(mut_db$mu_freq_fwr_s, na.rm = TRUE)
                shm_v_mean <- mean(shm_cdr_r + shm_cdr_s + shm_fwr_r + shm_fwr_s, na.rm = TRUE)
                message("Computed SHM frequencies (V-region, IMGT_V)")
            }, error = function(e) message("Warning: SHM computation failed: ", e$message))
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
}, error = function(e) message("Warning: clonality measures failed: ", e$message))

writeLines(paste0(
    '"SCOPER_HIERARCHICALCLONES":\n',
    '    scoper: ',   as.character(packageVersion("scoper")),   '\n',
    '    alakazam: ', as.character(packageVersion("alakazam")), '\n',
    '    shazam: ',   as.character(packageVersion("shazam")),   '\n',
    '    r-base: ',   paste(R.version$major, R.version$minor, sep = ".")
), con = "versions.yml")
