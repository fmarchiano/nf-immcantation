#!/usr/bin/env Rscript
# Usage: find_threshold.R <outname> <fixed_threshold> <drift_tol> <nproc> <tsv1> [tsv2 ...]
#
# Drift-check only. Computes the data-driven clonal threshold (distToNearest
# valley) for one patient (cloneby group) and compares it to the fixed
# threshold the pipeline actually clusters with. The detected value is NEVER
# applied to clustering; it only flags whether the fixed threshold should be
# re-evaluated for this patient.

suppressPackageStartupMessages({
    library(alakazam)
    library(shazam)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: find_threshold.R <outname> <fixed_threshold> <drift_tol> <nproc> <tsv1> [tsv2 ...]")
}

outname   <- args[1]
fixed_thr <- as.numeric(args[2])
drift_tol <- as.numeric(args[3])
nproc     <- as.integer(args[4])
tsv_files <- args[5:length(args)]

db <- do.call(rbind, lapply(tsv_files, function(f) {
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               check.names = FALSE, quote = "")
}))
message("Loaded ", nrow(db), " sequences from ", length(tsv_files), " file(s)")

# Length-normalised hamming, to match the scale of clonal_threshold
db <- distToNearest(db,
                    sequenceColumn = "junction",
                    vCallColumn    = "v_call",
                    jCallColumn    = "j_call",
                    model          = "ham",
                    normalize      = "len",
                    nproc          = nproc)

dist_vals <- db$dist_nearest[!is.na(db$dist_nearest)]

detected <- NA_real_
if (length(dist_vals) >= 2) {
    thr_obj <- tryCatch(
        findThreshold(dist_vals, method = "density"),
        error = function(e) { message("findThreshold(density) failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(thr_obj) && !is.na(thr_obj@threshold)) {
        detected <- thr_obj@threshold
    }
}

# Distance histogram (so the valley can be eyeballed against both thresholds)
pdf(paste0(outname, "_distToNearest.pdf"), width = 7, height = 5)
if (length(dist_vals) >= 1) {
    hist(dist_vals, breaks = 50, col = "grey80", border = "white",
         main = paste0("distToNearest — ", outname),
         xlab = "length-normalised hamming distance to nearest neighbour")
    abline(v = fixed_thr, col = "blue", lwd = 2, lty = 1)
    if (!is.na(detected)) abline(v = detected, col = "red", lwd = 2, lty = 2)
    legend("topright", bty = "n",
           legend = c(paste0("fixed (clustered): ", fixed_thr),
                      if (!is.na(detected)) paste0("detected valley: ", round(detected, 4)) else "detected: NA"),
           col    = c("blue", "red"), lwd = 2, lty = c(1, 2))
} else {
    plot.new(); title(paste0("distToNearest — ", outname, "\n(no distances computed)"))
}
invisible(dev.off())

# Drift verdict
drift  <- if (is.na(detected)) NA_real_ else abs(detected - fixed_thr)
status <- if (is.na(detected)) {
    "UNDETERMINED"
} else if (drift > drift_tol) {
    "REEVALUATE"
} else {
    "OK"
}

lines <- c(
    paste0("patient\t",          outname),
    paste0("n_sequences\t",      nrow(db)),
    paste0("n_distances\t",      length(dist_vals)),
    paste0("fixed_threshold\t",  fixed_thr),
    paste0("detected_threshold\t", ifelse(is.na(detected), "NA", round(detected, 6))),
    paste0("drift\t",            ifelse(is.na(drift), "NA", round(drift, 6))),
    paste0("drift_tol\t",        drift_tol),
    paste0("status\t",           status)
)
writeLines(lines, con = paste0(outname, "_threshold.txt"))

if (status == "REEVALUATE") {
    message("WARNING [", outname, "]: detected valley ", round(detected, 4),
            " drifts from fixed threshold ", fixed_thr,
            " by ", round(drift, 4), " (> tol ", drift_tol,
            ") — consider re-evaluating clonal_threshold for this patient.")
} else {
    message("Threshold drift-check [", outname, "]: status=", status)
}

writeLines(paste0(
    '"SHAZAM_DISTTONEAREST":\n',
    '    shazam: ',   as.character(packageVersion("shazam")),   '\n',
    '    alakazam: ', as.character(packageVersion("alakazam")), '\n',
    '    r-base: ',   paste(R.version$major, R.version$minor, sep = ".")
), con = "versions.yml")
