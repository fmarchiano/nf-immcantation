#!/usr/bin/env Rscript
# Usage: run_scoper.R <outname> <method> <linkage> <threshold> <nproc> <tsv1> [tsv2 ...]
#   method:   'nt' (hierarchicalClones nucleotide hamming) | 'novj' (spectralClones)
#   linkage:  'single' | 'average' | 'complete'  (used only when method='nt')
#   threshold: numeric distance cutoff             (used only when method='nt')

suppressPackageStartupMessages({
    library(alakazam)
    library(scoper)
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
    spectralClones(db, method = "novj", nproc = nproc, summarize_clones = FALSE)
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
                       summarize_clones = FALSE)
}

out_db <- tryCatch(as.data.frame(result), error = function(e) result@db)
message("Assigned ", length(unique(na.omit(out_db$clone_id))), " clones")

write.table(out_db,
            file      = paste0(outname, "_clone-pass.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

writeLines(paste0(
    '"SCOPER_HIERARCHICALCLONES":\n',
    '    scoper: ',   as.character(packageVersion("scoper")),   '\n',
    '    alakazam: ', as.character(packageVersion("alakazam")), '\n',
    '    r-base: ',   paste(R.version$major, R.version$minor, sep = ".")
), con = "versions.yml")
