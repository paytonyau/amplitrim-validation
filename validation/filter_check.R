# The two DADA2 steps that longread_filter_check.py compares with AmpliTrim's estimate, on one file of single-end long reads.
#   Rscript filter_check.R filter  <reads.fastq.gz> <maxEE values, comma-separated> <table to write>
#       filterAndTrim with ONLY maxEE acting (truncQ 0, minQ 0, maxN 0, no length limits, no PhiX removal), once for each maxEE:
#       writes one tab-separated line per maxEE: maxEE, reads_in, reads_out.
#   Rscript filter_check.R primers <reads.fastq.gz> <forward primer> <reverse primer as ordered> <table to write>
#       removePrimers with its defaults and orient = TRUE (the way DADA2's own long-read workflow calls it):
#       writes one line: reads_in, reads_out.
suppressMessages(library(dada2))
a <- commandArgs(TRUE)
mode <- a[1]
tmp <- file.path(tempdir(), "filter_check")
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

if (mode == "filter") {
  fn <- a[2]
  grid <- as.numeric(strsplit(a[3], ",")[[1]])
  rows <- lapply(grid, function(g) {
    out <- file.path(tmp, paste0("maxEE_", g, ".fastq.gz"))
    res <- suppressWarnings(filterAndTrim(fn, out, maxEE = g, truncQ = 0, minQ = 0, maxN = 0, minLen = 0, maxLen = Inf,
                                          rm.phix = FALSE, compress = TRUE, multithread = FALSE))
    if (file.exists(out)) file.remove(out)
    data.frame(maxEE = g, reads_in = res[1, 1], reads_out = res[1, 2])
  })
  write.table(do.call(rbind, rows), a[4], sep = "\t", quote = FALSE, row.names = FALSE)
} else if (mode == "primers") {
  fn <- a[2]
  out <- file.path(tmp, "primers_removed.fastq.gz")
  res <- suppressWarnings(removePrimers(fn, out, primer.fwd = a[3], primer.rev = dada2:::rc(a[4]), orient = TRUE, verbose = FALSE))
  write.table(data.frame(reads_in = res[1, 1], reads_out = res[1, 2]), a[5], sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  stop("The first argument must be 'filter' or 'primers'.")
}
