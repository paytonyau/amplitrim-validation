# DADA2's filterAndTrim with no truncLen on one pooled pair of FASTQs, run twice: as a pair, and with the forward reads alone. Prints one tab-separated line:
#   pairs in, pairs out, forward reads in, forward reads out
# Same settings as DADA2's ITS workflow and as AmpliTrim's ITS mode (maxN = 0, truncQ = 2, minLen); phiX is not removed, as in dada2_real.R.
# Usage: Rscript its_filter_check.R <R1> <R2> <minLen> <maxEE>
suppressMessages(library(dada2))
a <- commandArgs(TRUE); min_len <- as.integer(a[3]); max_ee <- as.numeric(a[4])
d <- file.path(tempdir(), "its_filt"); dir.create(d)
paired <- suppressWarnings(filterAndTrim(a[1], file.path(d, "p1.fastq.gz"), a[2], file.path(d, "p2.fastq.gz"), maxN = 0, maxEE = c(max_ee, max_ee), truncQ = 2,
                                         minLen = min_len, rm.phix = FALSE, compress = TRUE))
single <- suppressWarnings(filterAndTrim(a[1], file.path(d, "s1.fastq.gz"), maxN = 0, maxEE = max_ee, truncQ = 2, minLen = min_len, rm.phix = FALSE, compress = TRUE))
cat(sprintf("%d\t%d\t%d\t%d\n", paired[1, 1], paired[1, 2], single[1, 1], single[1, 2]))
