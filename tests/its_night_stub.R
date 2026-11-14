# A stand-in for the dada2 functions that its_night.R uses, for the tests only (ITS_NIGHT_STUB). It keeps the shapes of the real return values:
# dada() gives one object for one file and a list for several; a sample that keeps no read gets no filtered file; mergePairs() refuses a pair whose overlap is shorter
# than minOverlap; removeBimeraDenovo() removes the sequences whose name starts with CHIM. Every call writes its settings to ITS_NIGHT_LOG.
log_call <- function(...) { p <- Sys.getenv("ITS_NIGHT_LOG"); if (nzchar(p)) cat(paste0(paste(..., sep = "\t"), "\n"), file = p, append = TRUE) }
nreads <- function(f) length(readLines(f)) %/% 4
filterAndTrim <- function(fwd, filt, rev = NULL, filt.rev = NULL, truncLen = 0, maxN = 0, maxEE = Inf, truncQ = 2, minLen = 20, rm.phix = TRUE, compress = TRUE, ...) {
  log_call("filterAndTrim", if (is.null(rev)) "single" else "paired", paste(truncLen, collapse = ","), paste(maxEE, collapse = ","), truncQ, minLen)
  inn <- vapply(fwd, nreads, numeric(1))
  outn <- ifelse(grepl("empty", basename(fwd)), 0, pmax(0, inn - pmax(0, 6 - sum(maxEE))))       # a looser maxEE keeps more reads
  for (i in seq_along(fwd)) if (outn[i] > 0) { writeLines("x", filt[i]); if (!is.null(rev)) writeLines("x", filt.rev[i]) }
  m <- cbind(reads.in = inn, reads.out = outn); rownames(m) <- basename(fwd); m
}
learnErrors <- function(fls, multithread = FALSE, verbose = FALSE, ...) { log_call("learnErrors", length(fls)); list(n = length(fls)) }
dada <- function(derep, err, pool = FALSE, multithread = FALSE, verbose = FALSE, ...) {
  one <- function(f) structure(list(denoised = c(SEQA = 60L, SEQB = 30L, CHIM1 = 10L)), class = "dada")
  if (length(derep) == 1) one(derep[1]) else lapply(derep, one)
}
getUniques <- function(x, ...) x$denoised
mergePairs <- function(dadaF, derepF, dadaR, derepR, minOverlap = 12, maxMismatch = 0, returnRejects = FALSE, trimOverhang = FALSE, verbose = FALSE, ...) {
  log_call("mergePairs", minOverlap, maxMismatch)
  one <- function() {
    df <- data.frame(sequence = c("MERGEDA", "MERGEDB", "CHIM2", "MERGEDC"), abundance = c(50L, 25L, 10L, 15L), overlap = c(40, 25, 40, 15), stringsAsFactors = FALSE)
    df$accept <- df$overlap >= minOverlap
    if (returnRejects) df else df[df$accept, ]
  }
  if (inherits(dadaF, "dada")) one() else lapply(seq_along(dadaF), function(i) one())
}
makeSequenceTable <- function(samples, ...) {
  seqs <- sort(unique(unlist(lapply(samples, function(d) d$sequence))))
  st <- matrix(0L, nrow = length(samples), ncol = length(seqs), dimnames = list(names(samples), seqs))
  for (s in names(samples)) { d <- samples[[s]]; st[s, d$sequence] <- d$abundance }
  st
}
removeBimeraDenovo <- function(unqs, method = "consensus", multithread = FALSE, verbose = FALSE, ...) {
  log_call("removeBimeraDenovo", method)
  unqs[, !grepl("^CHIM", colnames(unqs)), drop = FALSE]
}
