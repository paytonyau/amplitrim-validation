# ITS overnight run (TEST_REGISTER.md, the section of the ITS overnight run): DADA2 arms on one dataset, with chimera removal in every arm.
# Usage: Rscript its_night.R <arms.tsv> <out_dir> [multithread: TRUE or FALSE, default TRUE]
# arms.tsv has one row per arm, in the order they are run: arm, mode, reads, truncF, truncR, maxEEF, maxEER, truncQ, minOverlap
#   mode paired : filterAndTrim, learnErrors and dada on both reads, then mergePairs (maxMismatch = 0); `reads` is a folder of <sample>_R1 and <sample>_R2 FASTQs
#   mode single : the same on the R1 files only (forward reads alone, or pairs merged before DADA2); truncR and maxEER are not used
#   mode fwd    : the denoised forward reads of the paired arm named in `reads`, without merging
#   mode overlap: mergePairs again on the denoised reads of the paired arm named in `reads`, with this row's minOverlap
# truncF and truncR of 0 mean no truncation (the ITS setting). minLen = 50 and maxN = 0 in every arm.
# Writes, for every arm, <out_dir>/<arm>/samples.tsv (sample, input, filtered, denoised, merged, nonchim, asvs), <out_dir>/<arm>/seqs.tsv (sample, sequence, abundance:
# the sequences after chimera removal) and <out_dir>/<arm>/DONE. An arm whose DONE is there is not run again, so a stopped night can be started again.
# denoised: for a paired arm the pairs of denoised reads that enter the merge step, for a single arm the denoised reads; merged: the pairs merged (single: = denoised);
# nonchim: the reads left after removeBimeraDenovo (method "consensus"); asvs: sequences holding at least 0.1% of a sample's reads after chimera removal.
# A paired arm keeps its filtered reads (<arm>/filtered, for the stricter estimates of the merge) and its denoised reads (<arm>/objects.rds, for fwd and overlap).
stub <- Sys.getenv("ITS_NIGHT_STUB")
if (nzchar(stub)) source(stub) else suppressMessages(library(dada2))      # the stand-in is for the tests only
a <- commandArgs(TRUE)
arms <- read.delim(a[1], stringsAsFactors = FALSE, colClasses = "character")
out_dir <- a[2]
mt <- if (length(a) >= 3) as.logical(a[3]) else TRUE
MIN_LEN <- 50
int <- function(x) as.integer(x)
num <- function(x) as.numeric(x)
sname <- function(f) sub("_R[12].*$", "", basename(f))
as_list <- function(x, cls) if (inherits(x, cls)) list(x) else x
reads_of <- function(dir, mate) sort(list.files(dir, pattern = paste0("_R", mate, ".*\\.f(ast)?q(\\.gz)?$"), full.names = TRUE))

finish <- function(dir, samples, input, filtered, denoised, merged, tables) {
  # tables: a list over `samples` of data.frames (sequence, abundance) or NULL; chimeras are removed over all samples together, as in the DADA2 workflow
  have <- which(vapply(tables, function(t) !is.null(t) && nrow(t) > 0, logical(1)))
  nonchim <- setNames(rep(0, length(samples)), samples)
  seqs <- data.frame(sample = character(0), sequence = character(0), abundance = numeric(0))
  if (length(have) > 0) {
    lst <- tables[have]; names(lst) <- samples[have]
    st <- makeSequenceTable(lst)
    if (ncol(st) > 0) {
      nc <- removeBimeraDenovo(st, method = "consensus", multithread = mt, verbose = FALSE)
      for (s in rownames(nc)) {
        v <- setNames(as.numeric(nc[s, , drop = FALSE]), colnames(nc)); v <- v[v > 0]      # names kept when one sequence is left (an isolate often has one)
        nonchim[s] <- sum(v)
        if (length(v) > 0) seqs <- rbind(seqs, data.frame(sample = s, sequence = names(v), abundance = as.numeric(v), row.names = NULL))
      }
    }
  }
  asvs <- vapply(samples, function(s) { v <- seqs$abundance[seqs$sample == s]; t <- sum(v); if (t > 0) sum(v >= 0.001 * t) else 0 }, numeric(1))
  tab <- data.frame(sample = samples, input = input, filtered = filtered, denoised = denoised, merged = merged, nonchim = as.numeric(nonchim), asvs = asvs, row.names = NULL)
  write.table(tab, file.path(dir, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(seqs, file.path(dir, "seqs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  file.create(file.path(dir, "DONE"))
  cat(sprintf("  %s: %d samples, %d input, %d merged, %d after chimera removal\n", basename(dir), length(samples), sum(input), sum(merged), sum(tab$nonchim)))
}

merge_and_finish <- function(o, min_ov, dir) {
  n <- length(o$fn1); samples <- sname(o$fn1); idx <- which(o$keep)
  denoised <- merged <- rep(0, n); tables <- vector("list", n)
  if (length(idx) > 0) {
    mm <- as_list(mergePairs(o$dF, o$f1[o$keep], o$dR, o$f2[o$keep], minOverlap = min_ov, maxMismatch = 0, returnRejects = TRUE, trimOverhang = FALSE, verbose = FALSE),
                  "data.frame")
    for (j in seq_along(idx)) {
      m <- mm[[j]]
      if (is.null(m) || nrow(m) == 0) next
      i <- idx[j]
      denoised[i] <- sum(m$abundance); acc <- m[m$accept, , drop = FALSE]; merged[i] <- sum(acc$abundance)
      if (nrow(acc) > 0) {
        s <- tapply(acc$abundance, as.character(acc$sequence), sum)
        tables[[i]] <- data.frame(sequence = names(s), abundance = as.integer(s), stringsAsFactors = FALSE)
      }
    }
  }
  finish(dir, samples, o$out[, 1], o$out[, 2], denoised, merged, tables)
}

single_finish <- function(fn1, out, keep, dF, dir) {
  n <- length(fn1); samples <- sname(fn1); idx <- which(keep)
  denoised <- rep(0, n); tables <- vector("list", n)
  for (j in seq_along(idx)) {
    u <- getUniques(dF[[j]])
    if (length(u) == 0) next
    i <- idx[j]; denoised[i] <- sum(u)
    tables[[i]] <- data.frame(sequence = names(u), abundance = as.integer(u), stringsAsFactors = FALSE)
  }
  finish(dir, samples, out[, 1], out[, 2], denoised, denoised, tables)
}

run_paired <- function(r, dir) {
  fn1 <- reads_of(r$reads, 1); fn2 <- reads_of(r$reads, 2)
  stopifnot(length(fn1) > 0, length(fn1) == length(fn2), identical(sname(fn1), sname(fn2)))
  fd <- file.path(dir, "filtered"); dir.create(fd, showWarnings = FALSE, recursive = TRUE)
  f1 <- file.path(fd, basename(fn1)); f2 <- file.path(fd, basename(fn2))
  out <- suppressWarnings(filterAndTrim(fn1, f1, fn2, f2, truncLen = c(int(r$truncF), int(r$truncR)), maxN = 0, maxEE = c(num(r$maxEEF), num(r$maxEER)),
                                        truncQ = int(r$truncQ), minLen = MIN_LEN, rm.phix = FALSE, compress = TRUE))
  keep <- file.exists(f1) & file.exists(f2)
  o <- list(fn1 = fn1, out = out, keep = keep, f1 = f1, f2 = f2, dF = list(), dR = list())
  if (any(keep)) {
    eF <- suppressMessages(learnErrors(f1[keep], multithread = mt, verbose = FALSE)); eR <- suppressMessages(learnErrors(f2[keep], multithread = mt, verbose = FALSE))
    o$dF <- as_list(suppressMessages(dada(f1[keep], err = eF, pool = FALSE, multithread = mt, verbose = FALSE)), "dada")
    o$dR <- as_list(suppressMessages(dada(f2[keep], err = eR, pool = FALSE, multithread = mt, verbose = FALSE)), "dada")
  }
  saveRDS(o, file.path(dir, "objects.rds"))
  merge_and_finish(o, int(r$minOverlap), dir)
}

run_single <- function(r, dir) {
  fn1 <- reads_of(r$reads, 1)
  stopifnot(length(fn1) > 0)
  fd <- file.path(dir, "filtered"); dir.create(fd, showWarnings = FALSE, recursive = TRUE)
  f1 <- file.path(fd, basename(fn1))
  out <- suppressWarnings(filterAndTrim(fn1, f1, truncLen = int(r$truncF), maxN = 0, maxEE = num(r$maxEEF), truncQ = int(r$truncQ), minLen = MIN_LEN, rm.phix = FALSE,
                                        compress = TRUE))
  keep <- file.exists(f1); dF <- list()
  if (any(keep)) {
    e <- suppressMessages(learnErrors(f1[keep], multithread = mt, verbose = FALSE))
    dF <- as_list(suppressMessages(dada(f1[keep], err = e, pool = FALSE, multithread = mt, verbose = FALSE)), "dada")
  }
  single_finish(fn1, out, keep, dF, dir)
}

objects_of <- function(arm) {
  p <- file.path(out_dir, arm, "objects.rds")
  if (!file.exists(p)) stop(sprintf("the arm '%s' must run before the arms that reuse it (no %s)", arm, p))
  readRDS(p)
}

for (k in seq_len(nrow(arms))) {
  r <- arms[k, ]
  dir <- file.path(out_dir, r$arm)
  if (file.exists(file.path(dir, "DONE"))) { cat(sprintf("  %s: done before\n", r$arm)); next }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  if (r$mode == "paired") run_paired(r, dir)
  else if (r$mode == "single") run_single(r, dir)
  else if (r$mode == "fwd") { o <- objects_of(r$reads); single_finish(o$fn1, o$out, o$keep, o$dF, dir) }
  else if (r$mode == "overlap") merge_and_finish(objects_of(r$reads), int(r$minOverlap), dir)
  else stop(sprintf("unknown mode '%s' for the arm '%s'", r$mode, r$arm))
}
