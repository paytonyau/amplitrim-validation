# Where do read pairs fail to merge in DADA2? Runs the standard steps on a folder of paired FASTQs at one truncation
# setting, then merges the same denoised pairs allowing 0, 1 and 2 mismatches and classifies the pairs that fail at 0.
# Usage: Rscript dada2_merge_check.R <input_dir> <truncLenF> <truncLenR> [max_samples=3] [maxEE=2 or F,R] [minOverlap=12]
suppressMessages(library(dada2))
a <- commandArgs(TRUE); tl <- c(as.integer(a[2]), as.integer(a[3]))
max_samples <- if (length(a) >= 4) as.integer(a[4]) else 3L
max_ee <- if (length(a) >= 5) as.numeric(strsplit(a[5], ",")[[1]]) else 2; if (length(max_ee) == 1) max_ee <- c(max_ee, max_ee); min_ov <- if (length(a) >= 6) as.integer(a[6]) else 12L
fn1 <- sort(list.files(a[1], pattern = "_R1.*\\.f(ast)?q(\\.gz)?$", full.names = TRUE)); fn2 <- sort(list.files(a[1], pattern = "_R2.*\\.f(ast)?q(\\.gz)?$", full.names = TRUE))
fn1 <- head(fn1, max_samples); fn2 <- head(fn2, max_samples)
d <- file.path(tempdir(), "f"); dir.create(d); f1 <- file.path(d, basename(fn1)); f2 <- file.path(d, basename(fn2))
out <- suppressWarnings(filterAndTrim(fn1, f1, fn2, f2, truncLen = tl, maxN = 0, maxEE = max_ee, truncQ = 2, rm.phix = FALSE, compress = TRUE))
keep <- file.exists(f1) & file.exists(f2); f1 <- f1[keep]; f2 <- f2[keep]
eF <- suppressMessages(learnErrors(f1, verbose = FALSE)); eR <- suppressMessages(learnErrors(f2, verbose = FALSE))
dF <- suppressMessages(dada(f1, err = eF, verbose = FALSE)); dR <- suppressMessages(dada(f2, err = eR, verbose = FALSE))
as_list <- function(x) if (is.data.frame(x) || inherits(x, "dada")) list(x) else x
reads_in <- function(x) sum(vapply(as_list(x), function(z) sum(getUniques(z)), 0))
n_in <- sum(out[, 1]); n_f <- sum(out[, 2])
cat(sprintf("%d samples, %d read pairs in, %.1f%% passed the filter; %.1f%% of them in denoised forward sequences, %.1f%% reverse\n",
            length(f1), n_in, 100 * n_f / n_in, 100 * reads_in(dF) / n_f, 100 * reads_in(dR) / n_f))
merge_at <- function(mm) as_list(suppressMessages(mergePairs(dF, f1, dR, f2, minOverlap = min_ov, maxMismatch = mm, returnRejects = TRUE, verbose = FALSE)))
for (mm in 0:2) {
  acc <- sum(vapply(merge_at(mm), function(m) sum(m$abundance[m$accept]), 0))
  cat(sprintf("merged with maxMismatch = %d: %.1f%% of the filtered read pairs\n", mm, 100 * acc / n_f))
}
# DADA2 does not report real overlap counts for pairs it rejects, so classify the failures by how many mismatches it takes
# to accept the same pair: merge again at 1, 2 ... 8 allowed mismatches and see where each pair first merges.
key <- function(m, id) paste(id, m$forward, m$reverse)
at0 <- merge_at(0)
pairs <- do.call(rbind, lapply(seq_along(at0), function(i) data.frame(key = key(at0[[i]], i), abundance = at0[[i]]$abundance, accept0 = at0[[i]]$accept)))
pairs$needs <- ifelse(pairs$accept0, 0, NA_real_)
for (mm in 1:8) {
  accepted <- do.call(c, lapply(seq_along(merge_at(mm)), function(i) { m <- merge_at(mm)[[i]]; key(m, i)[m$accept] }))
  pairs$needs[is.na(pairs$needs) & pairs$key %in% accepted] <- mm
}
fail <- pairs[!pairs$accept0, ]; tot <- sum(fail$abundance)
if (tot > 0) {
  share <- function(sel) 100 * sum(fail$abundance[sel]) / tot
  cat(sprintf("\nread pairs that fail at maxMismatch = 0: %.1f%% of the filtered read pairs. How many mismatches it takes to merge them:\n", 100 * tot / n_f))
  cat(sprintf("  1 mismatch:        %5.1f%%\n  2 mismatches:      %5.1f%%\n  3 or 4 mismatches: %5.1f%%\n  5 to 8 mismatches: %5.1f%%\n  more than 8 (or no usable overlap): %5.1f%%\n",
              share(fail$needs %in% 1), share(fail$needs %in% 2), share(fail$needs %in% 3:4), share(fail$needs %in% 5:8), share(is.na(fail$needs))))
}
