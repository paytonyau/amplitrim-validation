# ITS: DADA2 with no truncLen on a folder of paired FASTQs, then mergePairs on the SAME denoised reads twice (default, and trimOverhang = TRUE), with the rejected
# pairs counted by reason. Writes one row per sample (counts of read pairs):
#   sample, input, filtered, reach (pairs of denoised sequences that enter the merge step), merged (accepted by the default merge), rej_short (overlap shorter than
#   minOverlap), rej_mismatch (mismatches or indels in the overlap), rej_other (rejected with a clean overlap of at least minOverlap: an overhang), merged_trim_overhang
# Usage: Rscript its_merge_check.R <input_dir> <minLen> <maxEE> <minOverlap> <table to write>
suppressMessages(library(dada2))
a <- commandArgs(TRUE)
in_dir <- a[1]; min_len <- as.integer(a[2]); max_ee <- as.numeric(a[3]); min_ov <- as.integer(a[4]); out_path <- a[5]
fn1 <- sort(list.files(in_dir, pattern = "_R1.*\\.f(ast)?q(\\.gz)?$", full.names = TRUE))
fn2 <- sort(list.files(in_dir, pattern = "_R2.*\\.f(ast)?q(\\.gz)?$", full.names = TRUE))
d <- file.path(tempdir(), "filt"); dir.create(d)
f1 <- file.path(d, basename(fn1)); f2 <- file.path(d, basename(fn2))
out <- suppressWarnings(filterAndTrim(fn1, f1, fn2, f2, maxN = 0, maxEE = c(max_ee, max_ee), truncQ = 2, minLen = min_len, rm.phix = FALSE, compress = TRUE))
keep <- file.exists(f1) & file.exists(f2); idx <- which(keep)
eF <- suppressMessages(learnErrors(f1[keep], verbose = FALSE)); eR <- suppressMessages(learnErrors(f2[keep], verbose = FALSE))
dF <- suppressMessages(dada(f1[keep], err = eF, pool = FALSE, verbose = FALSE)); dR <- suppressMessages(dada(f2[keep], err = eR, pool = FALSE, verbose = FALSE))
as_frames <- function(m) if (is.data.frame(m)) list(m) else m
merge_with <- function(trim) as_frames(suppressMessages(mergePairs(dF, f1[keep], dR, f2[keep], minOverlap = min_ov, maxMismatch = 0, returnRejects = TRUE,
                                                                   trimOverhang = trim, verbose = FALSE)))
m0 <- merge_with(FALSE); m1 <- merge_with(TRUE)
count <- function(m) {
  if (is.null(m) || nrow(m) == 0) return(c(0, 0, 0, 0, 0))
  acc <- m$accept; short <- !acc & m$nmatch < min_ov; mism <- !acc & !short & (m$nmismatch > 0 | m$nindel > 0); other <- !acc & !short & !mism
  c(sum(m$abundance), sum(m$abundance[acc]), sum(m$abundance[short]), sum(m$abundance[mism]), sum(m$abundance[other]))
}
rows <- lapply(seq_along(fn1), function(i) {
  j <- match(i, idx); c0 <- if (is.na(j)) rep(0, 5) else count(m0[[j]]); c1 <- if (is.na(j)) rep(0, 5) else count(m1[[j]])
  data.frame(sample = sub("_R1.*$", "", basename(fn1[i])), input = out[i, 1], filtered = out[i, 2], reach = c0[1], merged = c0[2], rej_short = c0[3], rej_mismatch = c0[4],
             rej_other = c0[5], merged_trim_overhang = c1[2], row.names = NULL)
})
tab <- do.call(rbind, rows)
write.table(tab, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d samples, %d input pairs, %d merged, %d merged with trimOverhang\n", nrow(tab), sum(tab$input), sum(tab$merged), sum(tab$merged_trim_overhang)))
