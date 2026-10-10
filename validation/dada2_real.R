# Run DADA2 at one truncation setting on a folder of paired FASTQs and print one tab-separated line:
#   share passing filterAndTrim, share merged, share after chimera removal (all of the input reads), ASVs, then where
#   the reads went: share in denoised forward / reverse sequences, share of filtered forward reads that are singletons.
# Usage: Rscript dada2_real.R <input_dir> <truncLenF> <truncLenR> <maxEE, one number or F,R> <minOverlap> [minLen] [per-sample table to write] [ASV FASTA to write]
suppressMessages(library(dada2))
a <- commandArgs(TRUE); tl <- c(as.integer(a[2]), as.integer(a[3])); maxEE <- as.numeric(strsplit(a[4], ",")[[1]]); if (length(maxEE) == 1) maxEE <- c(maxEE, maxEE); minOv <- as.integer(a[5]); minLen <- if (length(a) >= 6) as.integer(a[6]) else 0L
track_path <- if (length(a) >= 7) a[7] else ""
fn1 <- sort(list.files(a[1], pattern = "_R1.*\\.f(ast)?q(\\.gz)?$", full.names = TRUE))
fn2 <- sort(list.files(a[1], pattern = "_R2.*\\.f(ast)?q(\\.gz)?$", full.names = TRUE))
d <- file.path(tempdir(), "filt"); dir.create(d)
f1 <- file.path(d, basename(fn1)); f2 <- file.path(d, basename(fn2))
out <- suppressWarnings(filterAndTrim(fn1, f1, fn2, f2, truncLen = tl, maxN = 0, maxEE = maxEE, truncQ = 2, minLen = minLen,
                                      rm.phix = FALSE, compress = TRUE))
keep <- file.exists(f1) & file.exists(f2)
eF <- suppressMessages(learnErrors(f1[keep], verbose = FALSE)); eR <- suppressMessages(learnErrors(f2[keep], verbose = FALSE))
# The defaults are DADA2's own (per-sample inference, consensus chimera removal). overnight.py sets these to test pooled inference and other chimera settings.
pool_text <- Sys.getenv("AMPLITRIM_DADA2_POOL", "FALSE"); pool_arg <- if (pool_text %in% c("TRUE", "FALSE")) as.logical(pool_text) else pool_text
dF <- suppressMessages(dada(f1[keep], err = eF, pool = pool_arg, verbose = FALSE)); dR <- suppressMessages(dada(f2[keep], err = eR, pool = pool_arg, verbose = FALSE))
reads_in <- function(x) if (inherits(x, "dada")) sum(getUniques(x)) else sum(vapply(x, function(d) sum(getUniques(d)), 0))
denoised_f <- reads_in(dF); denoised_r <- reads_in(dR)
drp <- suppressMessages(derepFastq(f1[keep])); drp <- if (inherits(drp, "derep")) list(drp) else drp
singleton_f <- sum(vapply(drp, function(d) sum(d$uniques == 1), 0)) / sum(vapply(drp, function(d) sum(d$uniques), 0))
mg <- mergePairs(dF, f1[keep], dR, f2[keep], minOverlap = minOv, maxMismatch = 0, verbose = FALSE)
if (!is.data.frame(mg)) mg <- if (length(mg) == 1) mg[[1]] else mg
n_in <- sum(out[, 1])
write_track <- function(merged_by_sample, nonchim_by_sample) {      # one row per sample, in counts of read pairs
  if (!nzchar(track_path)) return(invisible())
  names_all <- sub("_R1.*$", "", basename(fn1))
  idx <- which(keep)
  get <- function(v) { z <- rep(0, length(fn1)); z[idx] <- v; z }
  den <- function(x) vapply(if (inherits(x, "dada")) list(x) else x, function(d) sum(getUniques(d)), 0)
  tab_out <- data.frame(sample = names_all, input = out[, 1], filtered = out[, 2], denoised_f = get(den(dF)), denoised_r = get(den(dR)),
                        merged = get(merged_by_sample), nonchim = get(nonchim_by_sample))
  write.table(tab_out, track_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
merged <- if (is.data.frame(mg)) sum(mg$abundance) else sum(vapply(mg, function(m) sum(m$abundance), 0))
merged_by <- if (is.data.frame(mg)) sum(mg$abundance) else vapply(mg, function(m) sum(m$abundance), 0)
if (merged == 0) {
  write_track(rep(0, length(f1)), rep(0, length(f1)))     # nothing merged (for example the truncation leaves no overlap): report zeros instead of failing
  cat(sprintf("%.5f\t%.5f\t%.5f\t%d\t%.5f\t%.5f\t%.5f\n", sum(out[, 2]) / n_in, 0, 0, 0, denoised_f / n_in, denoised_r / n_in, singleton_f))
  quit(status = 0)
}
tab <- makeSequenceTable(mg); chim_args <- list(tab, method = Sys.getenv("AMPLITRIM_DADA2_CHIMERA", "consensus"), verbose = FALSE)
if (nzchar(Sys.getenv("AMPLITRIM_DADA2_MINFOLD"))) chim_args$minFoldParentOverAbundance <- as.numeric(Sys.getenv("AMPLITRIM_DADA2_MINFOLD"))
nc <- do.call(removeBimeraDenovo, chim_args)
write_track(merged_by, unname(rowSums(nc)))
asv_path <- if (length(a) >= 8) a[8] else ""
if (nzchar(asv_path)) {                                                    # the non-chimeric ASVs, most abundant first
  sz <- colSums(nc); o <- order(sz, decreasing = TRUE)
  writeLines(as.vector(rbind(sprintf(">ASV%d;size=%d", seq_along(o), as.integer(sz[o])), colnames(nc)[o])), asv_path)
}
cat(sprintf("%.5f\t%.5f\t%.5f\t%d\t%.5f\t%.5f\t%.5f\n", sum(out[, 2]) / n_in, merged / n_in, sum(nc) / n_in, ncol(nc),
            denoised_f / n_in, denoised_r / n_in, singleton_f))
