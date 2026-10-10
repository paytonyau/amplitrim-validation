suppressMessages({library(dada2); library(ShortRead)})
a <- commandArgs(TRUE)
fwd <- a[1]; rev <- a[2]; truth <- read.delim(a[3], stringsAsFactors = FALSE); outcsv <- a[4]
pairs <- a[5:length(a)]
ids <- as.character(id(readFastq(fwd)))
kind_of <- function(h) { p <- strsplit(sub(" .*$", "", h), "|", fixed = TRUE)[[1]]; if (length(p) >= 3) p[3] else "normal" }
kinds <- vapply(ids, kind_of, ""); N <- length(kinds); nk <- table(kinds)
res <- list()
for (p in pairs) {
  tt <- as.integer(strsplit(p, ",")[[1]]); ff <- tempfile(fileext = ".fastq.gz"); fr <- tempfile(fileext = ".fastq.gz")
  fl <- suppressWarnings(filterAndTrim(fwd, ff, rev, fr, truncLen = tt, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                                      rm.phix = FALSE, compress = TRUE, multithread = FALSE))
  if (fl[1, 2] == 0) { cat("no reads left for", p, "\n"); next }
  eF <- suppressMessages(learnErrors(ff, multithread = FALSE, verbose = FALSE))
  eR <- suppressMessages(learnErrors(fr, multithread = FALSE, verbose = FALSE))
  dF <- dada(ff, err = eF, multithread = FALSE, verbose = FALSE)
  dR <- dada(fr, err = eR, multithread = FALSE, verbose = FALSE)
  mg <- mergePairs(dF, ff, dR, fr, minOverlap = 12, maxMismatch = 0, verbose = FALSE)
  tab <- function(df) {                       # reads in ASVs that exactly equal a true template, by kind
    k <- truth$kind[match(df$sequence, truth$seq)]; k[is.na(k)] <- "other"
    v <- tapply(df$abundance, factor(k, c("normal", "repeat", "chimera", "other")), sum); v[is.na(v)] <- 0; v
  }
  t1 <- tab(mg)
  st <- makeSequenceTable(list(s = mg)); nc <- removeBimeraDenovo(st, method = "consensus", multithread = FALSE, verbose = FALSE)
  t2 <- tab(data.frame(sequence = colnames(nc), abundance = as.numeric(nc[1, ]), stringsAsFactors = FALSE))
  res[[p]] <- data.frame(trim_f = tt[1], trim_r = tt[2], n_in = fl[1, 1], filtered = fl[1, 2] / N,
    dada2_merged = sum(mg$abundance) / N, asvs = nrow(mg),
    exact_normal = t1[["normal"]] / N, exact_repeat = t1[["repeat"]] / N, exact_chimera = t1[["chimera"]] / N, other = t1[["other"]] / N,
    nochim_total = sum(nc) / N, nochim_valid = (t2[["normal"]] + t2[["repeat"]]) / N, nochim_chim_left = t2[["chimera"]] / N,
    recall_normal = t1[["normal"]] / nk[["normal"]],
    recall_repeat = if ("repeat" %in% names(nk)) t1[["repeat"]] / nk[["repeat"]] else NA)
  cat(sprintf("%s done: merged %.1f%%\n", p, 100 * sum(mg$abundance) / N))
}
write.csv(do.call(rbind, res), outcsv, row.names = FALSE)
