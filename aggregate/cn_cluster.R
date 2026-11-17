#!/usr/bin/env Rscript
# aggregate/cn_cluster.R — mtDNA copy-number clusters across a cohort (mitrix cohort).
#
# Samples of different DNA source (e.g. whole blood vs cell line / PBMC) or batch often form
# separate groups of mtDNA copy number. Fits Gaussian mixtures (mclust, univariate E and V
# models) to log(mtDNA_CN) and adds to cohort.tsv:
#   cn_cluster       1 = lowest mtDNA-CN group, 2 = next, ...; all 1 when one group fits
#   cn_cluster_prob  posterior probability of the assigned cluster
# and the qc_flag token "cn_cluster:uncertain" when that probability is below --min-prob.
# Clusters are possible sample-type or batch groups, not labelled tissues.
#
# K: "auto" picks K in 1..--max-k by ICL (Biernacki et al. 2000), which penalises
# overlapping components, so one skewed sample type is not split; a K is accepted only if
# every cluster holds >= max(--min-n, --min-frac * n) samples. A number fixes K (as
# Mclust(x, K)). Writes <prefix>.tsv (one row per cluster) and <prefix>_fit.tsv (ICL / BIC
# of every K and model).
#
#   Rscript cn_cluster.R cohort.tsv out_prefix [auto|K] [max_k 4] [min_frac 0.05] [min_n 10] [min_prob 0.9]

args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: cn_cluster.R cohort.tsv out_prefix [auto|K] [max_k] [min_frac] [min_n] [min_prob]")
cohort   <- args[1]
prefix   <- args[2]
mode     <- if (length(args) >= 3) args[3] else "auto"
max_k    <- if (length(args) >= 4) as.integer(args[4]) else 4L
min_frac <- if (length(args) >= 5) as.numeric(args[5]) else 0.05
min_n    <- if (length(args) >= 6) as.integer(args[6]) else 10L
min_prob <- if (length(args) >= 7) as.numeric(args[7]) else 0.9
models   <- c("E", "V")

suppressPackageStartupMessages(library(mclust))
msg <- function(...) message("cn_cluster: ", ...)

d <- read.delim(cohort, colClasses = "character", check.names = FALSE, na.strings = character())
if (!"mtDNA_CN" %in% names(d)) stop("no mtDNA_CN column in ", cohort)
cn <- suppressWarnings(as.numeric(d$mtDNA_CN))
ok <- is.finite(cn) & cn > 0 & !grepl("(^|;)mtDNA_CN:", d$qc_flag)
x <- log(cn[ok]); n <- length(x)
need <- max(min_n, ceiling(min_frac * n))

sizes_ok <- function(f) f$G == 1 || all(tabulate(f$classification, f$G) >= need)

fit_tab <- data.frame()
if (n < 2 * need || length(unique(x)) < 3) {
  msg(n, " samples with mtDNA_CN: too few to cluster; all in cluster 1")
  best <- NULL; chosen <- "too_few"
} else if (mode == "auto") {
  ks <- seq_len(max_k)
  icl <- mclustICL(x, G = ks, modelNames = models, verbose = FALSE)
  bic <- mclustBIC(x, G = ks, modelNames = models, verbose = FALSE)
  cand <- expand.grid(K = ks, model = models, stringsAsFactors = FALSE)
  cand$ICL <- mapply(function(k, m) icl[as.character(k), m], cand$K, cand$model)
  cand$BIC <- mapply(function(k, m) bic[as.character(k), m], cand$K, cand$model)
  cand <- cand[is.finite(cand$ICL) & !(cand$K == 1 & cand$model == "V"), ]   # K=1: E and V are the same
  cand <- cand[order(-cand$ICL), ]
  cand$min_cluster_n <- NA_integer_; cand$selected <- FALSE
  best <- NULL
  for (i in seq_len(nrow(cand))) {
    f <- if (cand$K[i] == 1) Mclust(x, G = 1, verbose = FALSE) else
      Mclust(x, G = cand$K[i], modelNames = cand$model[i], verbose = FALSE)
    if (is.null(f)) next
    cand$min_cluster_n[i] <- min(tabulate(f$classification, f$G))
    if (is.null(best) && sizes_ok(f)) { best <- f; cand$selected[i] <- TRUE }
  }
  fit_tab <- cand; chosen <- "ICL"
} else {
  k <- as.integer(mode)
  if (is.na(k) || k < 1) stop("K must be 'auto' or a positive integer: ", mode)
  best <- Mclust(x, G = k, modelNames = if (k == 1) NULL else models, verbose = FALSE)
  if (is.null(best)) stop("mclust could not fit K = ", k)
  fit_tab <- data.frame(K = k, model = best$modelName, ICL = icl(best), BIC = best$bic,
                        min_cluster_n = min(tabulate(best$classification, best$G)), selected = TRUE)
  chosen <- "fixed"
}

cl <- rep(NA_integer_, nrow(d)); pr <- rep(NA_real_, nrow(d))
if (is.null(best) || best$G == 1) {
  cl[ok] <- 1L; pr[ok] <- 1; K <- 1L; model <- if (is.null(best)) "none" else best$modelName
} else {
  K <- best$G; model <- best$modelName
  ord <- order(best$parameters$mean)                       # 1 = lowest mtDNA-CN
  z <- best$z[, ord, drop = FALSE]
  cl[ok] <- max.col(z, ties.method = "first"); pr[ok] <- apply(z, 1, max)
}
d$cn_cluster <- ifelse(is.na(cl), "NA", as.character(cl))
d$cn_cluster_prob <- ifelse(is.na(pr), "NA", sprintf("%.3f", pr))
unc <- !is.na(pr) & pr < min_prob
d$qc_flag[unc] <- ifelse(d$qc_flag[unc] == "PASS", "cn_cluster:uncertain",
                         paste0(d$qc_flag[unc], ";cn_cluster:uncertain"))
# keep qc_flag last
d <- d[, c(setdiff(names(d), "qc_flag"), "qc_flag")]
write.table(d, cohort, sep = "\t", quote = FALSE, row.names = FALSE)

sm <- do.call(rbind, lapply(seq_len(K), function(k) {
  v <- cn[ok][cl[ok] == k]
  data.frame(cluster = k, n = length(v), fraction = round(length(v) / n, 3),
             cn_median = round(median(v), 1), cn_geomean = round(exp(mean(log(v))), 1),
             log_sd = round(sd(log(v)), 3), uncertain_n = sum(unc[ok][cl[ok] == k]),
             K = K, model = model, chosen_by = chosen, min_cluster_n = need)
}))
write.table(sm, paste0(prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(fit_tab)) write.table(fit_tab, paste0(prefix, "_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
msg(n, " samples, K = ", K, " (", chosen, ", model ", model, "): ",
    paste(sprintf("cluster %d n=%d median CN %.0f", sm$cluster, sm$n, sm$cn_median), collapse = "; "),
    "; ", sum(unc), " uncertain (prob < ", min_prob, ")")
if (K > 1) msg("mtDNA-CN forms ", K, " groups: possible mixed DNA sources or batches; check cn_cluster before pooling")
