# biomarkers/coverage_shape/run.R — per-sample coverage-shape model of chrM depth (see run.sh)
#   Rscript run.R <corrected per-base bed.gz> <mt fasta> <out tsv> <bin> <exclude> <OL> <del START-END> <mt_len>
# Writes metric<TAB>value<TAB>unit<TAB>qc rows. Also usable as a library: source() then shape_fit().
read_mt <- function(fa, len) {
  s <- paste(readLines(fa)[-1], collapse = ""); toupper(substr(s, 1, len))
}
bin_depth <- function(pb, len, bin) {
  x <- read.table(gzfile(pb), col.names = c("chr", "s", "e", "d"), colClasses = c("character", "integer", "integer", "numeric"))
  dep <- numeric(len); for (i in seq_len(nrow(x))) { a <- x$s[i] + 1; b <- min(x$e[i], len); if (b >= a) dep[a:b] <- x$d[i] }
  dep
}
parse_ranges <- function(s) {  # "16024-576,3107-3107" -> list of c(start,end) (wrapping allowed)
  if (!nzchar(s)) return(list()); lapply(strsplit(strsplit(s, ",")[[1]], "-"), as.integer)
}
in_ranges <- function(pos, rr, len) {
  hit <- logical(length(pos))
  for (r in rr) hit <- hit | (if (r[1] <= r[2]) pos >= r[1] & pos <= r[2] else pos >= r[1] | pos <= r[2])
  hit
}
shape_fit <- function(dep, seq, bin = 100, exclude = "16024-576", OL = 5761, del = c(8470, 13447), len = 16569) {
  pos <- seq_len(len); b <- (pos - 1) %/% bin; nb <- max(b) + 1
  if (len %% bin < bin / 2) b[b == nb - 1] <- nb - 2      # merge a short last bin
  base <- strsplit(seq, "")[[1]]
  ok_pos <- !in_ranges(pos, parse_ranges(exclude), len) & base %in% c("A", "C", "G", "T")
  dd <- data.frame(b = b, d = dep, gc = base %in% c("G", "C"), ok = ok_pos, pos = pos)
  bins <- do.call(rbind, lapply(split(dd, dd$b), function(z) {
    k <- z[z$ok, ]; if (nrow(k) < 0.8 * nrow(z)) return(NULL)
    data.frame(mid = mean(z$pos), depth = mean(k$d), gc = mean(k$gc))
  }))
  bins <- bins[bins$depth > 0, ]
  med <- median(bins$depth)
  bins$y   <- log(bins$depth / med)
  bins$gcc <- (bins$gc - 0.44) * 10                       # per 10% GC, centred at mtDNA mean
  bins$e   <- ((bins$mid - OL) %% len) / len               # L-strand distance from O_L (0..1)
  bins$del <- as.numeric(bins$mid >= del[1] & bins$mid <= del[2])
  f <- MASS::rlm(y ~ gcc + I(gcc^2) + e + del, data = bins, maxit = 50)
  cf <- coef(f); r <- residuals(f); s <- mad(r)
  # GC-only and full fits for the R2 of each part
  r2 <- function(m) 1 - sum(residuals(m)^2) / sum((bins$y - mean(bins$y))^2)
  m_gc <- lm(y ~ gcc + I(gcc^2), bins); m_full <- lm(y ~ gcc + I(gcc^2) + e + del, bins)
  w <- which.min(r)
  list(bins = bins, fit = f, n_bins = nrow(bins), median_depth = med,
       gc_slope = unname(cf["gcc"]), gc_curv = unname(cf["I(gcc^2)"]),
       replication_gradient = unname(exp(-cf["e"])),     # depth just after O_L / just before (full circle)
       deletion_step = unname(exp(cf["del"])),           # depth inside common-deletion span vs outside
       r2_gc = r2(m_gc), r2_full = r2(m_full), residual_sd = s,
       local_dip = exp(r[w]), local_dip_pos = round(bins$mid[w]), local_dip_z = r[w] / s,
       gc_effect = function(g) { x <- (g - 0.44) * 10; cf["gcc"] * x + cf["I(gcc^2)"] * x^2 })
}
# GC-adjusted region ratio: raw depth ratio divided by the fitted GC effect difference
region_ratio <- function(dep, seq, fit, t, r) {
  gcof <- function(a, b) { z <- strsplit(substr(seq, a, b), "")[[1]]; mean(z %in% c("G", "C")) }
  raw <- mean(dep[t[1]:t[2]]) / mean(dep[r[1]:r[2]])
  c(raw = raw, gc_adj = unname(raw / exp(fit$gc_effect(gcof(t[1], t[2])) - fit$gc_effect(gcof(r[1], r[2])))))
}
if (sys.nframe() == 0) {
  a <- commandArgs(TRUE)
  len <- as.integer(a[8]); seq <- read_mt(a[2], len); dep <- bin_depth(a[1], len, as.integer(a[4]))
  del <- as.integer(strsplit(a[7], "-")[[1]])
  f <- shape_fit(dep, seq, as.integer(a[4]), a[5], as.integer(a[6]), del, len)
  nd <- region_ratio(dep, seq, f, c(12087, 12170), c(3485, 3553))
  mm <- region_ratio(dep, seq, f, del, c(3307, 5511))
  qc <- if (f$n_bins < 100 || f$median_depth < 50) "FAIL" else if (f$median_depth < 200) "LOW_DEPTH" else "PASS"
  out <- data.frame(metric = c("gc_slope", "gc_curvature", "replication_gradient", "deletion_step",
                               "shape_r2_gc", "shape_r2_full", "shape_residual_sd", "local_dip", "local_dip_pos", "local_dip_z",
                               "ND4_ND1_ratio_raw", "ND4_ND1_ratio_gc_adj", "MAJOR_MINOR_ratio_raw", "MAJOR_MINOR_ratio_gc_adj", "shape_bins"),
    value = c(f$gc_slope, f$gc_curv, f$replication_gradient, f$deletion_step, f$r2_gc, f$r2_full, f$residual_sd,
              f$local_dip, f$local_dip_pos, f$local_dip_z, nd, mm, f$n_bins),
    unit = c("log_per_10pctGC", "log_per_10pctGC2", "ratio", "ratio", "r2", "r2", "log", "ratio", "position", "z",
             "ratio", "ratio", "ratio", "ratio", "count"))
  out$value <- signif(out$value, 5); out$qc <- qc
  write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
