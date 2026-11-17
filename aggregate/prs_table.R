# aggregate/prs_table.R — long-format PRS table and optional validation (mitrix prs).
# Args: prs.sscore ancestry.tsv match.tsv min_overlap version out.tsv cohort.tsv|"" metadata|""
#       validation.tsv lead_spec ("score:file.sscore score:file.sscore ..." or "") covariates

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 11)
sscore <- args[1]; anc_file <- args[2]; match_file <- args[3]; min_ov <- as.numeric(args[4])
version <- args[5]; out_file <- args[6]; cohort_file <- args[7]; meta_file <- args[8]
val_file <- args[9]; lead_spec <- args[10]
meta_covars <- strsplit(args[11], ",", fixed = TRUE)[[1]]
fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.6g", x))

anc <- read.delim(anc_file, check.names = FALSE)
eur <- anc$sample[anc$ancestry == "EUR"]
read_sums <- function(f) {
  d <- read.delim(f, check.names = FALSE); names(d)[1] <- "sample"
  sums <- grep("_SUM$", names(d), value = TRUE)
  sums <- setdiff(sums, "NAMED_ALLELE_DOSAGE_SUM")
  out <- d[, c("sample", sums), drop = FALSE]
  names(out) <- c("sample", sub("_SUM$", "", sums)); out
}
sc <- read_sums(sscore)
if (nzchar(lead_spec)) for (x in strsplit(lead_spec, " ")[[1]]) {
  s <- sub(":.*", "", x); l <- read_sums(sub("^[^:]*:", "", x))
  names(l)[2] <- paste0("PRS_", s, "_lead")
  sc <- merge(sc, l, by = "sample", all.x = TRUE)
}
scores <- setdiff(names(sc), "sample")
mt <- read.delim(match_file)
low_ov <- mt$score[mt$fraction < min_ov]

# standardise within EUR samples (all samples if fewer than 2 EUR)
ref_set <- if (sum(sc$sample %in% eur) >= 2) sc$sample %in% eur else rep(TRUE, nrow(sc))
if (!all(ref_set) && sum(ref_set) < nrow(sc)) message("standardising within ", sum(ref_set), " EUR samples")
if (all(ref_set) && length(eur) < 2) message("warning: < 2 EUR samples; standardising within all samples")
z <- sc
for (s in scores) {
  m <- mean(sc[[s]][ref_set]); sdv <- sd(sc[[s]][ref_set])
  z[[s]] <- if (is.na(sdv) || sdv == 0) NA else (sc[[s]] - m) / sdv
}

ai <- match(sc$sample, anc$sample)
rows <- list()
for (i in seq_len(nrow(sc))) {
  s_id <- sc$sample[i]; a <- anc[ai[i], ]
  flag_s <- if (a$ancestry != "EUR") "non_EUR" else NULL
  for (s in scores) {
    base <- sub("_lead$", "", s)
    fl <- c(flag_s, if (base %in% low_ov) "low_overlap")
    q <- if (length(fl)) paste(fl, collapse = ";") else "PASS"
    rows[[length(rows) + 1]] <- c(s_id, s, fmt(z[[s]][i]), "z", q)
    rows[[length(rows) + 1]] <- c(s_id, paste0(s, "_sum"), fmt(sc[[s]][i]), "score", q)
  }
  rows[[length(rows) + 1]] <- c(s_id, "ancestry", a$ancestry, "label", "PASS")
  rows[[length(rows) + 1]] <- c(s_id, "ancestry_fraction", fmt(a$ancestry_fraction), "fraction", "PASS")
  for (k in 1:10) rows[[length(rows) + 1]] <- c(s_id, paste0("PC", k), fmt(a[[paste0("PC", k)]]), "PC", "PASS")
}
r <- do.call(rbind, rows)
out <- data.frame(sample = r[, 1], biomarker = "prs", metric = r[, 2], value = r[, 3],
                  unit = r[, 4], qc_flag = r[, 5], version = version)
write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("prs: ", nrow(sc), " samples, scores: ", paste(scores, collapse = ", "))

# ---- validation: each PRS vs the readout mitRix measures for the same trait
if (nzchar(cohort_file)) {
  co <- read.delim(cohort_file, check.names = FALSE, na.strings = "NA")
  readout <- c(PRS_CN = "mtDNA_CN", PRS_CN_raw = "mtDNA_CN", PRS_het = "age_accumulating", PRS_7S = "flap_7S", PRS_primer = "flap_primer")
  covs <- paste0("PC", 1:10)
  d <- merge(merge(z, anc[, c("sample", "ancestry", covs)], by = "sample"), co, by = "sample")
  if (nzchar(meta_file)) {
    me <- read.delim(meta_file, check.names = FALSE, na.strings = c("NA", ""))
    extra <- intersect(meta_covars, names(me))
    miss <- setdiff(meta_covars, names(me))
    if (length(miss)) message("validation: covariates not in metadata: ", paste(miss, collapse = ", "))
    if (length(extra)) { d <- merge(d, me[, c("sample", extra)], by = "sample"); covs <- c(covs, extra) }
    message("validation covariates: ", paste(covs, collapse = ", "))
  }
  d <- d[d$ancestry == "EUR", ]
  bt <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))
  val <- do.call(rbind, lapply(scores, function(s) {
    ro <- readout[[sub("_lead$", "", s)]]
    if (is.null(ro) || !(ro %in% names(d))) return(NULL)
    dd <- d[, c(ro, s, covs)]; dd[[ro]] <- suppressWarnings(as.numeric(dd[[ro]]))
    dd <- dd[complete.cases(dd), , drop = FALSE]
    if (nrow(dd) < 20) return(data.frame(score = s, readout = ro, n = nrow(dd), spearman_r = NA, p = NA,
                                         incremental_r2 = NA, note = "fewer than 20 EUR samples"))
    ct <- suppressWarnings(cor.test(dd[[s]], dd[[ro]], method = "spearman"))
    f0 <- reformulate(bt(covs), response = bt(ro)); f1 <- reformulate(bt(c(s, covs)), response = bt(ro))
    r2 <- summary(lm(f1, dd))$r.squared - summary(lm(f0, dd))$r.squared
    data.frame(score = s, readout = ro, n = nrow(dd), spearman_r = round(unname(ct$estimate), 4),
               p = signif(ct$p.value, 3), incremental_r2 = round(r2, 5), note = "")
  }))
  if (!is.null(val)) write.table(val, val_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
