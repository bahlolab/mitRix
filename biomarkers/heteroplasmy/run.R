# biomarkers/heteroplasmy/run.R — called by run.sh; base R only.
# Args: variants.tsv chrM.per-base.bed.gz blocklist.tsv annotation.tsv.gz hf_min hf_max
#       hf_split min_alt_reads exclude_filters sample version out_results.tsv out_variants.tsv
#       numt_alleles.tsv|none numt_k numt_max_hf nuclear_depth|NA

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 17)
numt_file <- args[14]; numt_k <- as.numeric(args[15]); numt_max_hf <- as.numeric(args[16])
nuclear_depth <- suppressWarnings(as.numeric(args[17]))
variants_file <- args[1]; coverage_file <- args[2]
blocklist_file <- args[3]; annotation_file <- args[4]
hf_min <- as.numeric(args[5]); hf_max <- as.numeric(args[6]); hf_split <- as.numeric(args[7])
min_alt <- as.integer(args[8])
exclude <- strsplit(args[9], ",", fixed = TRUE)[[1]]
sample <- args[10]; version <- args[11]
out_results <- args[12]; out_variants <- args[13]

log <- function(...) message("heteroplasmy: ", ...)
fmt <- function(x) sprintf("%.10g", x)

# ---- inputs
v  <- read.delim(variants_file, colClasses = "character", na.strings = character())
bl <- read.delim(blocklist_file, colClasses = "character", comment.char = "#")
bl_start <- as.integer(bl$start); bl_end <- as.integer(bl$end)
ann <- read.delim(gzfile(annotation_file), colClasses = "character", na.strings = character())
cov <- read.delim(gzfile(coverage_file), header = FALSE, colClasses = c("character", "integer", "integer", "numeric"))
depth_per_base <- rep(cov[[4]], cov[[3]] - cov[[2]])
mt_median_depth <- median(depth_per_base)

# ---- calling
pos <- as.integer(v$pos)
hf  <- as.numeric(v$hf)
alt_reads <- suppressWarnings(as.integer(v$alt_reads))
blocked  <- vapply(pos, function(p) any(p >= bl_start & p <= bl_end), logical(1))
filtered <- vapply(strsplit(v$filter, ";", fixed = TRUE), function(f) any(f %in% exclude), logical(1))
in_range <- !is.na(hf) & hf >= hf_min & hf <= hf_max
low_reads <- !is.na(alt_reads) & alt_reads < min_alt
# NUMT-like: an allele found in NUMT sequences at no more than the HF NUMT reads can give,
# ~ copies x nuclear depth / chrM depth; k NUMT-copy equivalents (default 1.5). Without a
# nuclear depth, a fixed ceiling (numt_max_hf) is used.
numt_hf <- if (!is.na(nuclear_depth) && nuclear_depth > 0 && mt_median_depth > 0)
  min(1, numt_k * nuclear_depth / mt_median_depth) else numt_max_hf
numt_src <- if (!is.na(nuclear_depth) && nuclear_depth > 0) "nuclear depth" else "fixed"
numt_like <- rep(FALSE, nrow(v))
if (numt_file != "none") {
  na_ <- read.delim(numt_file, colClasses = "character", comment.char = "#")
  numt_like <- paste(v$pos, v$ref, v$alt) %in% paste(na_$pos, na_$ref, na_$alt) & !is.na(hf) & hf < numt_hf
}
keep <- in_range & !blocked & !filtered & !low_reads & !numt_like
if (any(in_range & is.na(alt_reads)))
  log("warning: ", sum(in_range & is.na(alt_reads)), " variant(s) without depth; --min-alt-reads not applied to them")
log(sum(in_range), " variants with ", hf_min, " <= HF <= ", hf_max, "; removed ",
    sum(in_range & blocked), " blocklisted, ", sum(in_range & !blocked & filtered),
    " with ", paste(exclude, collapse = "/"), " and ", sum(in_range & !blocked & !filtered & low_reads),
    " with < ", min_alt, " ALT reads, ", sum(in_range & !blocked & !filtered & !low_reads & numt_like),
    " NUMT-like (HF < ", fmt(numt_hf), ", ", numt_src, "); ", sum(keep), " heteroplasmies")

het <- v[keep, c("sample", "pos", "ref", "alt", "hf", "depth", "alt_reads", "filter", "type")]
p <- as.integer(het$pos); h <- as.numeric(het$hf)
snv <- het$type == "SNV"

# ---- strand and substitution: rCRS is the light strand; REF C/A = light-strand change,
# REF G/T = heavy-strand change, written as its complement (Gupta et al., Nature 2026)
comp <- c(A = "T", C = "G", G = "C", T = "A")
strand <- ifelse(!snv, NA, ifelse(het$ref %in% c("C", "A"), "L", ifelse(het$ref %in% c("G", "T"), "H", NA)))
sub_ref <- ifelse(strand == "H", comp[het$ref], het$ref)
sub_alt <- ifelse(strand == "H", comp[het$alt], het$alt)
substitution <- ifelse(is.na(strand), NA, paste0(sub_ref, ">", sub_alt))
transition <- substitution %in% c("C>T", "A>G")
in_ori <- p >= 16172 | p <= 210
age_acc <- !is.na(strand) & !in_ori &
  ((strand == "H" & substitution %in% c("C>T", "A>G")) | (strand == "L" & substitution == "A>G"))
hf_bin <- ifelse(h <= hf_split, "low", "high")
low <- hf_bin == "low"

# ---- annotation: SNVs by pos/ref/alt; other variants get gene and region by position
key <- paste(ann$pos, ann$ref, ann$alt)
m <- match(paste(het$pos, het$ref, het$alt), key)
mp <- match(het$pos, ann$pos)
gene    <- ifelse(is.na(mp), ".", ann$gene[mp])
region  <- ifelse(is.na(mp), NA, ann$region[mp])
consequence <- ifelse(snv & !is.na(m), ann$consequence[m], NA)
mlc <- ifelse(snv & !is.na(m), as.numeric(ann$mlc_score[m]), NA)
nonsyn <- consequence %in% c("missense", "stop_gained", "start_lost", "stop_lost")

# ---- metrics
detection_limit <- max(hf_min, min_alt / mt_median_depth)
low_qc <- if (detection_limit > hf_min) "low_depth" else "PASS"
if (low_qc != "PASS")
  log("warning: median chrM depth ", fmt(mt_median_depth), " gives an HF detection limit of ",
      fmt(detection_limit), " > ", hf_min, "; low-HF metrics flagged low_depth")
s <- function(x) fmt(sum(x, na.rm = TRUE))
rows <- list(
  c("heteroplasmy_count",         nrow(het),                       "count", "PASS"),
  c("heteroplasmy_low_count",     sum(low),                        "count", low_qc),
  c("heteroplasmy_high_count",    sum(!low),                       "count", "PASS"),
  c("snv_count",                  sum(snv),                        "count", "PASS"),
  c("indel_count",                sum(het$type == "INDEL"),        "count", "PASS"),
  c("H_CtoT_count",               sum(strand %in% "H" & substitution == "C>T"), "count", "PASS"),
  c("H_AtoG_count",               sum(strand %in% "H" & substitution == "A>G"), "count", "PASS"),
  c("L_AtoG_count",               sum(strand %in% "L" & substitution == "A>G"), "count", "PASS"),
  c("L_CtoT_count",               sum(strand %in% "L" & substitution == "C>T"), "count", "PASS"),
  c("transversion_count",         sum(!is.na(substitution) & !transition),      "count", "PASS"),
  c("age_accumulating_count",     sum(age_acc),                    "count", "PASS"),
  c("age_accumulating_low_count", sum(age_acc & low),              "count", low_qc),
  c("MSS",                        s(mlc),                          "score", "PASS"),
  c("MSS_low",                    s(mlc[low]),                     "score", low_qc),
  c("MSS_high",                   s(mlc[!low]),                    "score", "PASS"),
  c("MSS_missing_count",          sum(is.na(mlc)),                 "count", "PASS"),
  c("nonsynonymous_count",        sum(nonsyn),                     "count", "PASS"),
  c("MSS_nonsynonymous",          s(mlc[nonsyn]),                  "score", "PASS"),
  c("mt_median_depth",            fmt(mt_median_depth),            "x",     "PASS"),
  c("hf_detection_limit",         fmt(detection_limit),            "fraction", "PASS"),
  c("numt_excluded_count",        sum(in_range & !blocked & !filtered & !low_reads & numt_like), "count", "PASS"),
  c("numt_hf_threshold",          if (numt_file == "none") "NA" else fmt(numt_hf), "fraction", "PASS")
)
r <- do.call(rbind, rows)
results <- data.frame(sample = sample, biomarker = "heteroplasmy", metric = r[, 1],
                      value = r[, 2], unit = r[, 3], qc_flag = r[, 4], version = version)
write.table(results, out_results, sep = "\t", quote = FALSE, row.names = FALSE)

na_dot <- function(x) ifelse(is.na(x), ".", as.character(x))
out <- data.frame(het[, c("sample", "pos", "ref", "alt", "hf", "depth", "alt_reads", "filter", "type")],
                  strand = na_dot(strand), substitution = na_dot(substitution), hf_bin = hf_bin,
                  age_accumulating = ifelse(age_acc, "yes", "no"), gene = gene,
                  region = na_dot(region), consequence = na_dot(consequence),
                  mlc_score = ifelse(is.na(mlc), ".", fmt(mlc)))
write.table(out, out_variants, sep = "\t", quote = FALSE, row.names = FALSE)
