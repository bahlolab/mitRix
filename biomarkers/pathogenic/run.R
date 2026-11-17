# biomarkers/pathogenic/run.R — match a sample's variant calls to pathogenic-variant lists.
# Args: variants.tsv coverage.bed.gz hf_min min_alt exclude_filters sample version
#       out.tsv out_variants.tsv NAME=list.tsv [NAME=list.tsv ...]
# list.tsv: pos ref alt gene disease source (normalised against rCRS by run.sh)

a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) >= 10)
hf_min <- as.numeric(a[3]); min_alt <- as.integer(a[4])
exclude <- strsplit(a[5], ",", fixed = TRUE)[[1]]
sample <- a[6]; version <- a[7]; out_file <- a[8]; var_file <- a[9]
specs <- a[10:length(a)]
fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.4g", x))
chr <- c(ref = "character", alt = "character", filter = "character")

v <- read.delim(a[1], colClasses = chr, na.strings = character())
v$hf <- as.numeric(v$hf); v$depth <- as.numeric(v$depth)
v$alt_reads <- suppressWarnings(as.numeric(v$alt_reads))
v$alt_reads <- ifelse(is.na(v$alt_reads), round(v$hf * v$depth), v$alt_reads)   # callers without AD
excluded <- vapply(strsplit(v$filter, ";", fixed = TRUE), function(f) any(f %in% exclude), logical(1))
ok <- !excluded & !is.na(v$hf) & v$hf >= hf_min & !is.na(v$alt_reads) & v$alt_reads >= min_alt
v <- v[ok, ]
key <- function(x) paste(x$pos, x$ref, x$alt)

cov <- read.delim(gzfile(a[2]), header = FALSE, colClasses = c("character", "integer", "integer", "numeric"))
depth_at <- function(p) { i <- findInterval(p - 1, cov$V2); ifelse(i > 0 & p <= cov$V3[pmax(i, 1)], cov$V4[pmax(i, 1)], 0) }
min_depth <- min_alt / hf_min

rows <- list(); det <- list()
for (s in specs) {
  nm <- sub("=.*", "", s); f <- sub("^[^=]*=", "", s)
  l <- read.delim(f, colClasses = c(ref = "character", alt = "character", gene = "character",
                                    disease = "character", source = "character"),
                  quote = "", na.strings = character())
  l <- l[!duplicated(key(l)), ]
  hit <- v[key(v) %in% key(l), ]
  hit <- hit[order(hit$pos), ]
  li <- match(key(hit), key(l))
  # positions where an absent call is uninformative: depth too low to see HF-MIN with MIN-ALT reads
  absent <- l[!(key(l) %in% key(hit)), ]
  low <- sum(depth_at(absent$pos) < min_depth)
  lab <- if (nrow(hit)) paste0("m.", hit$pos, hit$ref, ">", hit$alt, "(", fmt(hit$hf), ")", collapse = ";") else "none"
  q <- if (low > 0) "low_depth" else "PASS"
  rows[[length(rows) + 1]] <- rbind(
    c(paste0(nm, "_count"), nrow(hit), "count", q),
    c(paste0(nm, "_variants"), lab, "label", "PASS"),
    c(paste0(nm, "_max_hf"), if (nrow(hit)) fmt(max(hit$hf)) else "0", "fraction", "PASS"),
    c(paste0(nm, "_low_depth_sites"), low, "count", "PASS"))
  if (nrow(hit)) det[[length(det) + 1]] <- data.frame(
    sample = sample, list = nm, variant = paste0("m.", hit$pos, hit$ref, ">", hit$alt),
    pos = hit$pos, ref = hit$ref, alt = hit$alt, hf = hit$hf, depth = hit$depth, alt_reads = hit$alt_reads,
    filter = hit$filter, plasmy = ifelse(hit$hf >= 0.95, "homoplasmic", "heteroplasmic"),
    gene = l$gene[li], disease = l$disease[li], source = l$source[li])
  if (low > 0) message("pathogenic: list ", nm, ": ", low, " listed position(s) with depth < ", fmt(min_depth))
  message("pathogenic: list ", nm, ": ", nrow(hit), " detected", if (nrow(hit)) paste0(" (", lab, ")") else "")
}
r <- do.call(rbind, rows)
out <- data.frame(sample = sample, biomarker = "pathogenic", metric = r[, 1], value = r[, 2],
                  unit = r[, 3], qc_flag = r[, 4], version = version)
write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cols <- c("sample", "list", "variant", "pos", "ref", "alt", "hf", "depth", "alt_reads", "filter", "plasmy", "gene", "disease", "source")
d <- if (length(det)) do.call(rbind, det) else setNames(data.frame(matrix(nrow = 0, ncol = length(cols))), cols)
write.table(d, var_file, sep = "\t", quote = FALSE, row.names = FALSE)
