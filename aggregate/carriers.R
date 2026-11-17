# aggregate/carriers.R — dosages of listed variants and per-gene summaries (mitrix carriers).
# Args: gt.tsv (chr pos ref alt GT...) samples.txt version outdir NAME=list ...
# list: id chr pos ref alt gene label inheritance classification

a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) >= 5)
samples <- readLines(a[2]); version <- a[3]; out <- a[4]; specs <- a[5:length(a)]
ns <- length(samples)
g <- if (file.info(a[1])$size > 0) read.delim(a[1], header = FALSE, colClasses = "character", quote = "",
                                              na.strings = character()) else
  as.data.frame(matrix(character(), 0, 4 + ns))
names(g) <- c("chr", "pos", "ref", "alt", samples)
key <- function(x) paste(x$chr, x$pos, x$ref, x$alt, sep = ":")
g <- g[!duplicated(key(g)), ]
# GT -> ALT dosage (biallelic after splitting); missing -> NA
dose <- function(x) {
  al <- strsplit(gsub("[|/]", " ", x), " ", fixed = TRUE)
  vapply(al, function(y) if (!length(y) || any(y == ".")) NA_real_ else sum(as.integer(y) > 0), numeric(1))
}
D <- if (nrow(g)) sapply(samples, function(s) dose(g[[s]])) else matrix(numeric(), 0, ns, dimnames = list(NULL, samples))
if (is.null(dim(D))) D <- matrix(D, nrow = 1, dimnames = list(NULL, samples))
rownames(D) <- key(g)

rows <- list(); det <- list(); mt <- list()
for (sp in specs) {
  nm <- sub("=.*", "", sp)
  l <- read.delim(sub("^[^=]*=", "", sp), colClasses = "character", quote = "", na.strings = character())
  l <- l[!duplicated(key(l)), ]
  found <- key(l) %in% rownames(D)
  lf <- l[found, ]
  M <- D[key(lf), , drop = FALSE]                       # listed variants x samples
  genes <- sort(unique(l$gene))
  # severity classes per gene: sub-columns only for genes with more than one class
  if (!"severity" %in% names(l)) l$severity <- ""
  lf$severity <- l$severity[found]
  sev <- lapply(setNames(genes, genes), function(ge) {
    x <- unique(l$severity[l$gene == ge & nzchar(l$severity)])
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    if (length(x) > 1) sort(x) else character() })
  lf$severity <- gsub("[^A-Za-z0-9]+", "_", lf$severity)
  # modifiers (classification contains "modifier", e.g. APOE e4) get only their gene
  # columns: not counted in NAME_n, the NAME text or the carrier total
  mod <- if ("classification" %in% names(lf)) grepl("modifier", lf$classification, ignore.case = TRUE) else rep(FALSE, nrow(lf))
  carriers <- 0
  for (j in seq_len(ns)) {
    s <- samples[j]; d <- M[, j]
    car <- which(!is.na(d) & d > 0)
    cc <- car[!mod[car]]                                  # counted carriers (not modifiers)
    miss <- sum(is.na(d))
    q <- if (miss > 0) "missing_genotypes" else "PASS"
    if (length(cc)) carriers <- carriers + 1
    lab <- if (length(cc)) paste0(lf$gene[cc], ":", lf$label[cc], "(", ifelse(d[cc] >= 2, "hom", "het"), ")", collapse = ";") else "none"
    r <- rbind(c(paste0(nm, "_n"), length(cc), "count", q))
    for (ge in genes) {
      k <- lf$gene == ge
      r <- rbind(r, c(paste0(nm, "_", ge), sum(d[k], na.rm = TRUE), "count", q))
      for (sv in sev[[ge]]) {
        ks <- k & lf$severity == sv
        r <- rbind(r, c(paste0(nm, "_", ge, "_", sv), sum(d[ks], na.rm = TRUE), "count", q))
      }
    }
    r <- rbind(r, c(nm, lab, "label", q), c(paste0(nm, "_missing"), miss, "count", "PASS"))
    rows[[length(rows) + 1]] <- cbind(s, r)
    if (length(car)) det[[length(det) + 1]] <- data.frame(sample = s, list = nm, gene = lf$gene[car], label = lf$label[car],
      chr = lf$chr[car], pos = lf$pos[car], ref = lf$ref[car], alt = lf$alt[car], dosage = d[car],
      inheritance = lf$inheritance[car], classification = lf$classification[car], severity = lf$severity[car])
  }
  mt[[length(mt) + 1]] <- data.frame(list = nm, listed = nrow(l), found_in_vcf = sum(found),
                                     fraction = round(sum(found) / nrow(l), 4), carriers = carriers, samples = ns)
  if (sum(found) == 0) message("warning: list ", nm, ": no listed variant found in the VCF (build or chr naming?)")
}
r <- do.call(rbind, rows)
o <- data.frame(sample = r[, 1], biomarker = "carriers", metric = r[, 2], value = r[, 3], unit = r[, 4],
                qc_flag = r[, 5], version = version)
write.table(o, file.path(out, "carriers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cols <- c("sample", "list", "gene", "label", "chr", "pos", "ref", "alt", "dosage", "inheritance", "classification", "severity")
v <- if (length(det)) do.call(rbind, det) else setNames(data.frame(matrix(nrow = 0, ncol = length(cols))), cols)
write.table(v, file.path(out, "variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, mt), file.path(out, "match.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
