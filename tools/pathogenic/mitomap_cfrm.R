# tools/pathogenic/mitomap_cfrm.R — MITOMAP "Cfrm" (confirmed) variants as VCF-style
# pos/ref/alt (used by tools/make_pathogenic.sh).
#
# Args: coding_control.csv rna.csv rCRS.fa out.tsv
#   coding_control.csv  MITOMAP Foswiki export "Mutations Coding/Control" (Allele, Position,
#                       NucleotideChange, Locus, Disease, Status, ...)
#   rna.csv             MITOMAP Foswiki export "Mutations rRNA/tRNA" (Allele, Position, Locus, ...)
# Writes: pos ref alt gene disease mitomap_status (all rows whose Status starts with Cfrm).
#
# Allele parsing follows add_VAR_MitoMAP.R (Bahlo lab, ataxia project): SNVs/MNVs verbatim;
# insertions keep the REF base; deletions are anchored on the rCRS base upstream
# (VCF style); unparseable composite notations are dropped with a warning.
# MITOMAP data: CC BY 3.0 (https://www.mitomap.org), cite MITOMAP (PMID 25489354).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
fa <- readLines(args[3]); rcrs <- toupper(substr(paste(fa[!startsWith(fa, ">")], collapse = ""), 1, 16569))
stopifnot(nchar(rcrs) == 16569)
base_at <- function(p) { if (p < 1L) p <- 16569L; if (p > 16569L) p <- 1L; substr(rcrs, p, p) }
seq_at <- function(s, e) substr(rcrs, s, e)
v3 <- function(pos, ref, alt) c(pos, ref, alt)
NAv <- c(NA, NA, NA)

var_coding <- function(allele, pos, nuc) {
  pos <- suppressWarnings(as.integer(pos))
  if (is.na(pos) || pos < 1L || pos > 16569L) return(NAv)
  body <- sub("^m\\.", "", allele)
  if (grepl("del", body, ignore.case = TRUE)) {
    pre <- sub("del.*$", "", body)
    end <- if (grepl("_", pre)) as.integer(sub("^.*_", "", pre)) else pos
    np <- sub("-del$", "", nuc, ignore.case = TRUE)
    d <- if (grepl("^[ACGTN]+$", np) && nchar(np) == end - pos + 1L) np else seq_at(pos, end)
    a <- base_at(pos - 1L)
    return(v3(pos - 1L, paste0(a, d), a))
  }
  p <- strsplit(nuc, "-", fixed = TRUE)[[1]]
  if (length(p) == 2L && nzchar(p[1]) && nzchar(p[2])) return(v3(pos, p[1], p[2]))
  NAv
}

var_rna <- function(allele, pos) {
  pos <- suppressWarnings(as.integer(pos))
  if (is.na(pos) || pos < 1L || pos > 16569L) return(NAv)
  a <- gsub("\\s+", "", allele)
  rx <- function(p) regmatches(a, regexec(p, a, ignore.case = TRUE))[[1]]
  m <- rx("^([ACGT])([0-9]+):$")                      # C12325: -> 1-base deletion
  if (length(m) == 3L) { an <- base_at(pos - 1L); return(v3(pos - 1L, paste0(an, toupper(m[2])), an)) }
  if (grepl("[()/]|\\+", a)) return(NAv)               # composite notations
  m <- rx("^([ACGT])([0-9]+)ins([ACGT]+)$")
  if (length(m) == 4L) { r <- toupper(m[2]); return(v3(pos, r, paste0(r, toupper(m[4])))) }
  m <- rx("^([ACGT])?([0-9]+)del$")
  if (length(m) == 3L) { r <- toupper(m[2]); if (!nzchar(r)) r <- base_at(pos); an <- base_at(pos - 1L)
                         return(v3(pos - 1L, paste0(an, r), an)) }
  m <- rx("^([0-9]+)_([0-9]+)del([ACGT]*)$")
  if (length(m) == 4L) { s <- as.integer(m[2]); e <- as.integer(m[3]); b <- toupper(m[4])
                         if (!nzchar(b) || nchar(b) != e - s + 1L) b <- seq_at(s, e)
                         an <- base_at(s - 1L); return(v3(s - 1L, paste0(an, b), an)) }
  m <- rx("^([ACGT])([0-9]+)([ACGT]+)$")
  if (length(m) == 4L) return(v3(pos, toupper(m[2]), toupper(m[4])))
  NAv
}

rd <- function(f) read.csv(f, check.names = FALSE, colClasses = "character", encoding = "UTF-8")
co <- rd(args[1]); rn <- rd(args[2])
co <- co[grepl("^Cfrm", co$Status), ]; rn <- rn[grepl("^Cfrm", rn$Status), ]
vc <- t(mapply(var_coding, co$Allele, co$Position, co$NucleotideChange, USE.NAMES = FALSE))
vr <- t(mapply(var_rna, rn$Allele, rn$Position, USE.NAMES = FALSE))
clean <- function(x) gsub("[\t\r\n]+", " ", trimws(x))
out <- rbind(
  data.frame(pos = vc[, 1], ref = vc[, 2], alt = vc[, 3], gene = co$Locus, disease = clean(co$Disease), mitomap_status = co$Status),
  data.frame(pos = vr[, 1], ref = vr[, 2], alt = vr[, 3], gene = rn$Locus, disease = clean(rn$Disease), mitomap_status = rn$Status))
bad <- is.na(out$pos)
if (any(bad)) message("mitomap: skipped ", sum(bad), " unparseable allele(s)")
out <- out[!bad, ]
ok <- mapply(function(p, r) seq_at(as.integer(p), as.integer(p) + nchar(r) - 1L) == r, out$pos, out$ref)
if (any(!ok)) { message("mitomap: dropped ", sum(!ok), " allele(s) whose REF does not match rCRS: ",
                        paste(out$pos[!ok], out$ref[!ok], collapse = "; ")); out <- out[ok, ] }
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
message("mitomap: ", nrow(out), " Cfrm variants (coding/control ", nrow(co), ", rRNA/tRNA ", nrow(rn), ")")
