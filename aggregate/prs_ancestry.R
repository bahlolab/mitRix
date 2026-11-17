# aggregate/prs_ancestry.R — assign 1000 Genomes super-populations to projected samples.
# Args: proj.sscore (plink2 --score, PC*_AVG columns) ancestry_reference.tsv.gz out.tsv
# k-nearest neighbours (k = 20) among the 1000 Genomes founders on PC1-6; a sample gets
# the majority super-population if at least 80% of its neighbours agree, else "unassigned".

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
proj <- read.delim(args[1], check.names = FALSE)
ref  <- read.delim(gzfile(args[2]))
pcs  <- paste0("PC", 1:6)
names(proj)[1] <- "sample"
avg <- grep("^PC[0-9]+_AVG$", names(proj), value = TRUE)
names(proj)[match(avg, names(proj))] <- sub("_AVG$", "", avg)

k <- 20; min_frac <- 0.8
R <- as.matrix(ref[, pcs]); X <- as.matrix(proj[, pcs])
assign <- t(apply(X, 1, function(x) {
  if (any(is.na(x))) return(c(NA, NA))
  d <- colSums((t(R) - x)^2)
  nn <- ref$superpop[order(d)[seq_len(k)]]
  tab <- sort(table(nn), decreasing = TRUE)
  c(names(tab)[1], tab[[1]] / k)
}))
frac <- as.numeric(assign[, 2])
out <- data.frame(sample = proj$sample,
                  ancestry = ifelse(is.na(frac), "unassigned", ifelse(frac >= min_frac, assign[, 1], "unassigned")),
                  ancestry_fraction = round(frac, 2),
                  proj[, paste0("PC", 1:10)])
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
message("ancestry: ", paste(names(table(out$ancestry)), table(out$ancestry), sep = "=", collapse = " "))
