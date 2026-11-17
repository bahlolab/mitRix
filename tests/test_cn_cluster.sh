#!/usr/bin/env bash
# tests/test_cn_cluster.sh — unit test of aggregate/cn_cluster.R (mitrix cohort
# --cn-clusters) on simulated cohort tables: two mtDNA-CN groups -> K = 2 with 1 = lower
# group; one right-skewed group -> K = 1 (not split); a tiny third group -> rejected by
# the minimum cluster size; fixed K; samples with a mtDNA_CN qc flag left out (NA).
# Needs Rscript with mclust.

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d -t mitrix_test.XXXXXX); trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
Rscript -e 'quit(status = !requireNamespace("mclust", quietly = TRUE))' || { echo "SKIP: R package mclust not installed"; exit 0; }

Rscript - <<'EOF'
set.seed(1)
w <- function(f, cn, flag = rep("PASS", length(cn))) {
  d <- data.frame(sample = sprintf("S%03d", seq_along(cn)), mtDNA_CN = sprintf("%.2f", cn),
                  heteroplasmy = "1", qc_flag = flag)
  write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
}
two <- c(exp(rnorm(150, log(60), 0.15)), exp(rnorm(100, log(180), 0.15)))
flag <- rep("PASS", 250); flag[1] <- "mtDNA_CN:low_depth"; flag[2] <- "contaminated"
w("two.tsv", two, flag)
w("skew.tsv", exp(rgamma(300, shape = 3, rate = 6) + log(80)))     # one right-skewed group
w("tiny.tsv", c(exp(rnorm(300, log(100), 0.15)), exp(rnorm(5, log(400), 0.05))))
EOF

R="$HERE/aggregate/cn_cluster.R"
cp two.tsv two_fixed.tsv
Rscript "$R" two.tsv two_cl auto 2> two.log
Rscript "$R" skew.tsv skew_cl auto 2> skew.log
Rscript "$R" tiny.tsv tiny_cl auto 2> tiny.log
Rscript "$R" two_fixed.tsv fixed_cl 3 2> fixed.log

fail=0
check() {  # description, command
  if eval "$2"; then echo "ok   $1"; else echo "FAIL $1"; fail=1; fi
}
col() { awk -F'\t' -v c="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == c) k = i; next } { print $k }' "$1"; }
K() { awk -F'\t' 'NR == 2 { print $8 }' "$1"; }

check "two groups: K = 2"                     '[[ $(K two_cl.tsv) == 2 ]]'
check "two groups: cluster 1 = lower CN"      '[[ $(awk -F"\t" "NR==2{print \$4}" two_cl.tsv | cut -d. -f1) -lt 100 ]]'
check "two groups: S003-S150 in cluster 1"    '[[ $(col two.tsv cn_cluster | sed -n 3,150p | sort -u) == 1 ]]'
check "two groups: S151-S250 in cluster 2"    '[[ $(col two.tsv cn_cluster | sed -n 151,250p | sort -u) == 2 ]]'
check "mtDNA_CN-flagged sample is NA"         '[[ $(col two.tsv cn_cluster | sed -n 1p) == NA ]]'
check "other qc flags still clustered"        '[[ $(col two.tsv cn_cluster | sed -n 2p) == 1 ]]'
check "qc_flag stays the last column"         '[[ $(head -1 two.tsv | awk -F"\t" "{print \$NF}") == qc_flag ]]'
check "fit table written, one selected row"   '[[ $(awk -F"\t" "NR>1 && \$NF==\"TRUE\"" two_cl_fit.tsv | wc -l) == 1 ]]'
check "skewed single group: K = 1"            '[[ $(K skew_cl.tsv) == 1 ]]'
check "skewed single group: no uncertain"     '! grep -q cn_cluster:uncertain skew.tsv'
check "tiny third group rejected: K = 1"      '[[ $(K tiny_cl.tsv) == 1 ]]'
check "fixed K = 3"                           '[[ $(K fixed_cl.tsv) == 3 && $(awk -F"\t" "NR==2{print \$10}" fixed_cl.tsv) == fixed ]]'
check "uncertain samples flagged"             'grep -q "cn_cluster:uncertain" two_fixed.tsv'

if (( fail )); then for f in *.log; do echo "== $f"; cat "$f"; done; cat two_cl_fit.tsv skew_cl_fit.tsv tiny_cl_fit.tsv; exit 1; fi
echo "all cn_cluster tests passed"
