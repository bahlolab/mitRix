# lib/stats.sh — per-base extraction and coverage statistics.
# Values are printed with %.10g (full precision for the long tables); NA if empty.

# Expand per-base depths for BED intervals: one value per overlapping position.
# Works on RLE (mosdepth) or one-line-per-base files.
extract_from_perbase() {  # bed perbase_gz
  local bed=$1 perbase=$2
  [[ -e "$perbase" ]] || die "per-base output missing: $perbase"
  zcat "$perbase" | awk -v bedfile="$bed" '
    BEGIN { while ((getline < bedfile) > 0) { c[++n]=$1; s[n]=$2+0; e[n]=$3+0 } }
    {
      for (i=1;i<=n;i++) {
        if ($1==c[i]) {
          os = ($2+0 > s[i]) ? $2+0 : s[i]
          oe = ($3+0 < e[i]) ? $3+0 : e[i]
          if (oe > os) { for (p=os; p<oe; p++) print $4+0 }
        }
      }
    }'
}

calc_mean() {
  awk 'BEGIN{tot=0; n=0} {tot+=$1; n++}
       END{if(n>0) printf "%.10g", tot/n; else print "NA"}'
}

calc_median() {
  sort -n | awk '
    { a[NR]=$1 }
    END {
      if (NR==0) { print "NA"; exit }
      mid = (NR+1)/2
      if (NR%2==1) printf "%.10g", a[int(mid)]
      else          printf "%.10g", (a[int(mid)] + a[int(mid)+1]) / 2
    }'
}

calc_stat() {  # mean|median  (values on stdin)
  if [[ "$1" == median ]]; then calc_median; else calc_mean; fi
}

# Depth over a BED from a per-base file
region_depth() {  # stat bed perbase_gz
  extract_from_perbase "$2" "$3" | calc_stat "$1"
}

# Length-weighted mean of a mosdepth regions file (col 4)
weighted_mean() {  # regions.bed.gz
  [[ -e "$1" ]] || die "mosdepth output missing: $1"
  zcat "$1" | awk 'BEGIN{tot=0; len=0}
    { l = $3 - $2; if (l<=0) next; tot += l * $4; len += l }
    END { if (len > 0) printf "%.10g", tot/len; else print "NA" }'
}

# a/b: sets RATIO_VAL (or NA) and RATIO_QC (PASS, ZERO_DEPTH or FAIL).
# Call directly, not in $(...), so both variables survive.
ratio() {  # a b
  RATIO_QC=PASS
  if [[ "$1" == NA || "$2" == NA ]]; then RATIO_QC=FAIL; RATIO_VAL=NA; return 0; fi
  if awk -v b="$2" 'BEGIN{ exit !(b+0 == 0) }'; then RATIO_QC=ZERO_DEPTH; RATIO_VAL=NA; return 0; fi
  RATIO_VAL=$(awk -v a="$1" -v b="$2" 'BEGIN{ printf "%.10g", (a+0)/(b+0) }')
}
