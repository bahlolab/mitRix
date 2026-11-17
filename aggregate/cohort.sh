#!/usr/bin/env bash
# aggregate/cohort.sh — one row per sample from each <sample>.tsv (mitrix cohort).
#
# Writes <outdir>/cohort.tsv (wide: sample, one column per `cohort` entry of the
# modules' module.conf, qc_flag) and <outdir>/cohort_long.tsv (all rows of all samples).
# qc_flag is PASS, or lists columns whose source row is not PASS (column:FLAG; PRS
# columns share one entry, e.g. PRS:non_EUR) and
# "contaminated" when the haplocheck contamination level is above --max-contamination
# (default 0.02, as in Gupta et al., Nature 2026), and "high_heteroplasmy" when a sample
# has more than --max-heteroplasmy (default 20) heteroplasmies.

set -euo pipefail
LOG_PREFIX=cohort
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"

SS=""; OUT=""; ALLOW=0; MAXCONT=0.02; MAXHET=20; PRS=""; CARR=""
while (( $# )); do
  case "$1" in
    --samplesheet)      SS="${2:?}"; shift 2;;
    --outdir)           OUTDIR="${2:?}"; shift 2;;
    -o)                 OUT="${2:?}"; shift 2;;
    --allow-incomplete) ALLOW=1; shift;;
    --max-contamination) MAXCONT="${2:?}"; shift 2;;
    --max-heteroplasmy) MAXHET="${2:?}"; shift 2;;
    --prs)              PRS="${2:?}"; shift 2;;
    --carriers)         CARR="${2:?}"; shift 2;;
    -h|--help)          echo "usage: mitrix cohort --samplesheet SS --outdir DIR [-o cohort.tsv] [--allow-incomplete] [--max-contamination 0.02] [--max-heteroplasmy 20] [--prs DIR/prs/prs.tsv] [--carriers DIR/carriers/carriers.tsv]" >&2; exit 0;;
    *)                  die "unknown option: $1";;
  esac
done
[[ -s "$SS" ]] || die "--samplesheet not found: $SS"
[[ -z "$PRS" || -s "$PRS" ]] || die "--prs not found: $PRS (run mitrix prs first)"
[[ -z "$CARR" || -s "$CARR" ]] || die "--carriers not found: $CARR (run mitrix carriers first)"
[[ "$MAXCONT" =~ ^[0-9]*\.?[0-9]+$ ]] || die "--max-contamination must be a fraction: $MAXCONT"
[[ "$MAXHET" =~ ^[0-9]+$ ]] || die "--max-heteroplasmy must be a count: $MAXHET"
OUT="${OUT:-$OUTDIR/cohort.tsv}"
LONG="${OUT%.tsv}_long.tsv"

CHECK=$("$HERE/check.sh" --samplesheet "$SS" --outdir "$OUTDIR" 2>/dev/null) || {
  printf '%s\n' "$CHECK" | awk -F'\t' 'NR > 1 && $2 != "OK"' >&2
  (( ALLOW )) || die "not all samples are OK (see above; mitrix check); use --allow-incomplete to continue"
  log "warning: continuing with incomplete samples (their columns are NA)"
}

# cohort columns: column<TAB>metric, in this order for known columns, then the rest
PREFERRED="mtDNA_CN deletion 7S flap_7S flap_primer heteroplasmy heteroplasmy_low age_accumulating MSS pathogenic_n pathogenic pathogenic_max_hf haplogroup macro_haplogroup contamination contamination_level"
COLS=$(for f in "$MITRIX_HOME"/biomarkers/*/module.conf; do
         for e in $(conf_get "$f" cohort); do printf '%s\t%s\n' "${e%%:*}" "${e#*:}"; done
       done | awk -F'\t' -v pref="$PREFERRED" '
         BEGIN { n = split(pref, p, " "); for (i = 1; i <= n; i++) rank[p[i]] = i }
         { if ($1 in seen) { print "duplicate cohort column " $1 > "/dev/stderr"; exit 1 }
           seen[$1] = 1; print (($1 in rank) ? rank[$1] : 100 + NR) "\t" $0 }' \
       | sort -k1,1n | cut -f2-) || die "module.conf cohort entries clash"

# PRS columns (mitrix prs): standardised scores and ancestry
if [[ -n "$PRS" ]]; then
  PRSCOLS=$(awk -F'\t' 'NR > 1 && (($3 ~ /^PRS_/ && $3 !~ /_sum$/) || $3 == "ancestry") && !seen[$3]++ { print $3 "\t" $3 }' "$PRS")
  COLS=$(printf '%s\n%s\n' "$COLS" "$PRSCOLS" | awk 'NF')
fi

# carrier columns (mitrix carriers): NAME_n, NAME_<GENE>..., NAME; flags collapse to NAME:flag
CARGROUP=""
if [[ -n "$CARR" ]]; then
  CARCOLS=$(awk -F'\t' 'NR > 1 && $3 !~ /_missing$/ && !seen[$3]++ { print $3 "\t" $3 }' "$CARR")
  COLS=$(printf '%s\n%s\n' "$COLS" "$CARCOLS" | awk 'NF')
  CARGROUP=$(awk -F'\t' 'NR > 1 && $5 == "label" && !seen[$3]++ { print $3 }' "$CARR" | paste -sd' ')
fi

# all rows of all OK samples (plus their PRS and carrier rows)
{
  echo "$RESULT_HEADER"
  printf '%s\n' "$CHECK" | awk -F'\t' 'NR > 1 { print $1 }' | while read -r s; do
    if [[ -s "$OUTDIR/$s/$s.tsv" ]]; then tail -n +2 "$OUTDIR/$s/$s.tsv"; fi
  done
  if [[ -n "$PRS" ]]; then
    printf '%s\n' "$CHECK" | awk -F'\t' 'NR == FNR { if (FNR > 1) ok[$1] = 1; next } FNR > 1 && ($1 in ok)' - "$PRS"
  fi
  if [[ -n "$CARR" ]]; then
    printf '%s\n' "$CHECK" | awk -F'\t' 'NR == FNR { if (FNR > 1) ok[$1] = 1; next } FNR > 1 && ($1 in ok)' - "$CARR"
  fi
} > "$LONG"
# user variant lists of the pathogenic module (--list NAME=FILE): NAME_n, NAME, NAME_max_hf
LISTCOLS=$(awk -F'\t' 'NR > 1 && $2 == "pathogenic" && $3 ~ /_count$/ && $3 != "pathogenic_count" {
             n = substr($3, 1, length($3) - 6); if (!seen[n]++) print n }' "$LONG" \
           | while read -r n; do printf '%s_n\t%s_count\n%s\t%s_variants\n%s_max_hf\t%s_max_hf\n' "$n" "$n" "$n" "$n" "$n" "$n"; done)
if [[ -n "$LISTCOLS" ]]; then
  COLS=$(printf '%s\n%s\n' "$COLS" "$LISTCOLS" | awk -F'\t' 'NF && !seen[$1]++')
fi
if [[ -n "$PRS" ]]; then
  n_prs=$(awk -F'\t' 'NR > 1 && $3 == "ancestry"' "$LONG" | wc -l)
  log "PRS for $n_prs of $(( $(printf '%s\n' "$CHECK" | wc -l) - 1 )) samples"
fi

# wide table; a metric ending _mean falls back to _median (runs with --stat median)
COLS="$COLS" CHECK="$CHECK" CARGROUP="$CARGROUP" awk -F'\t' -v OFS='\t' -v maxc="$MAXCONT" -v maxh="$MAXHET" '
  function grp(c,   k) { for (k = 1; k <= ng; k++) if (c == cg[k] || index(c, cg[k] "_") == 1) return cg[k]; return c }
  BEGIN {
    nc = split(ENVIRON["COLS"], cl, "\n")
    for (i = 1; i <= nc; i++) { split(cl[i], x, "\t"); name[i] = x[1]; metric[i] = x[2] }
    ng = split(ENVIRON["CARGROUP"], cg, " ")
    ns = split(ENVIRON["CHECK"], rows, "\n")
    for (i = 2; i <= ns; i++) { split(rows[i], x, "\t"); samples[++n] = x[1]; status[x[1]] = x[2] }
  }
  NR > 1 { v[$1, $3] = $4; q[$1, $3] = $6 }
  END {
    printf "sample"; for (i = 1; i <= nc; i++) printf "\t%s", name[i]; print "\tqc_flag"
    for (j = 1; j <= n; j++) {
      s = samples[j]; flags = ""
      if (status[s] != "OK") flags = status[s]
      printf "%s", s
      for (i = 1; i <= nc; i++) {
        m = metric[i]
        if (!((s, m) in v) && m ~ /_mean$/) { alt = m; sub(/_mean$/, "_median", alt); if ((s, alt) in v) m = alt }
        if ((s, m) in v) {
          printf "\t%s", v[s, m]
          if (q[s, m] != "PASS") {
            # all PRS columns share one entry per flag (PRS:non_EUR), not one per column
            f = (name[i] ~ /^PRS_/ ? "PRS" : grp(name[i])) ":" q[s, m]
            if (!((s, f) in seen)) { seen[s, f] = 1; flags = flags (flags ? ";" : "") f }
          }
        } else printf "\tNA"
      }
      if (v[s, "contamination_status"] == "YES" && v[s, "contamination_level"] != "NA" && v[s, "contamination_level"] + 0 > maxc + 0) flags = flags (flags ? ";" : "") "contaminated"
      # many heteroplasmies (e.g. a mixture of two haplotypes) whatever haplocheck says
      if (((s, "heteroplasmy_count") in v) && v[s, "heteroplasmy_count"] ~ /^[0-9]+$/ && v[s, "heteroplasmy_count"] + 0 > maxh + 0) flags = flags (flags ? ";" : "") "high_heteroplasmy"
      print "\t" (flags ? flags : "PASS")
    }
  }' "$LONG" > "$OUT"

log "$(( $(wc -l < "$OUT") - 1 )) samples -> $OUT ($LONG)"
