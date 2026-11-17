# lib/regions.sh — region parsing, BED handling, presets.
# Needs lib/common.sh; append_range and preset lookups need load_ref_context.

PRESET_TSV=${PRESET_TSV:-$MITRIX_HOME/presets/thermo_4plex.tsv}

# Append a 1-based closed range as a BED line (0-based half-open).
# A wrap-around range (end < start) on the mtDNA contig is split at the origin.
append_range() {  # chr  start-end  outfile
  local chr=$1 range=$2 out=$3 s e
  if [[ ! "$range" =~ ^[0-9]+-[0-9]+$ ]]; then
    die "bad range '$range' (expected START-END, 1-based inclusive)"
  fi
  s=${range%-*}; e=${range#*-}
  (( s >= 1 )) || die "bad range '$range' (start must be >= 1)"
  if (( e < s )); then
    [[ "$chr" == "$MT_CHR" && -n "$MT_LEN" ]] \
      || die "bad range '$range' (end < start; wrap-around is only valid for the mtDNA contig)"
    (( e >= 1 )) || die "bad range '$range' (end must be >= 1 for a wrap-around range)"
    log "wrap-around range $range: auto-splitting into ${s}-${MT_LEN} + 1-${e}"
    awk -v c="$chr" -v s="$s" -v e="$MT_LEN" 'BEGIN{print c"\t"(s-1)"\t"e}' >> "$out"
    awk -v c="$chr" -v s=1    -v e="$e"       'BEGIN{print c"\t"0"\t"e}'      >> "$out"
    return
  fi
  awk -v c="$chr" -v s="$s" -v e="$e" 'BEGIN{print c"\t"(s-1)"\t"e}' >> "$out"
}

# Sort a BED file in place (chrom, then numeric start)
sort_bed_inplace() {
  LC_ALL=C sort -k1,1 -k2,2n "$1" -o "$1"
}

# Die if a BED file names chroms absent from the BAM header
verify_bed_chroms() {  # bed label
  local missing
  missing=$(awk 'NR==FNR{a[$1]=1; next} !($1 in a){print $1}' \
             <(cut -f1 "$CONTIGS_TSV") "$1" | sort -u || true)
  [[ -z "$missing" ]] \
    || die "$2 BED contains chrom(s) not present in BAM header: $(echo "$missing" | tr '\n' ' ')"
}

# BED file -> "chr:s-e,s-e;chr2:s-e" (1-based closed)
fmt_regions() {  # bed
  awk 'BEGIN{prev=""}
       {
         s=$2+1; e=$3;
         if ($1 != prev) {
           if (prev != "") printf ";"
           printf "%s:%d-%d", $1, s, e
           prev=$1
         } else {
           printf ",%d-%d", s, e
         }
       }
       END { printf "\n" }' "$1"
}

# Preset ranges for a region name: prints space-separated ranges, or nothing.
# Preset columns: region  build  contig(mt|chr15)  ranges(comma-separated, 1-based)
preset_ranges() {  # region [build]
  awk -F'\t' -v r="$1" -v b="${2:-any}" '
    /^#/ || NR==1 { next }
    $1==r && ($2=="any" || $2==b) { gsub(",", " ", $4); print $4; exit }' "$PRESET_TSV"
}

preset_contig() {  # region -> actual contig name (via ref context)
  local role
  role=$(awk -F'\t' -v r="$1" '/^#/ || NR==1 {next} $1==r {print $3; exit}' "$PRESET_TSV")
  case "$role" in
    mt)    printf '%s' "$MT_CHR";;
    chr15) printf '%s' "$CHR15";;
    *)     die "unknown contig role '$role' for preset region $1";;
  esac
}

# Build a sorted BED for one side: user BED, else user range on default chr,
# else preset ranges for the region name.
resolve_region() {  # out  name  default_chr  [user_range]  [user_bed]  [build]
  local out=$1 name=$2 chr=$3 urange=${4:-} ubed=${5:-} build=${6:-any} ranges r
  : > "$out"
  if [[ -n "$ubed" ]]; then
    [[ -e "$ubed" ]] || die "$name BED not found: $ubed"
    cp "$ubed" "$out"
  elif [[ -n "$urange" ]]; then
    append_range "$chr" "$urange" "$out"
  else
    ranges=$(preset_ranges "$name" "$build")
    [[ -n "$ranges" ]] || die "no preset region for $name (build $build)"
    for r in $ranges; do append_range "$chr" "$r" "$out"; done
  fi
  sort_bed_inplace "$out"
  [[ -s "$out" ]] || die "$name BED is empty"
  verify_bed_chroms "$out" "$name"
}
