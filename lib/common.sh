# lib/common.sh — logging, dependency checks, temp dir, standard module options.
# Sourced by core steps, biomarker modules and bin/mitrix. Set LOG_PREFIX first.

MITRIX_HOME=${MITRIX_HOME:-$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)}
MITRIX_VERSION="0.7.0"
LOG_PREFIX=${LOG_PREFIX:-mitrix}

# Java tools (tools/install.sh); override with the path to each .jar
MITRIX_HAPLOGREP3=${MITRIX_HAPLOGREP3:-$MITRIX_HOME/tools/haplogrep3/haplogrep3.jar}
MITRIX_HAPLOCHECK=${MITRIX_HAPLOCHECK:-$MITRIX_HOME/tools/haplocheck/haplocheck.jar}
HAPLOGREP_TREE=phylotree-fu-rcrs@1.3   # PhyloTree 17 Forensic Update

die() { printf '%s: error: %s\n' "$LOG_PREFIX" "$*" >&2; exit 2; }
log() { printf '%s: %s\n' "$LOG_PREFIX" "$*" >&2; }

need_tools() {  # tool...
  local t
  for t in "$@"; do
    command -v "$t" >/dev/null 2>&1 || die "$t not found in PATH"
  done
}

# Create $TMP, removed on exit unless KEEP_TMP=1
make_tmpdir() {
  TMP=$(mktemp -d -t mitrix.XXXXXX)
  trap 'if (( ${KEEP_TMP:-0} )); then log "kept temp dir: $TMP"; else rm -rf "$TMP"; fi' EXIT
}

# ----- standard options (docs/module-contract.md §2) -----
SAMPLE=""
BAM=""
FASTA=""
OUTDIR=""
THREADS=4
# Read filter (SAM flag bits to exclude): unmapped, secondary, QC-fail, supplementary.
# Duplicates are kept: at chrM depth, independent fragments often share start and
# end by chance and get marked duplicate (HG00096: 16% on chrM vs 10% nuclear),
# which would bias copy number down.
FLAG=2820          # mosdepth -F
STAT="mean"
KEEP_TMP=0
# mtDNA variant calls (optional; heteroplasmy, haplogroup and contamination are
# SKIPPED without): a VCF (MitoHPC/Mutect2, mity, other) or mtDNA-Server 2's
# variants.annotated.txt
MT_VARIANTS=""
VARIANTS_SAMPLE=""  # sample name inside that file; default: --sample or the BAM/CRAM name
CALLER=""           # mitohpc|mutect2|mity|mtdnaserver2|vcf; default: detected
NUCLEAR_DEPTH=""    # nuclear sequencing depth (heteroplasmy NUMT ceiling); optional

# Parse one standard option. Sets STD_SHIFT to the number of args consumed;
# returns 1 if "$1" is not a standard option.
std_opt() {
  STD_SHIFT=2
  case "$1" in
    --sample)          SAMPLE="${2:?--sample needs a value}";;
    --bam)             BAM="${2:?--bam needs a value}";;
    --reference-fasta) FASTA="${2:?--reference-fasta needs a value}";;
    --outdir)          OUTDIR="${2:?--outdir needs a value}";;
    --threads)         THREADS="${2:?--threads needs a value}";;
    --flag)            FLAG="${2:?--flag needs a value}";;
    --stat)            STAT="${2:?--stat needs a value}";;
    --mean)            STAT="mean";   STD_SHIFT=1;;
    --median)          STAT="median"; STD_SHIFT=1;;
    --keep-tmp)        KEEP_TMP=1;    STD_SHIFT=1;;
    --mt-variants)     MT_VARIANTS="${2:?--mt-variants needs a value}";;
    --variants-sample) VARIANTS_SAMPLE="${2:?--variants-sample needs a value}";;
    --caller)          CALLER="${2:?--caller needs a value}";;
    --nuclear-depth)   NUCLEAR_DEPTH="${2:?--nuclear-depth needs a value}";;
    # deprecated in 0.4.0 (removed in a later release)
    --mt-vcf)          MT_VARIANTS="${2:?--mt-vcf needs a value}";;
    --mitohpc-name)    VARIANTS_SAMPLE="${2:?--mitohpc-name needs a value}";;
    --haplogroup-tab|--haplocheck-tab)
                       die "$1 was removed in 0.4.0: haplogroup and contamination are now computed from --mt-variants";;
    *)                 return 1;;
  esac
}

# Validate standard options; set SAMTOOLS_REF_OPT / MOSDEPTH_FASTA_OPT
check_std_args() {
  [[ -n "$SAMPLE" ]] || die "--sample is required"
  # variants-only mode: --mt-variants without --bam (coverage-based steps are skipped)
  [[ -n "$BAM" || -n "$MT_VARIANTS" ]] || die "--bam is required (or --mt-variants alone for variants-only mode)"
  [[ -n "$OUTDIR" ]] || die "--outdir is required"
  [[ -z "$BAM" || -e "$BAM" ]] || die "BAM/CRAM not found: $BAM"
  [[ "$FLAG"    =~ ^[0-9]+$ ]] || die "--flag must be a non-negative integer: $FLAG"
  [[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer: $THREADS"
  [[ "$STAT" == mean || "$STAT" == median ]] || die "--stat must be mean or median: $STAT"
  if [[ "$BAM" == *.cram && -z "$FASTA" ]]; then
    die "CRAM input requires --reference-fasta"
  fi
  [[ -z "$MT_VARIANTS" || -e "$MT_VARIANTS" ]] || die "file not found: $MT_VARIANTS"
  [[ -z "${NUCLEAR_DEPTH:-}" || "$NUCLEAR_DEPTH" =~ ^[0-9]*\.?[0-9]+$ ]] || die "--nuclear-depth must be a number: $NUCLEAR_DEPTH"
  case "$CALLER" in ""|mitohpc|mutect2|mity|mtdnaserver2|vcf) ;;
    *) die "--caller must be mitohpc, mutect2, mity, mtdnaserver2 or vcf: $CALLER";; esac
  SAMTOOLS_REF_OPT=()
  MOSDEPTH_FASTA_OPT=()
  if [[ -n "$FASTA" ]]; then
    [[ -e "$FASTA" ]] || die "reference FASTA not found: $FASTA"
    SAMTOOLS_REF_OPT=(--reference "$FASTA")
    MOSDEPTH_FASTA_OPT=(--fasta "$FASTA")
  fi
}

core_dir()   { printf '%s/%s/core' "$OUTDIR" "$SAMPLE"; }
module_dir() { printf '%s/%s/%s' "$OUTDIR" "$SAMPLE" "$1"; }

# Load core/ref_context.sh output: MT_CHR MT_LEN CHR15 CHR15_LEN BUILD, CONTIGS_TSV
load_ref_context() {
  local f k v
  f="$(core_dir)/${SAMPLE}.ref_context.tsv"
  MT_CHR=""; MT_LEN=""; CHR15=""; CHR15_LEN=""; BUILD=""
  # variants-only mode: no BAM header. A precomputed core/ (tools/from_perbase.sh) still
  # provides the context.
  [[ -z "$BAM" && ! -s "$f" ]] && return 0
  [[ -s "$f" ]] || die "ref context missing: $f (run core/ref_context.sh first)"
  while IFS=$'\t' read -r k v; do
    case "$k" in
      mt_contig)    MT_CHR="$v";;
      mt_length)    MT_LEN="$v";;
      chr15_contig) CHR15="$v";;
      chr15_length) CHR15_LEN="$v";;
      build)        BUILD="$v";;
    esac
  done < "$f"
  CONTIGS_TSV="$(core_dir)/${SAMPLE}.contigs.tsv"
}

# ----- module output (docs/module-contract.md §3) -----
RESULT_HEADER=$'sample\tbiomarker\tmetric\tvalue\tunit\tqc_flag\tversion'

write_row() {  # metric value unit qc_flag  (uses SAMPLE, MODULE, MODULE_VERSION)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SAMPLE" "$MODULE" "$1" "$2" "$3" "$4" "$MODULE_VERSION"
}

write_param() {  # key value
  printf '%s\t%s\n' "$1" "$2"
}
