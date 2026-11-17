# lib/conf.sh — read module.conf files (docs/module-conf.md). Parsed, never sourced.

CONF_KEYS=" name version description entry language requires_core requires_modules requires_inputs tools r_packages options cohort "

# Print the value of KEY (empty if absent)
conf_get() {  # file key
  awk -v k="$2" '
    { sub(/#.*/, "") }
    /=/ {
      key = $0; sub(/=.*/, "", key); gsub(/^[ \t]+|[ \t]+$/, "", key)
      val = $0; sub(/^[^=]*=/, "", val); gsub(/^[ \t]+|[ \t]+$/, "", val)
      if (key == k) { print val; exit }
    }' "$1"
}

# Die on syntax errors, unknown keys or missing required keys
conf_validate() {  # file
  local f=$1 key bad
  [[ -s "$f" ]] || die "module.conf not found: $f"
  bad=$(awk '{ sub(/#.*/, "") } /^[ \t]*$/ { next } !/=/ { print "line " NR ": no =" }' "$f")
  [[ -z "$bad" ]] || die "$f: $bad"
  while read -r key; do
    [[ "$CONF_KEYS" == *" $key "* ]] || die "$f: unknown key '$key'"
  done < <(awk '{ sub(/#.*/, "") } /=/ { k=$0; sub(/=.*/, "", k); gsub(/[ \t]/, "", k); print k }' "$f")
  for key in name version description entry; do
    [[ -n "$(conf_get "$f" "$key")" ]] || die "$f: missing required key '$key'"
  done
  [[ "$(conf_get "$f" name)" == "$(basename "$(dirname "$f")")" ]] \
    || die "$f: name must match the module directory"
}
