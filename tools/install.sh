#!/usr/bin/env bash
# tools/install.sh — download the Java tools used by the haplogroup and contamination
# modules into tools/ (or DEST):
#
#   HaploGrep3 3.3.2 (MIT)  tools/haplogrep3/haplogrep3.jar  + PhyloTree 17 FU tree
#   haplocheck 1.3.3 (MIT)  tools/haplocheck/haplocheck.jar
#
# Needs curl, tar, unzip and Java >= 11. Run once on a node with internet access;
# the tree is downloaded here so later runs work offline. To use copies installed
# elsewhere, set MITRIX_HAPLOGREP3 / MITRIX_HAPLOCHECK to the .jar files instead.
#
# usage: tools/install.sh [DEST]

set -euo pipefail
LOG_PREFIX=install
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

DEST=$(readlink -f "${1:-$MITRIX_HOME/tools}")
HG3_VERSION=3.3.2
HG3_URL=https://github.com/genepi/haplogrep3/releases/download/v$HG3_VERSION/haplogrep3-$HG3_VERSION-linux.tar.gz
HG3_SHA256=957165143ddd7d6ad9f34cc3af56b784fc236081e031c99bea4a574d6ad70959
HC_VERSION=1.3.3
HC_URL=https://github.com/genepi/haplocheck/releases/download/v$HC_VERSION/haplocheck.zip
HC_SHA256=c9a903925b7eddc6118aa2d06a9c86508bdc791a8e3ff49563a90b24d3a929a3

need_tools curl tar unzip sha256sum java
JAVA_MAJOR=$(java -version 2>&1 | sed -n 's/.*version "\([0-9]*\).*/\1/p' | head -1)
(( ${JAVA_MAJOR:-0} >= 11 )) || die "Java >= 11 needed (found: ${JAVA_MAJOR:-none})"
make_tmpdir

fetch() {  # url sha256 out
  log "downloading $1"
  curl -fsSL -o "$3" "$1" || die "download failed: $1"
  echo "$2  $3" | sha256sum -c --quiet - || die "checksum mismatch: $1"
}

# HaploGrep3
fetch "$HG3_URL" "$HG3_SHA256" "$TMP/hg3.tgz"
rm -rf "$DEST/haplogrep3"; mkdir -p "$DEST/haplogrep3"
tar -xzf "$TMP/hg3.tgz" -C "$DEST/haplogrep3"
log "installing tree $HAPLOGREP_TREE"
# HaploGrep3 downloads a tree into <jar dir>/trees on first use; classify the bundled
# example once so the tree is in place for offline runs
( cd "$DEST/haplogrep3" && java -Xmx2G -jar haplogrep3.jar classify --tree "$HAPLOGREP_TREE" \
    --in data/examples/evaluation-data.hsd --out "$TMP/example.txt" >/dev/null ) \
  || die "HaploGrep3 failed to install tree $HAPLOGREP_TREE"

# haplocheck
fetch "$HC_URL" "$HC_SHA256" "$TMP/hc.zip"
rm -rf "$DEST/haplocheck"; mkdir -p "$DEST/haplocheck"
unzip -q "$TMP/hc.zip" -d "$DEST/haplocheck"

log "installed: $DEST/haplogrep3/haplogrep3.jar (HaploGrep3 $HG3_VERSION, tree $HAPLOGREP_TREE)"
log "installed: $DEST/haplocheck/haplocheck.jar (haplocheck $HC_VERSION)"
