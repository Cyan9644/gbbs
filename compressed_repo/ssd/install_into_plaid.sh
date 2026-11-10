#!/usr/bin/env bash
# Copies the GBBS rewrite into a PLAID checkout so the SSD backend can be
# built next to EM-graph (PLAID's graph/ submodule). Run on baldr:
#
#   gbbs/compressed_repo/ssd/install_into_plaid.sh ~/PLAID
#   cd ~/PLAID && bazel build -c opt //gbbs/ssd:ssd_bench //graph:graph_main
#
# What it does:
#   1. compressed_repo/  ->  PLAID/gbbs/   (its //gbbs/... labels then resolve;
#      unit_tests/ is left out because it needs gtest and gbbs' build macros)
#   2. the unchanged GBBS BFS.h and PageRank.h  ->  PLAID/gbbs/ssd/apps/
#   3. makes EM-graph's Bazel targets visible to other packages (one line at
#      the top of graph/BUILD.bazel), so //gbbs can depend on //graph.
# Safe to re-run: it overwrites PLAID/gbbs and adds the visibility line once.
set -euo pipefail

if [[ $# -ne 1 || ! -f "$1/WORKSPACE.bazel" || ! -d "$1/graph" ]]; then
  echo "usage: $0 <PLAID root, with the graph/ submodule checked out>" >&2
  exit 1
fi
plaid="$(cd "$1" && pwd)"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/.." && pwd)"          # compressed_repo
gbbs_root="$(cd "$repo/.." && pwd)"     # the gbbs repository

marker="$plaid/gbbs/.installed_by_install_into_plaid"
if [[ -e "$plaid/gbbs" && ! -e "$marker" ]]; then
  echo "$plaid/gbbs exists and was not made by this script; not touching it" >&2
  exit 1
fi
rm -rf "$plaid/gbbs"
mkdir -p "$plaid/gbbs"
(cd "$repo" && tar --exclude=./unit_tests -cf - .) | (cd "$plaid/gbbs" && tar -xf -)

touch "$marker"
mkdir -p "$plaid/gbbs/ssd/apps"
cp "$gbbs_root/benchmarks/BFS/NonDeterministicBFS/BFS.h" "$plaid/gbbs/ssd/apps/"
cp "$gbbs_root/benchmarks/PageRank/PageRank.h" "$plaid/gbbs/ssd/apps/"

build="$plaid/graph/BUILD.bazel"
if ! grep -q 'default_visibility' "$build"; then
  { echo 'package(default_visibility = ["//visibility:public"])'; echo; cat "$build"; } > "$build.tmp"
  mv "$build.tmp" "$build"
  echo "note: made //graph targets public (local edit to graph/BUILD.bazel)"
fi

echo "installed into $plaid/gbbs"
echo "next: cd $plaid && bazel build -c opt //gbbs/ssd:ssd_bench //graph:graph_main"
