#!/usr/bin/env bash
# Regenerate expected/ golden files with the Python 2 reference (write_results_go.py)
# in the phylowgs container. The ONLY way expected files may change.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
IMG=ghcr.io/mskcc-omics-workflows/phylowgs:v1.5-msk
cd "$HERE"
python3 make_cases.py

ref() { docker run -i --rm --platform linux/amd64 -v "$HERE:$HERE" -w "$HERE" "$IMG" "$@"; }

ref python2 -c 'import scipy,numpy;print scipy.__version__, numpy.__version__' > SCIPY_VERSION

# run_one <case> <flagset> <flags...>
run_one() {
  local case=$1 fs=$2; shift 2
  local exp="$HERE/$case/$fs/expected" tmp="$HERE/.tmp_${case}_${fs}"
  rm -rf "$exp" "$tmp"; mkdir -p "$exp" "$tmp"
  if ref python2 write_results_go.py "$@" --include-ssm-names "$case" \
        "$case/trees.zip" "$case/mutlist.json" \
        "$tmp/summ.json.gz" "$tmp/muts.json.gz" "$tmp/mutass.zip" \
        >"$tmp/stdout" 2>"$tmp/stderr"; then
    ref python2 dump_expected.py "$tmp" "$exp"
  else
    # Controller ruling 2: error.txt holds the bare exception message.
    tail -n 1 "$tmp/stderr" | sed -E 's/^[A-Za-z_][A-Za-z0-9_.]*: //' > "$exp/error.txt"
  fi
  rm -rf "$tmp"
}

PROD=(--max-multiprimary 1.0 --include-multiprimary)
for c in linear identical single superclone cnvonly; do
  run_one $c default
  run_one $c prod "${PROD[@]}"
done
run_one multiprimary default
run_one multiprimary prod "${PROD[@]}"
run_one sim default
run_one sim prod "${PROD[@]}"
run_one sim keepsc --keep-superclones
run_one sim mincount --min-ssms 5

# nonames: the reference ignores --include-ssm-names, so expected = default minus ssms/*/name.
rm -rf sim/nonames; mkdir -p sim/nonames/expected
cp sim/default/expected/summ.json sim/default/expected/mutass.json sim/nonames/expected/
ref python2 - <<'PY'
import json
d = "sim"
m = json.load(open(d + "/default/expected/muts.json"))
for s in m["ssms"].values():
    s.pop("name", None)
open(d + "/nonames/expected/muts.json", "w").write(
    json.dumps(m, sort_keys=True, indent=1, separators=(',', ': ')) + '\n')
PY
