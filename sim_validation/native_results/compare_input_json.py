#!/usr/bin/env python3
"""Compare the *structure* of two neoantigen *_input.json files.

Walks both documents and records every key path (list indices collapsed to
"[]") with the JSON type found there: object, array, string, int, float, bool,
null. int and float are different types (strict), so a Go writer that emits
1 where Python wrote 1.0 shows up as a difference. A key missing from some
instances of its object (but not all) is recorded as type "optional".

Prints only paths, types and counts -- never values. Map keys that are data
(mutation ids, clone ids) would leak into paths, so any object whose keys are
all digit strings, or which has more than --max-keys keys, is treated as a map
and its keys collapse to "{}".

Exit 0 when both files have the same set of (path, type) pairs, 1 when they
differ, 2 on usage/IO errors or an empty walk.
"""
import argparse
import json
import sys
from collections import Counter, defaultdict


def jtype(v):
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "float"
    if isinstance(v, str):
        return "string"
    if isinstance(v, list):
        return "array"
    if isinstance(v, dict):
        return "object"
    raise TypeError(type(v).__name__)


def walk(v, path, out, max_keys, seen):
    out[(path, jtype(v))] += 1
    if isinstance(v, dict):
        as_map = len(v) > max_keys or (v and all(k.isdigit() for k in v))
        if not as_map:
            seen[path][None] += 1          # instances of this object path
            for k in v:
                seen[path][k] += 1         # instances carrying key k
        for k, child in v.items():
            walk(child, path + ("{}" if as_map else "." + k), out, max_keys, seen)
    elif isinstance(v, list):
        for child in v:
            walk(child, path + "[]", out, max_keys, seen)


def load(p, max_keys):
    with open(p) as f:
        doc = json.load(f)
    out, seen = Counter(), defaultdict(Counter)
    walk(doc, "$", out, max_keys, seen)
    if len(out) < 2:
        sys.stderr.write("error: %s: walk found %d path(s); not a populated document\n" % (p, len(out)))
        sys.exit(2)
    # A key present in only some instances of its object is "optional":
    # record that as its own pseudo-type so always-vs-sometimes differences show.
    for path, c in seen.items():
        for k, n in c.items():
            if k is not None and n < c[None]:
                out[(path + "." + k, "optional")] += c[None] - n
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("reference", help="Python-pipeline *_input.json")
    ap.add_argument("candidate", help="Go-pipeline *_input.json")
    ap.add_argument("--max-keys", type=int, default=50,
                    help="objects with more keys than this are treated as maps (default 50)")
    a = ap.parse_args()
    ref, cand = load(a.reference, a.max_keys), load(a.candidate, a.max_keys)

    only_ref = sorted(set(ref) - set(cand))
    only_cand = sorted(set(cand) - set(ref))
    print("paths: reference=%d candidate=%d shared=%d"
          % (len(ref), len(cand), len(set(ref) & set(cand))))
    for label, keys, src in (("only in reference", only_ref, ref),
                             ("only in candidate", only_cand, cand)):
        for path, t in keys:
            print("%s: %s <%s> x%d" % (label, path, t, src[(path, t)]))
    if only_ref or only_cand:
        print("RESULT: DIFFER (%d only in reference, %d only in candidate)"
              % (len(only_ref), len(only_cand)))
        return 1
    print("RESULT: SAME STRUCTURE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
