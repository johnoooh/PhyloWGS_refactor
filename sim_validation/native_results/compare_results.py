#!/usr/bin/env python3
"""Compare decoded write-results outputs (summ/muts/mutass .json) of Go vs Python.

usage: compare_results.py [--loose-int-float] LABEL GO_DIR PY_DIR

Rules mirror diffPythonJSON (results_compare_test.go): same keys, same number
class (unless --loose-int-float), exact values, except tree_densities (rel
1e-9) and mut_assignments ssms/cnvs lists (compared sorted).

Ruling 8: tree_densities differences are `densities_undefined` (counted, not a
failure) when the reference trees have <= 2 distinct
(clustering_index, branching/(branching+linearity+1e-4)) points.

Patient-data safe: prints only LABEL, file:top-level-key names and counts.
Never prints values, or keys below the top level. Exit 1 on any failure.
stdlib only.
"""
import json
import math
import os
import sys

FILES = ("summ.json", "muts.json", "mutass.json")


def is_float_lit(s):
    return any(c in s for c in ".eE") or s in ("NaN", "Infinity", "-Infinity")


def load(path):
    with open(path) as f:
        # keep numeric literals verbatim so int/float class survives
        return json.load(f, parse_float=lambda s: ("F", s), parse_int=lambda s: ("I", s))


def kind(v):
    if isinstance(v, tuple):
        return "number"
    if v is None:
        return "null"
    return type(v).__name__


class Differ:
    def __init__(self, loose):
        self.loose = loose
        self.count = 0
        self.dens_only = True  # all diffs seen are under tree_densities

    def add(self, path):
        self.count += 1
        if not (path and path[0] == "tree_densities"):
            self.dens_only = False

    def number(self, path, g, w):
        (gc, gs), (wc, ws) = g, w
        if gc != wc:
            if self.loose and float(gs) == float(ws):
                return
            return self.add(path)
        if gc == "I":
            if gs != ws:
                self.add(path)
            return
        gf, wf = float(gs), float(ws)
        if gf == wf or (math.isnan(gf) and math.isnan(wf)):
            return
        if path and path[0] == "tree_densities" and \
                abs(gf - wf) <= 1e-9 * max(abs(gf), abs(wf)):
            return
        self.add(path)

    def walk(self, path, g, w):
        if kind(g) != kind(w):
            return self.add(path)
        if isinstance(w, dict):
            for k in sorted(set(g) | set(w)):
                if k not in g or k not in w:
                    self.add(path + [k])
                else:
                    self.walk(path + [k], g[k], w[k])
        elif isinstance(w, list):
            n = len(path)
            if n >= 3 and path[n - 3] == "mut_assignments" and path[n - 1] in ("ssms", "cnvs"):
                g, w = sorted(g, key=str), sorted(w, key=str)
            if len(g) != len(w):
                return self.add(path)
            for i, (a, b) in enumerate(zip(g, w)):
                self.walk(path + [str(i)], a, b)
        elif isinstance(w, tuple):
            self.number(path, g, w)
        elif g != w:
            self.add(path)


def distinct_points(summ):
    pts = set()
    for t in summ.get("trees", {}).values():
        c = float(t["clustering_index"][1])
        b = float(t["branching_index"][1])
        l = float(t["linearity_index"][1])
        pts.add((c, b / (b + l + 1e-4)))
    return len(pts)


def compare(label, go_dir, py_dir, loose):
    fails = {}
    dens = False
    for fn in FILES:
        gp, pp = os.path.join(go_dir, fn), os.path.join(py_dir, fn)
        if not (os.path.isfile(gp) and os.path.isfile(pp)):
            fails[fn] = fails.get(fn, 0) + 1
            continue
        g, w = load(gp), load(pp)
        top = sorted(set(g) | set(w)) if isinstance(w, dict) and isinstance(g, dict) else [None]
        for k in top:
            d = Differ(loose)
            if k is None:
                d.walk([], g, w)
            elif k not in g or k not in w:
                d.add([k])
            else:
                d.walk([k], g[k], w[k])
            if d.count == 0:
                continue
            if fn == "summ.json" and k == "tree_densities" and d.dens_only \
                    and distinct_points(w) <= 2:
                dens = True
                continue
            name = fn[:-5] + (":" + k if k is not None and fn != "mutass.json" else "")
            fails[name] = fails.get(name, 0) + d.count
    if fails:
        print("%s FAIL %s" % (label, " ".join("%s=%d" % kv for kv in sorted(fails.items()))))
        return "fail"
    if dens:
        print("%s densities_undefined" % label)
        return "densities_undefined"
    print("%s ok" % label)
    return "ok"


def main(argv):
    loose = "--loose-int-float" in argv
    args = [a for a in argv if a != "--loose-int-float"]
    if len(args) != 3:
        sys.exit(__doc__)
    r = compare(args[0], args[1], args[2], loose)
    sys.exit(1 if r == "fail" else 0)


if __name__ == "__main__":
    main(sys.argv[1:])
