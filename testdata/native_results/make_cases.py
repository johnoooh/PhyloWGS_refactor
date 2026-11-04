#!/usr/bin/env python3
"""Write hand-made Go-format trees.zip + mutlist.json cases (synthetic data only).

Floats are always Python floats (json.dumps prints 1.0), counts are ints.
Run from anywhere: python3 make_cases.py
"""
import json
import os
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))


def mutlist(n_ssm, n_cnv):
    ssms = {}
    for i in range(n_ssm):
        ssms["s%d" % i] = {
            "ref_reads": [120 + (i % 30)], "total_reads": [200],
            "expected_ref_in_ref": 0.999, "expected_ref_in_variant": 0.5,
            "name": "SYN_%d" % i}
    cnvs = {}
    for i in range(n_cnv):
        cnvs["c%d" % i] = {
            "ref_reads": [4000 + 100 * i], "total_reads": [10000],
            "physical_cnvs": [{"chrom": "1", "start": 1000 * i,
                               "end": 1000 * i + 500, "major_cn": 2,
                               "minor_cn": 1, "cell_prev": [0.5]}],
            "ssms": []}
    return {"ssms": ssms, "cnvs": cnvs, "dataset_name": "synthetic"}


def tree(llh, nodes, n_ssm, n_cnv):
    """nodes: list of (idx, parent, cp, n_ssms, n_cnvs); returns an entry dict."""
    structure, pops, mutass = {}, {"0": {"cellular_prevalence": [1.0],
                                         "num_ssms": 0, "num_cnvs": 0}}, {}
    s_next = c_next = 0
    for idx, parent, cp, ns, nc in nodes:
        structure.setdefault(str(parent), []).append(idx)
        pops[str(idx)] = {"cellular_prevalence": [float(cp)],
                          "num_ssms": ns, "num_cnvs": nc}
        mutass[str(idx)] = {
            "ssms": ["s%d" % j for j in range(s_next, s_next + ns)],
            "cnvs": ["c%d" % j for j in range(c_next, c_next + nc)]}
        s_next, c_next = s_next + ns, c_next + nc
    assert (s_next, c_next) == (n_ssm, n_cnv), (s_next, c_next)
    return {"llh": float(llh), "structure": structure, "populations": pops,
            "mut_assignments": mutass}


def write_case(name, trees, n_ssm, n_cnv):
    d = os.path.join(HERE, name)
    os.makedirs(d, exist_ok=True)
    with zipfile.ZipFile(os.path.join(d, "trees.zip"), "w", zipfile.ZIP_DEFLATED) as z:
        for i, t in enumerate(trees):
            z.writestr("tree_%d_%r" % (i, t["llh"]), json.dumps(t))
    with open(os.path.join(d, "mutlist.json"), "w") as f:
        json.dump(mutlist(n_ssm, n_cnv), f, indent=2)


def linear():
    # chains 1->2->3(->4) with different SSM splits => varying clustering index
    splits = [(30, 20, 10), (40, 10, 10), (20, 20, 20), (25, 25, 10), (45, 10, 5)]
    out = []
    for i, (a, b, c) in enumerate(splits):
        out.append(tree(-1000.0 - 10 * i, [
            (1, 0, 0.99, a, 1), (2, 1, 0.7, b, 1), (3, 2, 0.4, c, 1)], 60, 3))
    return out


def identical():
    nodes = [(1, 0, 0.99, 30, 1), (2, 1, 0.6, 20, 1), (3, 1, 0.3, 10, 1)]
    return [tree(-500.0 - i, nodes, 60, 3) for i in range(3)]


def single():
    return [tree(-750.0, [(1, 0, 0.98, 40, 1), (2, 1, 0.5, 20, 2)], 60, 3)]


def superclone():
    out = []
    for i, (a, b) in enumerate([(10, 40), (11, 39), (9, 41), (20, 30)]):
        cps = [(0.95, 0.9), (0.96, 0.92), (0.94, 0.9), (0.95, 0.9)][i]
        out.append(tree(-800.0 - 5 * i, [
            (1, 0, cps[0], a, 1), (2, 1, cps[1], b, 1)], 50, 2))
    return out


def cnvonly():
    out = []
    for i in range(3):
        out.append(tree(-1200.0 - 7 * i, [
            (1, 0, 0.99, 100, 0), (2, 1, 0.5, 98, 1), (3, 1, 0.4, 0, 2),
            (4, 1, 0.2, 1, 1), (5, 1, 0.1, 1, 0)], 200, 4))
    return out


def multiprimary():
    out = []
    for i in range(4):
        out.append(tree(-900.0 - 3 * i, [
            (1, 0, 0.6, 35, 1), (2, 0, 0.4, 25, 1)], 60, 2))
    out.append(tree(-910.0, [(1, 0, 0.98, 40, 1), (2, 1, 0.4, 20, 1)], 60, 2))
    return out


if __name__ == "__main__":
    for name, fn, ns, nc in [("linear", linear, 60, 3), ("identical", identical, 60, 3),
                             ("single", single, 60, 3), ("superclone", superclone, 50, 2),
                             ("cnvonly", cnvonly, 200, 4), ("multiprimary", multiprimary, 60, 2)]:
        write_case(name, fn(), ns, nc)
