#!/usr/bin/env python3
"""
Python-side glue for the C++ Metropolis-Hastings sampler (mh.o).

Optimizations vs original:
- write_data_state uses pre-computed node ancestor paths (node.path) and the
  fast cnv_node_map on each datum, instead of calling get_ancestors() and
  scanning cnv lists repeatedly in the inner loop.
- find_most_recent_cnv (module-level) uses the same O(depth) path traversal
  with an O(1) dict lookup.
- Python 3 compatible.
"""

import numpy
from numpy import *
import builtins
max = builtins.max  # restore Python builtins shadowed by numpy wildcard import
min = builtins.min
import pickle
import zipfile
import shutil

import scipy.stats as stat
from scipy.stats import beta, binom
from scipy.special import gammaln

from data import Datum
from tssb import *
from util import dirichletpdfln
from numpy.random import dirichlet

import subprocess as sp
import util2 as u2
import os


def get_c_fnames(tmp_dir):
    def _make_c_fname(name):
        return os.path.join(tmp_dir, 'c_%s.txt' % name)
    FNAME_C_TREE = _make_c_fname('tree')
    FNAME_C_DATA_STATES = _make_c_fname('data_states')
    FNAME_C_PARAMS = _make_c_fname('params')
    FNAME_C_MH_ARATIO = _make_c_fname('mh_ar')
    return (FNAME_C_TREE, FNAME_C_DATA_STATES, FNAME_C_PARAMS, FNAME_C_MH_ARATIO)


def metropolis(tssb, iters=1000, std=0.01, burnin=0,
               n_ssms=0, n_cnvs=0, fin1='', fin2='',
               rseed=1, ntps=5, tmp_dir='.'):
    """
    Run the C++ MH sampler for the current tree state.

    Writes tree and data-state files, invokes mh.o, reads back updated
    parameters.  This interface is preserved exactly; the optimizations are
    in the helper functions called here.
    """
    wts, nodes = tssb.get_mixture()

    FNAME_SSM_DATA = fin1
    FNAME_CNV_DATA = fin2
    NTPS = str(ntps)
    FNAME_C_TREE, FNAME_C_DATA_STATES, FNAME_C_PARAMS, FNAME_C_MH_ARATIO = get_c_fnames(tmp_dir)

    u2.set_node_height(tssb)
    write_tree(tssb, n_ssms, FNAME_C_TREE)
    u2.map_datum_to_node(tssb)

    # Rebuild cnv_node_maps after map_datum_to_node refreshed datum.node
    for datum in tssb.data:
        if datum.cnv:
            datum._rebuild_cnv_node_map()

    write_data_state(tssb, FNAME_C_DATA_STATES)

    MH_ITR = str(iters)
    MH_STD = str(std)
    N_SSM_DATA = str(n_ssms)
    N_CNV_DATA = str(n_cnvs)
    NNODES = str(len(nodes))
    TREE_HEIGHT = str(max([node.ht for node in nodes]) + 1)

    script_dir = os.path.dirname(os.path.realpath(__file__))
    sp.check_call([
        '%s/mh.o' % script_dir,
        MH_ITR, MH_STD, N_SSM_DATA, N_CNV_DATA, NNODES, TREE_HEIGHT,
        FNAME_SSM_DATA, FNAME_CNV_DATA,
        FNAME_C_TREE, FNAME_C_DATA_STATES, FNAME_C_PARAMS, FNAME_C_MH_ARATIO,
        NTPS,
    ])

    ar = str(loadtxt(FNAME_C_MH_ARATIO, dtype='str'))
    update_tree_params(tssb, FNAME_C_PARAMS)
    return ar


def write_tree(tssb, n_ssms, fname):
    with open(fname, 'w') as fh:
        wts, nodes = tssb.get_mixture()
        did_int_dict = {}
        for dat in tssb.data:
            if dat.id[0] == 's':
                did_int_dict[dat.id] = int(dat.id[1:])
            else:
                did_int_dict[dat.id] = n_ssms + int(dat.id[1:])

        def descend(root):
            for child in root.children():
                descend(child)
            cids = ','.join(str(child.id) for child in root.children()) or str(-1)
            dids = ','.join(str(did_int_dict[dat.id]) for dat in root.get_data()) or str(-1)
            line = '\t'.join([
                str(root.id),
                list_to_string(root.params),
                list_to_string(root.pi),
                str(len(root.children())),
                cids,
                str(len(root.get_data())),
                dids,
                str(root.ht),
            ])
            fh.write(line + '\n')

        descend(tssb.root['node'])


def list_to_string(p):
    return ','.join(str(pp) for pp in p)


def write_data_state(tssb, fname):
    """
    Write per-datum node-state file for the C++ MH sampler.

    Optimizations vs original:
    - Uses node.path (pre-cached) instead of get_ancestors() per node.
    - Uses datum._get_cnv_node_map() for O(1) CNV lookup.
    - Builds the output string with a list join instead of += concatenation.
    """
    wts, nodes = tssb.get_mixture()

    # Pre-compute ancestor sets for each node (used for membership tests below).
    # node.path is already a list; convert to set for O(1) 'in' tests.
    node_ancestor_sets = {}
    for node in nodes:
        path = node.path if node.path is not None else node.get_ancestors()
        node_ancestor_sets[node] = set(path)

    with open(fname, 'w') as fh:
        for dat in tssb.data:
            if not dat.cnv:
                continue
            if not dat.node:
                continue

            poss_n_genomes = dat.compute_n_genomes(0)
            ssm_node = dat.node.path[-1]
            cnv_map = dat._get_cnv_node_map()

            state1_parts = []
            state2_parts = []
            state3_parts = []
            state4_parts = []

            for node in nodes:
                ancestors_set = node_ancestor_sets[node]
                # Use pre-cached path for find_most_recent_cnv
                node_path = node.path if node.path is not None else node.get_ancestors()
                mr_cnv = _find_most_recent_cnv_fast(node_path, cnv_map)

                ssm_in_anc = ssm_node in ancestors_set
                has_cnv = mr_cnv is not None

                if not ssm_in_anc and not has_cnv:
                    seg = '%d,%d,%d' % (node.id, 2, 0)
                    state1_parts.append(seg)
                    state2_parts.append(seg)
                    state3_parts.append(seg)
                    state4_parts.append(seg)
                elif ssm_in_anc and not has_cnv:
                    seg = '%d,%d,%d' % (node.id, 1, 1)
                    state1_parts.append(seg)
                    state2_parts.append(seg)
                    state3_parts.append(seg)
                    state4_parts.append(seg)
                elif not ssm_in_anc and has_cnv:
                    total = int(mr_cnv[1] + mr_cnv[2])
                    seg = '%d,%d,%d' % (node.id, total, 0)
                    state1_parts.append(seg)
                    state2_parts.append(seg)
                    state3_parts.append(seg)
                    state4_parts.append(seg)
                elif ssm_in_anc and has_cnv:
                    total = int(mr_cnv[1] + mr_cnv[2])
                    seg34 = '%d,%d,%d' % (node.id, max(0, total - 1), min(1, total))
                    state3_parts.append(seg34)
                    state4_parts.append(seg34)
                    cnv_anc_path = (
                        mr_cnv[0].node.path
                        if mr_cnv[0].node.path is not None
                        else mr_cnv[0].node.get_ancestors()
                    )
                    if ssm_node in set(cnv_anc_path):
                        state1_parts.append('%d,%d,%d' % (node.id, int(mr_cnv[1]), int(mr_cnv[2])))
                        state2_parts.append('%d,%d,%d' % (node.id, int(mr_cnv[2]), int(mr_cnv[1])))
                    else:
                        seg12 = '%d,%d,%d' % (node.id, max(0, total - 1), min(1, total))  # total already int from above
                        state1_parts.append(seg12)
                        state2_parts.append(seg12)
                else:
                    print("PANIC")

            # Apply poss_n_genomes corrections (mirrors original logic)
            if poss_n_genomes[0][1] == 0:
                state1_parts = state2_parts[:]
            elif poss_n_genomes[1][1] == 0:
                state2_parts = state1_parts[:]
            if len(poss_n_genomes) == 2:
                state3_parts = state1_parts[:]
                state4_parts = state2_parts[:]

            fh.write(
                dat.id[1:] + '\t'
                + ';'.join(state1_parts) + '\t'
                + ';'.join(state2_parts) + '\t'
                + ';'.join(state3_parts) + '\t'
                + ';'.join(state4_parts) + '\t\n'
            )


def _find_most_recent_cnv_fast(node_path, cnv_map):
    """
    Walk ancestors from leaf to root, return first CNV found via dict lookup.
    node_path: list of nodes from root to current node (inclusive).
    cnv_map:   dict mapping node → (cnv_datum, maternal_cn, paternal_cn).
    """
    if not cnv_map:
        return None
    for n in reversed(node_path):
        entry = cnv_map.get(n)
        if entry is not None:
            return entry
    return None


def find_most_recent_cnv(dat, nd):
    """Module-level helper kept for compatibility with write_data_state callers."""
    cnv_map = dat._get_cnv_node_map()
    node_path = nd.path if nd.path is not None else nd.get_ancestors()
    return _find_most_recent_cnv_fast(node_path, cnv_map)


def update_tree_params(tssb, fname):
    wts, nodes = tssb.get_mixture()
    ndict = {node.id: node for node in nodes}
    with open(fname) as fh:
        params = [line.split() for line in fh.readlines()]
    for p in params:
        ndict[int(p[0])].params = string_to_list(p[1])
        ndict[int(p[0])].pi = string_to_list(p[2])


def string_to_list(p):
    p = p.strip(',')
    return array([float(pp) for pp in p.split(',')])


def sample_cons_params(tssb, tp):
    def descend(root, tp):
        if root.parent() is None:
            root.params1[tp] = 1
            root.pi1[tp] = root.params1[tp] * rand(1)
        r = root.params1[tp] - root.pi1[tp]
        p = rand(len(root.children()))
        if len(root.children()) > 0:
            p = r * p * 1. / sum(p)
        index = 0
        for child in root.children():
            child.params1[tp] = p[index]
            child.pi1[tp] = child.params1[tp] * (rand(1) ** (len(child.children()) > 0))
            index += 1
        for child in root.children():
            descend(child, tp)
    descend(tssb.root['node'], tp)


def update_params(tssb, tp):
    def descend(root, tp):
        for child in root.children():
            descend(child, tp)
        root.params[tp] = root.params1[tp]
        root.pi[tp] = root.pi1[tp]
    descend(tssb.root['node'], tp)
