#!/usr/bin/env python3
"""
I/O utilities, tree serialisation, and tree-metadata helpers for PhyloWGS.

Optimizations vs original:
- Python 3 (pickle, queue, print()).
- remove_empty_nodes, TreeReader, TreeWriter, StateManager, BackupManager
  unchanged in semantics but ported to Python 3 syntax.
- logmsg writes to stdout by default (unchanged).
"""

import numpy
from numpy import *
import pickle
import zipfile
import shutil
import os
import sys
import csv

import scipy.stats as stat
from scipy.stats import beta, binom
from scipy.special import gammaln
from math import exp, log
from datetime import datetime
from collections import defaultdict

csv.field_size_limit(2147483647)

from data import Datum
from tssb import *


# ---------------------------------------------------------------------------
# Numeric helpers
# ---------------------------------------------------------------------------

def log_factorial(n):
    return gammaln(n + 1)


def log_bin_coeff(n, k):
    return log_factorial(n) - log_factorial(k) - log_factorial(n - k)


def log_binomial_likelihood(x, n, mu):
    return x * log(mu) + (n - x) * log(1 - mu)


def log_beta(a, b):
    return gammaln(a) + gammaln(b) - gammaln(a + b)


def logsumexp(X, axis=None):
    import scipy.special
    return scipy.special.logsumexp(X, axis=axis)


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_physical_cnvs(pcnvs):
    physical_cnvs = []
    for physical_cnv in pcnvs.split(';'):
        fields = physical_cnv.split(',')
        cnv = dict([F.split('=', 1) for F in fields])
        for key in ('start', 'end', 'major_cn', 'minor_cn'):
            cnv[key] = int(cnv[key])
        cnv['cell_prev'] = [float(C) for C in cnv['cell_prev'].split('|')]
        physical_cnvs.append(cnv)
    return physical_cnvs


def load_data(fname1, fname2):
    """Load SSM and CNV data from tab-separated files."""
    # --- SSMs ---
    reader = csv.DictReader(open(fname1, 'r', newline=''), delimiter='\t')
    data = {}
    for row in reader:
        name = row['gene']
        id_ = row['id']
        a = [int(x) for x in row['a'].split(',')]
        d = [int(x) for x in row['d'].split(',')]
        mu_r = mu_v = 0
        if 'mu_r' in row:
            mu_r = float(row['mu_r'])
            mu_v = float(row['mu_v'])
        data[id_] = Datum(name, id_, a, d, mu_r, mu_v)

    n_ssms = len(data)

    # --- CNVs ---
    reader = csv.DictReader(open(fname2, 'r', newline=''), delimiter='\t')
    cnv_logical_physical_mapping = {}
    for row in reader:
        name = row['cnv']
        id_ = row['cnv']
        cnv_logical_physical_mapping[id_] = parse_physical_cnvs(row['physical_cnvs'])
        a = [int(x) for x in row['a'].split(',')]
        d = [int(x) for x in row['d'].split(',')]
        data[id_] = Datum(name, id_, a, d, 0.999, 0.5)

        ssms_field = row.get('ssms', '')
        if ssms_field:
            for ssm in ssms_field.split(';'):
                tok = ssm.split(',')
                data[tok[0]].cnv.append((data[id_], int(tok[1]), int(tok[2])))

    n_cnvs = len(data) - n_ssms

    return (
        [data[key] for key in data.keys()],
        n_ssms,
        n_cnvs,
        cnv_logical_physical_mapping,
    )


# ---------------------------------------------------------------------------
# Tree-metadata helpers (called once per MCMC iteration, not per datum)
# ---------------------------------------------------------------------------

def set_node_height(tssb):
    tssb.root['node'].ht = 0

    def descend(root, ht):
        for child in root.children():
            child.ht = ht
            descend(child, ht + 1)

    descend(tssb.root['node'], 1)


def set_path_from_root_to_node(tssb):
    """
    Pre-compute and cache each node's ancestor path.
    We call node._invalidate_ancestors_cache() first so that get_ancestors()
    re-computes fresh paths after any structural change.
    """
    # Invalidate stale caches from root downward.
    tssb.root['node']._invalidate_ancestors_cache()
    nodes = tssb.get_nodes()
    for node in nodes:
        node.path = node.get_ancestors()


def map_datum_to_node(tssb):
    nodes = tssb.get_nodes()
    for node in nodes:
        for datum in node.get_data():
            datum.node = node


# ---------------------------------------------------------------------------
# Empty-node removal (used in post-processing)
# ---------------------------------------------------------------------------

def remove_empty_nodes(root, parent=None):
    for child in list(root['children']):
        remove_empty_nodes(child, root)
    if root['node'].get_data() == []:
        if root['children'] == []:
            if parent is not None:
                ind = parent['children'].index(root)
                parent['children'].remove(root)
                root['node'].kill()
                parent['sticks'] = delete(parent['sticks'], ind, 0)
            return
        else:
            if parent is not None:
                parent_ = root['node'].parent()
                ind = parent['children'].index(root)
                for i, child in enumerate(list(root['children'])):
                    parent['children'].append(child)
                    toappend = zeros((1, 1))
                    toappend[0] = root['sticks'][i]
                    parent['sticks'] = append(parent['sticks'], toappend, 0)
                    root['children'].remove(child)
                for child in list(root['node'].children()):
                    child._parent = parent_
                    parent_.add_child(child)
                    child._invalidate_ancestors_cache()
                    root['node'].remove_child(child)
                parent['children'].remove(root)
                parent['sticks'] = delete(parent['sticks'], ind, 0)
                root['node'].kill()


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------

def rm_safely(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno == 2:
            pass
        else:
            raise


# ---------------------------------------------------------------------------
# Backup / state / tree I/O classes
# ---------------------------------------------------------------------------

class CorruptZipFileError(Exception):
    pass


class BackupManager(object):
    def __init__(self, filenames):
        self._filenames = filenames
        self._backup_filenames = [os.path.realpath(fn) + '.backup' for fn in self._filenames]

    def save_backup(self):
        for fn, backup_fn in zip(self._filenames, self._backup_filenames):
            shutil.copy2(fn, backup_fn)

    def restore_backup(self):
        for fn, backup_fn in zip(self._filenames, self._backup_filenames):
            shutil.copy2(backup_fn, fn)

    def remove_backup(self):
        for backup_fn in self._backup_filenames:
            try:
                os.remove(backup_fn)
            except OSError:
                pass


class StateManager(object):
    default_last_state_fn = 'state.last.pickle'
    default_initial_state_fn = 'state.initial.pickle'

    def __init__(self):
        self._initial_state_fn = StateManager.default_initial_state_fn
        self._last_state_fn = StateManager.default_last_state_fn

    def _write_state(self, state, state_fn):
        with open(state_fn, 'wb') as state_file:
            pickle.dump(state, state_file, protocol=pickle.HIGHEST_PROTOCOL)

    def write_state(self, state):
        self._write_state(state, self._last_state_fn)

    def load_state(self):
        with open(self._last_state_fn, 'rb') as state_file:
            return pickle.load(state_file)

    def load_initial_state(self):
        with open(self._initial_state_fn, 'rb') as state_file:
            return pickle.load(state_file)

    def write_initial_state(self, state):
        self._write_state(state, self._initial_state_fn)

    def delete_state_file(self):
        rm_safely(self._last_state_fn)

    def state_exists(self):
        return os.path.isfile(self._last_state_fn)


class TreeWriter(object):
    default_archive_fn = 'trees.zip'

    def __init__(self, resume_run=False):
        self._archive_fn = TreeWriter.default_archive_fn
        if resume_run:
            self._ensure_archive_is_valid()
        else:
            rm_safely(self._archive_fn)

    def add_extra_file(self, filename, data):
        self._open_archive()
        self._archive.writestr(filename, data)
        self._close_archive()

    def _ensure_archive_is_valid(self):
        with zipfile.ZipFile(self._archive_fn) as zipf:
            if zipf.testzip() is not None:
                raise CorruptZipFileError('Corrupt zip file: %s' % self._archive_fn)

    def _open_archive(self):
        self._archive = zipfile.ZipFile(
            self._archive_fn, 'a', compression=zipfile.ZIP_DEFLATED, allowZip64=True
        )

    def _close_archive(self):
        self._archive.close()

    def write_trees(self, serialized_trees):
        self._open_archive()
        for serialized_tree, idx, llh in serialized_trees:
            is_burnin = idx < 0
            prefix = 'burnin' if is_burnin else 'tree'
            treefn = '%s_%s_%s' % (prefix, idx, llh)
            self._archive.writestr(treefn, serialized_tree)
        self._close_archive()


class TreeReader(object):
    def __init__(self, archive_fn):
        self._archive = zipfile.ZipFile(archive_fn)
        infolist = self._archive.infolist()
        tree_info = [t for t in infolist if t.filename.startswith('tree_')]
        burnin_info = [t for t in infolist if t.filename.startswith('burnin_')]

        tree_info.sort(key=lambda tinfo: self._extract_metadata(tinfo)[0])
        burnin_info.sort(key=lambda tinfo: self._extract_burnin_idx(tinfo))

        self._trees = []
        self._burnin_trees = []

        for info in tree_info:
            idx, llh = self._extract_metadata(info)
            assert idx == len(self._trees)
            self._trees.append((idx, llh, info))
        for info in burnin_info:
            idx = self._extract_burnin_idx(info)
            if not len(burnin_info) + idx == len(self._burnin_trees):
                print('Burnin not finished, exiting', file=sys.stderr)
                sys.exit(1)
            self._burnin_trees.append((idx, info))
        assert len(tree_info) > 0

    def read_extra_file(self, filename):
        return self._archive.read(filename)

    def num_trees(self):
        return len(self._trees)

    def close(self):
        self._archive.close()

    def _extract_metadata(self, zinfo):
        tokens = zinfo.filename.split('_')
        idx = int(tokens[1])
        llh = float(tokens[2])
        return (idx, llh)

    def _extract_burnin_idx(self, zinfo):
        idx = int(zinfo.filename.split('_')[1])
        return idx

    def _parse_tree(self, zinfo, remove_empty_vertices=False):
        pickled = self._archive.read(zinfo)
        tree = pickle.loads(pickled)
        if remove_empty_vertices:
            remove_empty_nodes(tree.root)
        return tree

    def load_tree(self, idx, remove_empty_vertices=False):
        tidx, llh, zinfo = self._trees[idx]
        assert tidx == idx
        return self._parse_tree(zinfo, remove_empty_vertices)

    def load_trees(self, num_trees=None, remove_empty_vertices=False):
        for idx, llh, tree in self.load_trees_and_metadata(num_trees, remove_empty_vertices):
            yield tree

    def load_trees_and_burnin(self, remove_empty_vertices=False):
        for tidx, zinfo in self._burnin_trees:
            tree = self._parse_tree(zinfo, remove_empty_vertices)
            yield (tidx, tree)
        for tidx, llh, zinfo in self._trees:
            tree = self._parse_tree(zinfo, remove_empty_vertices)
            yield (tidx, tree)

    def load_trees_and_metadata(self, num_trees=None, remove_empty_vertices=False):
        trees = sorted(self._trees, key=lambda t: t[1], reverse=True)
        if num_trees is not None:
            num_trees = min(num_trees, len(trees))
            trees = trees[:num_trees]
        for tidx, llh, zinfo in trees:
            tree = self._parse_tree(zinfo, remove_empty_vertices)
            yield (tidx, llh, tree)


def logmsg(msg, fd=sys.stdout):
    print('[%s] %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), msg), file=fd)
