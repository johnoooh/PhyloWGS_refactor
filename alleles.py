#!/usr/bin/env python3
"""
alleles.Node — the clone/node class holding per-timepoint cellular frequency
parameters (phi / pi).

Optimization vs original:
- logprob() calls _log_likelihood(update_tree=False) — the caller
  (resample_assignments in tssb.py) is responsible for ensuring tree metadata
  is fresh before the loop begins.  This removes O(N_nodes) tree traversals
  that the original triggered on *every datum* during assignment resampling.
- Python 3 compatible.
"""

from scipy.stats import beta, binom
import scipy.stats as stat
from scipy.special import comb
from util import *
from numpy import *
from node import *

import util2 as _u2


class alleles(Node):

    init_mean = 0.5
    min_conc = 0.01
    max_conc = 0.1

    def __init__(self, parent=None, tssb=None, conc=0.1, ntps=5):
        super(alleles, self).__init__(parent=parent, tssb=tssb)

        if tssb is not None:
            ntps = len(tssb.data[0].a)

        self.pi = zeros(ntps)
        self.params = zeros(ntps)
        self.params1 = zeros(ntps)
        self.pi1 = zeros(ntps)   # used in MH to store old state

        self.path = None   # set of nodes from root to this node
        self.ht = 0

        if parent is None:
            self._conc = conc
            self.pi = 1.0 * ones(ntps)
            self.params = 1.0 * ones(ntps)
        else:
            self.pi = rand(1) * parent.pi
            parent.pi = parent.pi - self.pi
            self.params = self.pi

    def conc(self):
        if self.parent() is None:
            return self._conc
        else:
            return self.parent().conc()

    def kill(self):
        if self._parent is not None:
            self._parent._children.remove(self)
        self._parent.pi = self._parent.pi + self.pi
        self._parent = None
        self._children = None
        self._ancestors_cache = None

    def logprob(self, x):
        """
        Log probability of data item x under this node's params.

        update_tree=False: tree metadata (heights, ancestor paths, datum→node
        map) must already be current — set_node_height / set_path_from_root /
        map_datum_to_node should have been called once before the enclosing
        loop, not on every datum.
        """
        return x[0]._log_likelihood(self.params, update_tree=False)

    def complete_logprob(self):
        return sum([self.logprob([data]) for data in self.get_data()])
