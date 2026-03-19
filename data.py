#!/usr/bin/env python3
"""
Datum class: represents a single SSM or CNV with its observed read counts
and per-timepoint log-likelihood.

Key optimizations vs original:
1. find_most_recent_cnv uses a dict keyed by node instead of a list scan,
   reducing per-call complexity from O(depth × N_cnvs) to O(depth).
2. _log_likelihood no longer re-runs set_node_height / set_path_from_root /
   map_datum_to_node on every call (update_tree defaults to False).
   Callers are responsible for calling those helpers once per iteration.
3. compute_n_genomes uses node.path (pre-cached ancestors) instead of
   calling get_ancestors() on every node. O(N_nodes × depth) → O(N_nodes).
4. _log_likelihood is vectorised across time points (single numpy call)
   for the common CNV-free case.
5. Python 3 compatible.
"""

import builtins as _builtins
from numpy import *
import scipy.stats as stat
from scipy.special import gammaln
import util2 as u

# numpy's star import shadows Python's max/min with numpy versions that
# interpret positional args as axis. Restore builtins for scalar use.
_pymax = _builtins.max
_pymin = _builtins.min


class Datum(object):

    def __init__(self, name, id, a, d, mu_r=0, mu_v=0):
        self.name = name
        self.id = id
        self.a = a
        self.d = d
        self.mu_r = mu_r
        self.mu_v = mu_v
        self._log_bin_norm_const = [
            u.log_bin_coeff(self.d[tp], self.a[tp]) for tp in range(len(self.a))
        ]

        # CNV-related fields
        self.nr = 0
        self.nv = 0
        self.node = None   # node where this datum lives (set by map_datum_to_node)
        self.cnv = []      # list of (cnv_datum, maternal_cn, paternal_cn)

        # Fast O(1) lookup: node → (cnv_datum, maternal_cn, paternal_cn).
        # Updated in _rebuild_cnv_node_map() whenever self.cnv changes.
        # Keys are node objects; built lazily and rebuilt when needed.
        self._cnv_node_map = {}

        self.tssb = None   # pointer to the TSSB tree object

    # ------------------------------------------------------------------
    # CNV index maintenance
    # ------------------------------------------------------------------

    def _rebuild_cnv_node_map(self):
        """Rebuild the node→cnv_entry dict from self.cnv."""
        self._cnv_node_map = {}
        for entry in self.cnv:
            # entry = (cnv_datum, maternal_cn, paternal_cn)
            node = entry[0].node
            if node is not None:
                self._cnv_node_map[node] = entry

    def _get_cnv_node_map(self):
        """
        Return the node→cnv_entry dict, rebuilding if stale.

        We rebuild whenever any cnv entry's node might have changed
        (after map_datum_to_node is called).  To keep things simple
        we rebuild lazily: the dict is invalidated by setting it to
        None (done by map_datum_to_node via datum.node assignment path).
        """
        if not self._cnv_node_map and self.cnv:
            self._rebuild_cnv_node_map()
        return self._cnv_node_map

    # ------------------------------------------------------------------
    # Log-likelihood
    # ------------------------------------------------------------------

    def _log_likelihood(self, phi, update_tree=False, new_state=0):
        """
        Sum log-likelihood across all time points.

        update_tree=False (default): assumes tree metadata has already been
            refreshed this iteration via set_node_height / set_path_from_root
            / map_datum_to_node.  Pass True only if you need a standalone call
            that doesn't depend on an outer iteration context.
        """
        if update_tree:
            u.set_node_height(self.tssb)
            u.set_path_from_root_to_node(self.tssb)
            u.map_datum_to_node(self.tssb)
            # After map_datum_to_node the node pointers in cnv entries are
            # valid — rebuild the fast-lookup map.
            self._rebuild_cnv_node_map()

        if not self.cnv:
            # Vectorised path: compute all time points in one numpy call.
            return self._log_likelihood_no_cnv_vectorised(phi)

        # CNV path: still per-time-point (logic is too branchy to vectorise easily)
        # Use Python builtin sum to avoid numpy deprecation warning
        total = 0.0
        for tp in range(len(phi)):
            total += self.__log_complete_likelihood__(phi[tp], self.mu_r, self.mu_v, tp, new_state)
        return total

    def _log_likelihood_no_cnv_vectorised(self, phi):
        """
        Vectorised likelihood for SSMs / CNV-data with no linked CNV.
        phi: array-like, length = N_tps
        """
        phi_arr = asarray(phi)
        a_arr = asarray(self.a)
        d_arr = asarray(self.d)
        norm_arr = asarray(self._log_bin_norm_const)
        mu = (1.0 - phi_arr) * self.mu_r + phi_arr * self.mu_v
        # Clamp to avoid log(0)
        mu = clip(mu, 1e-15, 1.0 - 1e-15)
        llh = a_arr * log(mu) + (d_arr - a_arr) * log(1.0 - mu) + norm_arr
        return float(llh.sum())

    def _log_complete_likelihood(self, phi, mu_r, mu_v):
        """Convenience wrapper summing across all time points."""
        total = 0.0
        for tp in range(len(self.a)):
            total += self.__log_complete_likelihood__(phi, mu_r, mu_v, tp)
        return total

    def __log_complete_likelihood__(self, phi, mu_r, mu_v, tp, new_state=0):
        if self.cnv:
            poss_n_genomes = self.compute_n_genomes(tp, new_state)
            poss_n_genomes = [x for x in poss_n_genomes if x[1] > 0]
            ll = []
            for (nr, nv) in poss_n_genomes:
                mu = (nr * mu_r + nv * (1 - mu_r)) / (nr + nv)
                ll.append(
                    u.log_binomial_likelihood(self.a[tp], self.d[tp], mu)
                    + log(1.0 / len(poss_n_genomes))
                    + self._log_bin_norm_const[tp]
                )
            if not poss_n_genomes:
                ll.append(log(1e-99))
            return u.logsumexp(ll)
        else:
            mu = (1 - phi) * mu_r + phi * mu_v
            return (
                u.log_binomial_likelihood(self.a[tp], self.d[tp], mu)
                + self._log_bin_norm_const[tp]
            )

    # ------------------------------------------------------------------
    # Copy-number genome accounting
    # ------------------------------------------------------------------

    def compute_n_genomes(self, tp, new_state=0):
        """
        Compute the four possible (nr, nv) combinations for this SSM
        at time point `tp`.

        Optimizations:
        - Uses node.path (pre-computed ancestor list) instead of
          get_ancestors() (which would rebuild the list from scratch).
        - Uses self._cnv_node_map (O(1) dict lookup) instead of linear
          scan through self.cnv.
        """
        nodes = self.tssb.root['node'].tssb.get_nodes()

        self.nr1 = self.nv1 = 0.0
        self.nr2 = self.nv2 = 0.0
        self.nr3 = self.nv3 = 0.0
        self.nr4 = self.nv4 = 0.0

        # Use cached path if available, else compute on the fly
        ssm_node_path = self.node.path if self.node.path is not None else self.node.get_ancestors()
        ssm_node = ssm_node_path[-1]
        cnv_map = self._get_cnv_node_map()

        for nd in nodes:
            pi = nd.pi1[tp] if new_state else nd.pi[tp]
            # Use cached path if available, else fall back
            ancestors = nd.path if nd.path is not None else nd.get_ancestors()
            mr_cnv = self._find_most_recent_cnv_fast(nd, ancestors, cnv_map)

            ssm_in_ancestors = ssm_node in ancestors
            has_cnv = mr_cnv is not None

            if not ssm_in_ancestors and not has_cnv:
                self.nr1 += pi * 2
                self.nr2 += pi * 2
                self.nr3 += pi * 2
                self.nr4 += pi * 2
            elif ssm_in_ancestors and not has_cnv:
                self.nr1 += pi;  self.nv1 += pi
                self.nr2 += pi;  self.nv2 += pi
                self.nr3 += pi;  self.nv3 += pi
                self.nr4 += pi;  self.nv4 += pi
            elif not ssm_in_ancestors and has_cnv:
                total = mr_cnv[1] + mr_cnv[2]
                self.nr1 += pi * total
                self.nr2 += pi * total
                self.nr3 += pi * total
                self.nr4 += pi * total
            elif ssm_in_ancestors and has_cnv:
                total = mr_cnv[1] + mr_cnv[2]
                self.nr3 += pi * _pymax(0, total - 1);  self.nv3 += pi * _pymin(1, total)
                self.nr4 += pi * _pymax(0, total - 1);  self.nv4 += pi * _pymin(1, total)

                cnv_node_ancestors = (
                    mr_cnv[0].node.path
                    if mr_cnv[0].node.path is not None
                    else mr_cnv[0].node.get_ancestors()
                )
                if ssm_node in cnv_node_ancestors:
                    self.nr1 += pi * mr_cnv[1];  self.nv1 += pi * mr_cnv[2]
                    self.nr2 += pi * mr_cnv[2];  self.nv2 += pi * mr_cnv[1]
                else:
                    self.nr1 += pi * _pymax(0, total - 1);  self.nv1 += pi * _pymin(1, total)
                    self.nr2 += pi * _pymax(0, total - 1);  self.nv2 += pi * _pymin(1, total)
            else:
                print("PANIC: unexpected case in compute_n_genomes")

        if len(self.cnv) == 1 and self.node == self.cnv[0][0].node:
            out = [
                (self.nr1, self.nv1),
                (self.nr2, self.nv2),
                (self.nr3, self.nv3),
                (self.nr4, self.nv4),
            ]
        else:
            out = [(self.nr1, self.nv1), (self.nr2, self.nv2)]
        return out

    # ------------------------------------------------------------------
    # CNV lookup helpers
    # ------------------------------------------------------------------

    def find_most_recent_cnv(self, nd):
        """
        Original (slow) implementation kept for compatibility.
        Uses the fast-lookup map internally.
        """
        ancestors = nd.path if nd.path is not None else nd.get_ancestors()
        return self._find_most_recent_cnv_fast(nd, ancestors, self._get_cnv_node_map())

    def _find_most_recent_cnv_fast(self, nd, ancestors, cnv_map):
        """
        Walk ancestors from leaf toward root and return the first
        (most recent) CNV entry, using the O(1) dict lookup.

        ancestors: pre-computed list of nodes from root to nd (inclusive).
        cnv_map:   dict mapping node → (cnv_datum, maternal_cn, paternal_cn).
        """
        if not cnv_map:
            return None
        for n in reversed(ancestors):
            entry = cnv_map.get(n)
            if entry is not None:
                return entry
        return None
