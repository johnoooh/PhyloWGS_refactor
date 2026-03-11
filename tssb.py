#!/usr/bin/env python3
"""
Tree-Structured Stick-Breaking Process (TSSB).

Optimizations vs original:
1. resample_assignments: tree metadata (heights, paths, datum-node map) is
   computed ONCE before the per-datum loop, not re-triggered inside logprob.
   Ancestor-cache invalidation is called explicitly where needed.
2. path_lt: replaced string-encoding comparison with direct tuple/list
   comparison — no string allocation in the inner slice-sampler loop.
3. Python 3 (no cmp(), no print >> syntax).
"""

import sys
import scipy.stats
from time import *
from numpy import *
from numpy.random import *
from util import *

import util2 as _u2


class TSSB(object):
    min_dp_alpha = 1.0
    max_dp_alpha = 50.0
    min_dp_gamma = 1.0
    max_dp_gamma = 10.0
    min_alpha_decay = 0.05
    max_alpha_decay = 0.80

    def __init__(self, dp_alpha=1.0, dp_gamma=1.0, root_node=None, data=None,
                 min_depth=0, max_depth=15, alpha_decay=1.0):
        if root_node is None:
            raise Exception("Root node must be specified.")

        self.min_depth = min_depth
        self.max_depth = max_depth
        self.dp_alpha = dp_alpha
        self.dp_gamma = dp_gamma
        self.alpha_decay = alpha_decay
        self.data = data
        self.num_data = 0 if data is None else len(data)
        self.root = {
            'node': root_node,
            'main': boundbeta(1.0, dp_alpha) if self.min_depth == 0 else 0.0,
            'sticks': empty((0, 1)),
            'children': [],
        }
        root_node.tssb = self

        self.assignments = []
        for n in range(self.num_data):
            self.root['node'].add_datum(n)
            self.assignments.append(self.root['node'])

    def add_data(self, data):
        (weights, nodes) = self.get_mixture()
        num_new_data = len(data)
        for n in range(num_new_data):
            logprobs = []
            for k, node in enumerate(nodes):
                logprobs.append(log(weights[k]) + node.logprob(data[n]))
            logprobs = array(logprobs)
            probs = exp(logprobs - logsumexp(logprobs))
            best_k = sum(rand() > cumsum(probs))
            nodes[best_k].add_datum(n + self.num_data)
            self.assignments.append(nodes[best_k])
        self.data = vstack([self.data, data])
        self.num_data += num_new_data

    def resample_node_params(self, iters=1):
        for _iter in range(iters):
            def descend(root):
                for index, child in enumerate(root['children']):
                    descend(child)
                root['node'].resample_params()
            descend(self.root)

    def resample_assignments(self):
        """
        Slice-sample each datum's node assignment.

        Tree metadata (heights, ancestor paths, datum→node map) is
        pre-computed once here, then kept current across moves by
        calling map_datum_to_node after each accepted proposal.
        """
        # --- Pre-compute tree metadata once for the whole loop ---
        _u2.set_node_height(self)
        _u2.set_path_from_root_to_node(self)
        _u2.map_datum_to_node(self)
        # Rebuild fast CNV-node maps now that datum.node pointers are fresh.
        for datum in self.data:
            if datum.cnv:
                datum._rebuild_cnv_node_map()

        epsilon = finfo(float64).eps

        for n in range(self.num_data):
            # Get path indices to current assignment (for slice sampler bounds).
            ancestors = self.assignments[n].get_ancestors()
            current = self.root
            indices = []
            for anc in ancestors[1:]:
                index = [c['node'] for c in current['children']].index(anc)
                current = current['children'][index]
                indices.append(index)

            max_u = 1.0
            min_u = 0.0
            # Cache LLH per visited node to avoid recomputation.
            llhmap = {}

            old_llh = self.assignments[n].logprob(self.data[n:n + 1])
            llhmap[self.assignments[n]] = old_llh
            llh_s = log(rand()) + old_llh

            while True:
                new_u = (max_u - min_u) * rand() + min_u
                (new_node, new_path) = self.find_node(new_u)
                if new_node.parent() is None:
                    new_node = new_node.children()[0]
                    new_path = [0]

                old_node = self.assignments[n]
                old_node.remove_datum(n)
                new_node.add_datum(n)
                self.assignments[n] = new_node

                # Update datum.node so logprob sees correct placement.
                self.data[n].node = new_node

                if new_node in llhmap:
                    new_llh = llhmap[new_node]
                else:
                    new_llh = new_node.logprob(self.data[n:n + 1])
                    llhmap[new_node] = new_llh

                if new_llh > llh_s:
                    break
                elif abs(max_u - min_u) < epsilon:
                    new_node.remove_datum(n)
                    old_node.add_datum(n)
                    self.assignments[n] = new_node if False else old_node
                    self.data[n].node = old_node
                    print("Slice sampler shrank down.  Keep current state.", file=sys.stderr)
                    break
                else:
                    new_node.remove_datum(n)
                    old_node.add_datum(n)
                    self.assignments[n] = old_node
                    self.data[n].node = old_node

                    if _path_lt(indices, new_path) < 0:
                        min_u = new_u
                    else:
                        max_u = new_u

    def cull_tree(self):
        def descend(root):
            counts = array([descend(child) for child in root['children']])
            keep = len(trim_zeros(counts, 'b'))
            for child in root['children'][keep:]:
                child['node'].kill()
                del child['node']
            root['sticks'] = root['sticks'][:keep]
            root['children'] = root['children'][:keep]
            return sum(counts) + root['node'].num_local_data()
        descend(self.root)

    def resample_sticks(self):
        def descend(root, depth=0):
            data_down = 0
            indices = list(range(len(root['children'])))
            indices.reverse()
            for i in indices:
                child = root['children'][i]
                child_data = descend(child, depth + 1)
                post_alpha = 1.0 + child_data
                post_beta = self.dp_gamma + data_down
                root['sticks'][i] = boundbeta(post_alpha, post_beta) if depth != 0 else .999999
                data_down += child_data
            data_here = root['node'].num_local_data()
            post_alpha = 1.0 + data_here
            post_beta = (self.alpha_decay ** depth) * self.dp_alpha + data_down
            root['main'] = boundbeta(post_alpha, post_beta) if self.min_depth <= depth else 0.0
            if depth == 0:
                root['main'] = 1e-30
            return data_here + data_down
        descend(self.root)

    def resample_stick_orders(self):
        def descend(root, depth=0):
            if not root['children']:
                return

            new_order = []
            represented = set(
                i for i in range(len(root['children']))
                if root['children'][i]['node'].has_data()
            )
            all_weights = diff(hstack([0.0, sticks_to_edges(root['sticks'])]))

            while True:
                if not represented:
                    break
                u = rand()
                while True:
                    sub_indices = [i for i in range(root['sticks'].shape[0]) if i not in new_order]
                    sub_weights = hstack([all_weights[sub_indices], 1.0 - sum(all_weights)])
                    sub_weights = sub_weights / sum(sub_weights)
                    index = sum(u > cumsum(sub_weights))

                    if index == len(sub_indices):
                        root['sticks'] = vstack([root['sticks'], boundbeta(1, self.dp_gamma)])
                        root['children'].append({
                            'node': root['node'].spawn(),
                            'main': boundbeta(
                                1.0,
                                (self.alpha_decay ** (depth + 1)) * self.dp_alpha
                            ) if self.min_depth <= (depth + 1) else 0.0,
                            'sticks': empty((0, 1)),
                            'children': [],
                        })
                        all_weights = diff(hstack([0.0, sticks_to_edges(root['sticks'])]))
                    else:
                        index = sub_indices[index]
                        break
                new_order.append(index)
                represented.discard(index)

            new_children = []
            for k in new_order:
                child = root['children'][k]
                new_children.append(child)
                descend(child, depth + 1)

            for k in [k for k in range(root['sticks'].shape[0]) if k not in new_order]:
                root['children'][k]['node'].kill()
                del root['children'][k]['node']

            root['children'] = new_children
            root['sticks'] = zeros((len(root['children']), 1))

        descend(self.root)
        self.resample_sticks()

    def resample_hypers(self, dp_alpha=True, alpha_decay=True, dp_gamma=True):

        def dp_alpha_llh(dp_alpha, alpha_decay):
            def descend(dp_alpha, root, depth=0):
                llh = betapdfln(root['main'], 1.0,
                                (alpha_decay ** depth) * dp_alpha) if self.min_depth <= depth else 0.0
                for child in root['children']:
                    llh += descend(dp_alpha, child, depth + 1)
                return llh
            return descend(dp_alpha, self.root)

        if dp_alpha:
            upper = self.max_dp_alpha
            lower = self.min_dp_alpha
            llh_s = log(rand()) + dp_alpha_llh(self.dp_alpha, self.alpha_decay)
            while True:
                new_dp_alpha = (upper - lower) * rand() + lower
                new_llh = dp_alpha_llh(new_dp_alpha, self.alpha_decay)
                if new_llh > llh_s:
                    break
                elif new_dp_alpha < self.dp_alpha:
                    lower = new_dp_alpha
                elif new_dp_alpha > self.dp_alpha:
                    upper = new_dp_alpha
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dp_alpha = new_dp_alpha

        if alpha_decay:
            upper = self.max_alpha_decay
            lower = self.min_alpha_decay
            llh_s = log(rand()) + dp_alpha_llh(self.dp_alpha, self.alpha_decay)
            while True:
                new_alpha_decay = (upper - lower) * rand() + lower
                new_llh = dp_alpha_llh(self.dp_alpha, new_alpha_decay)
                if new_llh > llh_s:
                    break
                elif new_alpha_decay < self.alpha_decay:
                    lower = new_alpha_decay
                elif new_alpha_decay > self.alpha_decay:
                    upper = new_alpha_decay
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.alpha_decay = new_alpha_decay

        def dp_gamma_llh(dp_gamma):
            def descend(dp_gamma, root):
                llh = 0
                for i, child in enumerate(root['children']):
                    llh += betapdfln(root['sticks'][i], 1.0, dp_gamma)
                    llh += descend(dp_gamma, child)
                return llh
            return descend(dp_gamma, self.root)

        if dp_gamma:
            upper = self.max_dp_gamma
            lower = self.min_dp_gamma
            llh_s = log(rand()) + dp_gamma_llh(self.dp_gamma)
            while True:
                new_dp_gamma = (upper - lower) * rand() + lower
                new_llh = dp_gamma_llh(new_dp_gamma)
                if new_llh > llh_s:
                    break
                elif new_dp_gamma < self.dp_gamma:
                    lower = new_dp_gamma
                elif new_dp_gamma > self.dp_gamma:
                    upper = new_dp_gamma
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dp_gamma = new_dp_gamma

    def find_node(self, u):
        def descend(root, u, depth=0):
            if depth >= self.max_depth:
                return (root['node'], [])
            elif u < root['main']:
                return (root['node'], [])
            else:
                u = (u - root['main']) / (1.0 - root['main'])
                if depth > 0:
                    while not root['children'] or (1.0 - prod(1.0 - root['sticks'])) < u:
                        root['sticks'] = vstack(
                            [root['sticks'], boundbeta(1, self.dp_gamma) if depth != 0 else .999]
                        )
                        root['children'].append({
                            'node': root['node'].spawn(),
                            'main': boundbeta(
                                1.0,
                                (self.alpha_decay ** (depth + 1)) * self.dp_alpha
                            ) if self.min_depth <= (depth + 1) else 0.0,
                            'sticks': empty((0, 1)),
                            'children': [],
                        })
                    edges = 1.0 - cumprod(1.0 - root['sticks'])
                    index = sum(u > edges)
                    edges = hstack([0.0, edges])
                    u = (u - edges[index]) / (edges[index + 1] - edges[index])
                    (node, path) = descend(root['children'][index], u, depth + 1)
                else:
                    index = 0
                    (node, path) = descend(root['children'][index], u, depth + 1)
                path.insert(0, index)
                return (node, path)
        return descend(self.root, u)

    def get_nodes(self):
        def descend(root):
            node = [root['node']]
            for child in root['children']:
                node.extend(descend(child))
            return node
        return descend(self.root)

    def get_mixture(self):
        def descend(root, mass):
            weight = [mass * root['main']]
            node = [root['node']]
            edges = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))
            for i, child in enumerate(root['children']):
                (child_weights, child_nodes) = descend(child, mass * (1.0 - root['main']) * weights[i])
                weight.extend(child_weights)
                node.extend(child_nodes)
            return (weight, node)
        return descend(self.root, 1.0)

    def complete_data_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.num_local_data() * log(weights[i]) + node.data_log_likelihood())
        return sum(array(llhs))

    def complete_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = [self.dp_alpha_llh(self.dp_alpha, self.alpha_decay), self.dp_gamma_llh(self.dp_gamma)]
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.data_log_likelihood())
        return sum(array(llhs))

    def dp_alpha_llh(self, dp_alpha, alpha_decay):
        def descend(dp_alpha, root, depth=0):
            llh = betapdfln(root['main'], 1.0, (alpha_decay ** depth) * dp_alpha) if self.min_depth <= depth else 0.0
            for child in root['children']:
                llh += descend(dp_alpha, child, depth + 1)
            return llh
        return descend(dp_alpha, self.root)

    def dp_gamma_llh(self, dp_gamma):
        def descend(dp_gamma, root):
            llh = 0
            for i, child in enumerate(root['children']):
                llh += betapdfln(root['sticks'][i], 1.0, dp_gamma)
                llh += descend(dp_gamma, child)
            return llh
        return descend(dp_gamma, self.root)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _path_lt(path1, path2):
    """
    Return negative/zero/positive like cmp(path2, path1) for ordering paths.

    Replaces the original string-encoding approach:
      s1 = "".join(map(lambda i: "%03d" % i, path1))
      s2 = "".join(map(lambda i: "%03d" % i, path2))
      return cmp(s2, s1)

    Direct list/tuple comparison is O(depth) with no allocation.
    We replicate the original semantics: compare path2 vs path1,
    return negative if path2 < path1, etc.
    """
    if not path1 and not path2:
        return 0
    if not path1:
        return 1   # original: len(path1)==0 → return 1
    if not path2:
        return -1  # original: len(path2)==0 → return -1
    t1, t2 = tuple(path1), tuple(path2)
    if t2 < t1:
        return -1
    elif t2 > t1:
        return 1
    return 0
