#!/usr/bin/env python3
"""
Node base class for TSSB tree nodes.

Optimizations vs original:
- get_ancestors() result is cached on the node; cache is invalidated
  recursively when parent changes. This eliminates repeated O(depth)
  list allocations throughout the hot MCMC paths.
- Python 3 compatible.
"""

import sys
from numpy.random import *
from numpy import *


class Node(object):

    def __init__(self, parent=None, tssb=None):
        self.data = set([])
        self._children = []
        self.tssb = tssb
        self._ancestors_cache = None  # cached result of get_ancestors()

        if parent is not None:
            parent.add_child(self)
            self._parent = parent
        else:
            self._parent = None

    def kill(self):
        if self._parent is not None:
            self._parent._children.remove(self)
        self._parent = None
        self._children = None
        self._ancestors_cache = None

    def spawn(self):
        return self.__class__(parent=self, tssb=self.tssb)

    def has_data(self):
        if len(self.data):
            return True
        else:
            for child in self._children:
                if child.has_data():
                    return True
        return False

    def num_data(self):
        return sum([c.num_data() for c in self._children], len(self.data))

    def num_local_data(self):
        return len(self.data)

    def add_datum(self, id):
        self.data.add(id)

    def remove_datum(self, id):
        self.data.remove(id)

    def resample_params(self):
        pass

    def add_child(self, child):
        self._children.append(child)

    def remove_child(self, child):
        self._children.remove(child)

    def children(self):
        return self._children

    def get_data(self):
        ids = list(self.data)
        return [self.tssb.data[id] for id in ids]

    def logprob(self, x):
        raise NotImplementedError

    def data_log_likelihood(self):
        return self.complete_logprob()

    def sample(self, num_data=1):
        return rand(num_data, 2)

    def parent(self):
        return self._parent

    def _set_parent(self, new_parent):
        """Change parent and invalidate ancestor caches recursively."""
        self._parent = new_parent
        self._invalidate_ancestors_cache()

    def _invalidate_ancestors_cache(self):
        """Recursively invalidate cached ancestors for self and all descendants."""
        self._ancestors_cache = None
        if self._children:
            for child in self._children:
                child._invalidate_ancestors_cache()

    def get_ancestors(self):
        """
        Return list of nodes from root down to (and including) self.

        Result is cached and reused until the tree structure changes.
        Cache is invalidated via _invalidate_ancestors_cache() when a
        node's parent is reassigned.
        """
        if self._ancestors_cache is not None:
            return self._ancestors_cache
        if self._parent is None:
            self._ancestors_cache = [self]
        else:
            self._ancestors_cache = self._parent.get_ancestors() + [self]
        return self._ancestors_cache

    def global_param(self, key):
        if self.parent() is None:
            return self.__dict__[key]
        else:
            return self.parent().global_param(key)
