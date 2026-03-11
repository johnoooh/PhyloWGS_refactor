#!/usr/bin/env python3
"""
Math utility functions for PhyloWGS.

Optimizations vs original:
- Python 3 (no cPickle, no Python-2-only syntax).
- logsumexp delegates to scipy.special for accuracy.
- Removed dead code (hotshot profiler, Python-2 bin() shim).
"""

import os
import sys
import pickle
import numpy
import scipy.special
import scipy.stats
from functools import reduce


def bucket(edges, value):
    return numpy.sum(value > edges)


def sticks_to_edges(sticks):
    return 1.0 - numpy.cumprod(1.0 - sticks)


def normpdfln(x, m, prec):
    return numpy.sum(
        -0.5 * numpy.log(2 * numpy.pi) + 0.5 * numpy.log(prec) - 0.5 * prec * (x - m) ** 2,
        axis=1,
    )


def gammaln(x):
    return scipy.special.gammaln(x)


def gammapdfln(x, a, b):
    return -gammaln(a) + a * numpy.log(b) + (a - 1.0) * numpy.log(x) - b * x


def exp_gammapdfln(y, a, b):
    return a * numpy.log(b) - gammaln(a) + a * y - b * numpy.exp(y)


def betapdfln(x, a, b):
    return (
        gammaln(a + b)
        - gammaln(a)
        - gammaln(b)
        + (a - 1.0) * numpy.log(x)
        + (b - 1.0) * numpy.log(1.0 - x)
    )


def boundbeta(a, b):
    return (1.0 - numpy.finfo(numpy.float64).eps) * (numpy.random.beta(a, b) - 0.5) + 0.5


def lnbetafunc(a):
    return numpy.sum(gammaln(a)) - gammaln(numpy.sum(a))


def dirichletpdfln(p, a):
    return -lnbetafunc(a) + numpy.sum((a - 1) * numpy.log(p))


def logsumexp(X, axis=None):
    """Numerically stable log-sum-exp using scipy."""
    return scipy.special.logsumexp(X, axis=axis)


def merge(l):
    return [item for sublist in l for item in sublist]


def slice_sample(
    init_x,
    logprob,
    sigma=1.0,
    step_out=True,
    max_steps_out=1000,
    compwise=False,
    verbose=False,
):
    def direction_slice(direction, init_x):
        def dir_logprob(z):
            return logprob(direction * z + init_x)

        r = numpy.random.rand()
        upper = init_x + (1 - r) * sigma
        lower = init_x - r * sigma
        llh_s = numpy.log(numpy.random.rand()) + dir_logprob(0.0)
        l_steps_out = 0
        u_steps_out = 0
        if step_out:
            while dir_logprob(lower) > llh_s and l_steps_out < max_steps_out:
                l_steps_out += 1
                lower -= sigma
            while dir_logprob(upper) > llh_s and u_steps_out < max_steps_out:
                u_steps_out += 1
                upper += sigma

        steps_in = 0
        while True:
            steps_in += 1
            new_z = (upper - lower) * numpy.random.rand()
            new_llh = dir_logprob(new_z)
            if numpy.isnan(new_llh):
                raise Exception("Slice sampler got a NaN")
            if new_llh > llh_s:
                break
            elif new_z < 0:
                lower = new_z
            elif new_z > 0:
                upper = new_z
            else:
                raise Exception("Slice sampler shrank to zero!")

        if verbose:
            print("Steps Out:", l_steps_out, u_steps_out, " Steps In:", steps_in)

        return new_z * direction + init_x

    dims = 1
    if compwise:
        ordering = list(range(dims))
        numpy.random.shuffle(ordering)
        cur_x = init_x + 0.0
        for d in ordering:
            direction = numpy.zeros((dims))
            direction[d] = 1.0
            cur_x = direction_slice(direction, cur_x)
        return cur_x
    else:
        direction = numpy.random.randn(dims)
        direction = direction / numpy.sqrt(numpy.sum(direction ** 2))
        return direction_slice(direction, init_x)
