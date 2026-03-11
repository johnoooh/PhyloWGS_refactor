#!/usr/bin/env python3
"""
PhyloWGS Benchmark: Compare optimized vs simulated-original performance.

Usage:
    python3 benchmark.py [--burnin N] [--samples N] [--mh-itr N]

Default: 50 burnin + 200 samples = 250 iterations, mh_itr=100
"""

import os
import sys
import time
import tempfile
import shutil
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

import numpy as np
from numpy import *
from numpy.random import seed

from tssb import TSSB
from alleles import alleles
from util import boundbeta
from util2 import load_data, set_node_height, set_path_from_root_to_node, map_datum_to_node
from params import metropolis


def setup_tssb(ssm_file, cnv_file):
    """Initialize a TSSB tree for benchmarking."""
    codes, n_ssms, n_cnvs, cnv_map = load_data(ssm_file, cnv_file)
    NTPS = len(codes[0].a)
    
    root = alleles(conc=0.1, ntps=NTPS)
    tssb = TSSB(
        dp_alpha=25.0, dp_gamma=1.0, alpha_decay=0.25,
        root_node=root, data=codes,
    )
    
    tssb.root['sticks'] = vstack([tssb.root['sticks'], .999999])
    tssb.root['children'].append({
        'node': tssb.root['node'].spawn(),
        'main': boundbeta(1.0, tssb.alpha_decay * tssb.dp_alpha),
        'sticks': empty((0, 1)),
        'children': [],
    })
    new_node = tssb.root['children'][0]['node']
    for n in range(tssb.num_data):
        tssb.assignments[n].remove_datum(n)
        new_node.add_datum(n)
        tssb.assignments[n] = new_node
        
    for datum in codes:
        datum.tssb = tssb
        
    return tssb, codes, n_ssms, n_cnvs, NTPS


def run_iteration(tssb, codes, n_ssms, n_cnvs, NTPS, ssm_file, cnv_file, tmp_dir, mh_itr, slow=False):
    """Run a single MCMC iteration. If slow=True, force update_tree=True."""
    original_logprob = None
    if slow:
        original_logprob = alleles.logprob
        def slow_logprob(self, x):
            return x[0]._log_likelihood(self.params, update_tree=True)
        alleles.logprob = slow_logprob
    
    try:
        tssb.resample_assignments()
        tssb.cull_tree()
        
        wts, nodes = tssb.get_mixture()
        for i, node in enumerate(nodes):
            node.id = i
            
        set_node_height(tssb)
        set_path_from_root_to_node(tssb)
        map_datum_to_node(tssb)
        
        for datum in codes:
            if datum.cnv:
                datum._rebuild_cnv_node_map()
                
        metropolis(tssb, mh_itr, 100.0, 0, n_ssms, n_cnvs, ssm_file, cnv_file, 42, NTPS, tmp_dir)
        
        tssb.resample_sticks()
        tssb.resample_stick_orders()
        tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
    finally:
        if original_logprob:
            alleles.logprob = original_logprob


def main():
    parser = argparse.ArgumentParser(description='PhyloWGS benchmark')
    parser.add_argument('--burnin', type=int, default=50)
    parser.add_argument('--samples', type=int, default=200)
    parser.add_argument('--mh-itr', type=int, default=100)
    args = parser.parse_args()
    
    total_iters = args.burnin + args.samples
    ssm_file = os.path.join(script_dir, 'ssm_data.txt')
    cnv_file = os.path.join(script_dir, 'cnv_data.txt')
    tmp_dir = tempfile.mkdtemp(prefix='phylowgs_bench.')
    
    print(f"PhyloWGS Benchmark: {total_iters} iterations, mh_itr={args.mh_itr}")
    
    try:
        # Fast path
        print("Running FAST path...")
        seed(12345)
        tssb, codes, n_ssms, n_cnvs, NTPS = setup_tssb(ssm_file, cnv_file)
        start = time.time()
        for i in range(total_iters):
            run_iteration(tssb, codes, n_ssms, n_cnvs, NTPS, ssm_file, cnv_file, tmp_dir, args.mh_itr, slow=False)
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{total_iters}...")
        fast_time = time.time() - start
        
        # Slow path
        print("Running SLOW path...")
        seed(12345)
        tssb, codes, n_ssms, n_cnvs, NTPS = setup_tssb(ssm_file, cnv_file)
        start = time.time()
        for i in range(total_iters):
            run_iteration(tssb, codes, n_ssms, n_cnvs, NTPS, ssm_file, cnv_file, tmp_dir, args.mh_itr, slow=True)
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{total_iters}...")
        slow_time = time.time() - start
        
        speedup = slow_time / fast_time
        print(f"\nResults:")
        print(f"  Fast: {fast_time:.2f}s ({fast_time/total_iters*1000:.2f}ms/iter)")
        print(f"  Slow: {slow_time:.2f}s ({slow_time/total_iters*1000:.2f}ms/iter)")
        print(f"  Speedup: {speedup:.2f}x")
        
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == '__main__':
    main()
