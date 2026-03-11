#!/usr/bin/env python3
"""
PhyloWGS MCMC driver.

Optimizations vs original:
- Python 3 throughout (pickle, print(), threading, argparse).
- do_mcmc calls set_node_height / set_path_from_root_to_node /
  map_datum_to_node once per iteration (before metropolis()), not per datum.
  The tssb.resample_assignments() method handles its own pre-computation.
- After map_datum_to_node, cnv_node_maps are rebuilt on all data so that
  the fast O(1) lookup dict stays current.
"""

import os
import sys
import pickle

from numpy import *
from numpy.random import *
from tssb import *
from alleles import *
from util import *

import numpy.random

from util2 import *
from params import *
from printo import *

import argparse
import signal
import tempfile
import threading
import traceback
import time
import json
from datetime import datetime
from collections import OrderedDict


def start_new_run(state_manager, backup_manager, safe_to_exit, run_succeeded,
                  config, ssm_file, cnv_file, params_file,
                  burnin_samples, num_samples, mh_itr, mh_std,
                  write_state_every, write_backups_every, rand_seed, tmp_dir):
    state = {}

    try:
        state['rand_seed'] = int(rand_seed)
    except TypeError:
        try:
            with open('random_seed.txt') as seedf:
                state['rand_seed'] = int(seedf.read().strip())
        except (TypeError, IOError):
            state['rand_seed'] = randint(2 ** 32)

    seed(state['rand_seed'])
    with open('random_seed.txt', 'w') as seedf:
        seedf.write('%s\n' % state['rand_seed'])

    state['working_dir'] = os.getcwd()
    state['ssm_file'] = ssm_file
    state['cnv_file'] = cnv_file
    state['tmp_dir'] = tmp_dir
    state['write_state_every'] = write_state_every
    state['write_backups_every'] = write_backups_every

    codes, n_ssms, n_cnvs, cnv_logical_physical_mapping = load_data(
        state['ssm_file'], state['cnv_file']
    )
    if len(codes) == 0:
        logmsg('No SSMs or CNVs provided. Exiting.', sys.stderr)
        return
    NTPS = len(codes[0].a)
    state['glist'] = [datum.name for datum in codes if len(datum.name) > 0]

    state['burnin'] = burnin_samples
    state['num_samples'] = num_samples
    state['dp_alpha'] = 25.0
    state['dp_gamma'] = 1.0
    state['alpha_decay'] = 0.25

    state['mh_burnin'] = 0
    state['mh_itr'] = mh_itr
    state['mh_std'] = mh_std

    state['cd_llh_traces'] = zeros((state['num_samples'], 1))
    state['burnin_cd_llh_traces'] = zeros((state['burnin'], 1))

    root = alleles(conc=0.1, ntps=NTPS)
    state['tssb'] = TSSB(
        dp_alpha=state['dp_alpha'],
        dp_gamma=state['dp_gamma'],
        alpha_decay=state['alpha_decay'],
        root_node=root,
        data=codes,
    )

    # Initialise tree with one child node containing all data.
    if True:
        depth = 0
        state['tssb'].root['sticks'] = vstack([
            state['tssb'].root['sticks'],
            boundbeta(1, state['tssb'].dp_gamma) if depth != 0 else .999999,
        ])
        state['tssb'].root['children'].append({
            'node': state['tssb'].root['node'].spawn(),
            'main': boundbeta(1.0, (state['tssb'].alpha_decay ** (depth + 1)) * state['tssb'].dp_alpha)
                    if state['tssb'].min_depth <= (depth + 1) else 0.0,
            'sticks': empty((0, 1)),
            'children': [],
        })
        new_node = state['tssb'].root['children'][0]['node']
        for n in range(state['tssb'].num_data):
            state['tssb'].assignments[n].remove_datum(n)
            new_node.add_datum(n)
            state['tssb'].assignments[n] = new_node

    for datum in codes:
        datum.tssb = state['tssb']

    tree_writer = TreeWriter()
    tree_writer.add_extra_file('cnv_logical_physical_mapping.json',
                               json.dumps(cnv_logical_physical_mapping))

    if params_file is not None:
        with open(params_file) as F:
            params = json.load(F)
    else:
        params = {}
    tree_writer.add_extra_file('params.json', json.dumps(params))

    state_manager.write_initial_state(state)
    logmsg("Starting MCMC run...")
    state['last_iteration'] = -state['burnin'] - 1

    with open('mcmc_samples.txt', 'w') as mcmcf:
        mcmcf.write('Iteration\tLLH\tTime\n')

    do_mcmc(state_manager, backup_manager, safe_to_exit, run_succeeded,
            config, state, tree_writer, codes, n_ssms, n_cnvs, NTPS, tmp_dir)


def resume_existing_run(state_manager, backup_manager, safe_to_exit,
                        run_succeeded, config):
    try:
        state = state_manager.load_state()
        os.chdir(state['working_dir'])
        tree_writer = TreeWriter(resume_run=True)
    except Exception:
        logmsg('Restoring state failed:', sys.stderr)
        traceback.print_exc()
        logmsg('Restoring from backup and trying again.', sys.stderr)
        backup_manager.restore_backup()
        state = state_manager.load_state()
        os.chdir(state['working_dir'])
        tree_writer = TreeWriter(resume_run=True)

    numpy.random.set_state(state['rand_state'])
    codes, n_ssms, n_cnvs, cnv_logical_physical_mapping = load_data(
        state['ssm_file'], state['cnv_file']
    )
    NTPS = len(codes[0].a)

    do_mcmc(state_manager, backup_manager, safe_to_exit, run_succeeded,
            config, state, tree_writer, codes, n_ssms, n_cnvs, NTPS,
            state['tmp_dir'])


def do_mcmc(state_manager, backup_manager, safe_to_exit, run_succeeded,
            config, state, tree_writer, codes, n_ssms, n_cnvs, NTPS,
            tmp_dir_parent):
    start_iter = state['last_iteration'] + 1
    unwritten_trees = []
    mcmc_sample_times = []
    last_mcmc_sample_time = time.time()

    config['tmp_dir'] = tempfile.mkdtemp(prefix='pwgsdataexchange.', dir=tmp_dir_parent)

    for iteration in range(start_iter, state['num_samples']):
        sys.stdout.flush()
        safe_to_exit.set()

        status = OrderedDict()
        status['iteration'] = iteration
        status['trees_sampled'] = state['burnin'] + iteration
        status['total_trees'] = state['burnin'] + state['num_samples']

        tssb = state['tssb']

        # resample_assignments pre-computes tree metadata internally.
        tssb.resample_assignments()
        tssb.cull_tree()

        wts, nodes = tssb.get_mixture()
        for i, node in enumerate(nodes):
            node.id = i

        # These three helpers prepare per-node height/path/datum info used
        # by metropolis() and the likelihood computations that follow.
        set_node_height(tssb)
        set_path_from_root_to_node(tssb)
        map_datum_to_node(tssb)

        # Keep cnv_node_maps fresh after datum→node assignment changed.
        for datum in codes:
            if datum.cnv:
                datum._rebuild_cnv_node_map()

        state['mh_acc'] = metropolis(
            tssb,
            state['mh_itr'],
            state['mh_std'],
            state['mh_burnin'],
            n_ssms,
            n_cnvs,
            state['ssm_file'],
            state['cnv_file'],
            state['rand_seed'],
            NTPS,
            config['tmp_dir'],
        )

        if float(state['mh_acc']) < 0.08 and state['mh_std'] < 10000:
            state['mh_std'] = state['mh_std'] * 2.0
            logmsg("Shrinking MH proposals. Now %f" % state['mh_std'])
        if float(state['mh_acc']) > 0.5 and float(state['mh_acc']) < 0.99:
            state['mh_std'] = state['mh_std'] / 2.0
            logmsg("Growing MH proposals. Now %f" % state['mh_std'])

        tssb.resample_sticks()
        tssb.resample_stick_orders()
        tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)

        last_llh = tssb.complete_data_log_likelihood()
        status['llh'] = last_llh

        if iteration >= 0:
            state['cd_llh_traces'][iteration] = last_llh
            weights, nodes = tssb.get_mixture()
            status['nodes'] = len(nodes)
            status['mh_acc'] = state['mh_acc']
            status['dp_alpha'] = tssb.dp_alpha
            status['dp_gamma'] = tssb.dp_gamma
            status['alpha_decay'] = tssb.alpha_decay
        else:
            state['burnin_cd_llh_traces'][iteration + state['burnin']] = last_llh

        logmsg(' '.join(['%s=%s' % (K, V) for K, V in status.items()]))

        serialized = pickle.dumps(tssb, protocol=pickle.HIGHEST_PROTOCOL)
        unwritten_trees.append((serialized, iteration, last_llh))
        state['tssb'] = tssb
        state['rand_state'] = numpy.random.get_state()
        state['last_iteration'] = iteration

        if len([C for C in state['tssb'].root['children'] if C['node'].has_data()]) > 1:
            logmsg('Polyclonal tree detected with %s clones.' % len(state['tssb'].root['children']))

        new_mcmc_sample_time = time.time()
        mcmc_sample_times.append(new_mcmc_sample_time - last_mcmc_sample_time)
        last_mcmc_sample_time = new_mcmc_sample_time

        safe_to_exit.clear()
        should_write_backup = (
            iteration % state['write_backups_every'] == 0 and iteration != start_iter
        )
        should_write_state = iteration % state['write_state_every'] == 0
        is_last_iteration = iteration == state['num_samples'] - 1

        if should_write_backup or should_write_state or is_last_iteration:
            with open('mcmc_samples.txt', 'a') as mcmcf:
                lines = '\n'.join(
                    '%s\t%s\t%s' % (itr, llh, itr_time)
                    for (_, itr, llh), itr_time in zip(unwritten_trees, mcmc_sample_times)
                )
                mcmcf.write(lines + '\n')
            tree_writer.write_trees(unwritten_trees)
            state_manager.write_state(state)
            unwritten_trees = []
            mcmc_sample_times = []
            if should_write_backup:
                backup_manager.save_backup()

    backup_manager.remove_backup()
    safe_to_exit.clear()

    freq = {g: [] for g in state['glist']}
    safe_to_exit.set()
    run_succeeded.set()


def create_argparser():
    parser = argparse.ArgumentParser(
        description='Run PhyloWGS to infer subclonal composition from SSMs and CNVs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-O', '--output-dir', dest='output_dir',
                        help='Path to output directory')
    return parser


def switch_working_dir():
    parser = create_argparser()
    args, other_args = parser.parse_known_args()
    working_dir = args.output_dir
    orig_working_dir = os.getcwd()
    if working_dir is not None:
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        os.chdir(working_dir)
    return orig_working_dir


def create_argparser_with_all_args():
    parser = create_argparser()
    parser.add_argument('-b', '--write-backups-every', dest='write_backups_every',
                        default=100, type=int)
    parser.add_argument('-S', '--write-state-every', dest='write_state_every',
                        default=10, type=int)
    parser.add_argument('-B', '--burnin-samples', dest='burnin_samples',
                        default=1000, type=int)
    parser.add_argument('-s', '--mcmc-samples', dest='mcmc_samples',
                        default=2500, type=int)
    parser.add_argument('-i', '--mh-iterations', dest='mh_iterations',
                        default=5000, type=int)
    parser.add_argument('-r', '--random-seed', dest='random_seed', type=int)
    parser.add_argument('-t', '--tmp-dir', dest='tmp_dir')
    parser.add_argument('-p', '--params', dest='params_file')
    parser.add_argument('ssm_file')
    parser.add_argument('cnv_file')
    return parser


def parse_args():
    parser = create_argparser_with_all_args()
    return parser.parse_args()


def print_help():
    parser = create_argparser_with_all_args()
    parser.print_help()
    sys.exit()


def run(safe_to_exit, run_succeeded, config):
    orig_working_dir = switch_working_dir()
    state_manager = StateManager()
    backup_manager = BackupManager(
        [StateManager.default_last_state_fn, TreeWriter.default_archive_fn]
    )

    if state_manager.state_exists():
        logmsg('Resuming existing run. Ignoring command-line parameters (except --output-dir).')
        resume_existing_run(state_manager, backup_manager, safe_to_exit, run_succeeded, config)
    else:
        args = parse_args()
        paths = {
            'ssm_file': args.ssm_file,
            'cnv_file': args.cnv_file,
            'params_file': args.params_file,
            'tmp_dir': args.tmp_dir,
        }
        for K, V in paths.items():
            if V is not None:
                paths[K] = os.path.normpath(os.path.join(orig_working_dir, V))

        try:
            ssm_file = open(paths['ssm_file'])
            cnv_file = open(paths['cnv_file'])
            ssm_file.close()
            cnv_file.close()
        except IOError as e:
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        start_new_run(
            state_manager,
            backup_manager,
            safe_to_exit,
            run_succeeded,
            config,
            paths['ssm_file'],
            paths['cnv_file'],
            paths['params_file'],
            burnin_samples=args.burnin_samples,
            num_samples=args.mcmc_samples,
            mh_itr=args.mh_iterations,
            mh_std=100,
            write_state_every=args.write_state_every,
            write_backups_every=args.write_backups_every,
            rand_seed=args.random_seed,
            tmp_dir=paths['tmp_dir'],
        )


def remove_tmp_files(tmp_dir):
    if tmp_dir is None:
        return
    tmp_filenames = get_c_fnames(tmp_dir)
    for tmpfn in tmp_filenames:
        try:
            os.remove(tmpfn)
        except OSError:
            pass
    try:
        os.rmdir(tmp_dir)
    except OSError:
        pass


def main():
    if '-h' in sys.argv or '--help' in sys.argv:
        print_help()

    safe_to_exit = threading.Event()
    run_succeeded = threading.Event()
    config = {'tmp_dir': None}

    def sigterm_handler(_signo, _stack_frame):
        logmsg('Signal %s received.' % _signo, sys.stderr)
        safe_to_exit.wait()
        remove_tmp_files(config['tmp_dir'])
        logmsg('Exiting now.')
        sys.exit(3)

    signal.signal(signal.SIGTERM, sigterm_handler)
    signal.signal(signal.SIGINT, sigterm_handler)

    run_thread = threading.Thread(target=run, args=(safe_to_exit, run_succeeded, config))
    run_thread.daemon = True
    run_thread.start()

    while True:
        if not run_thread.is_alive():
            break
        run_thread.join(10)

    remove_tmp_files(config['tmp_dir'])
    if run_succeeded.is_set():
        logmsg('Run succeeded.')
        sys.exit(0)
    else:
        logmsg('Run failed.')
        sys.exit(1)


if __name__ == '__main__':
    main()
