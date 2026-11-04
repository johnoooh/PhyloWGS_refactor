#!/usr/bin/env python2
"""Run inside the reference container: decode write_results outputs to pretty JSON.

usage: dump_expected.py <outdir containing summ.json.gz muts.json.gz mutass.zip> <expected dir>
Uses json.dumps(sort_keys=True, indent=1) so Python's float literals are kept.
"""
import gzip
import json
import os
import sys
import zipfile


def dump(obj, path):
    with open(path, 'w') as f:
        f.write(json.dumps(obj, sort_keys=True, indent=1, separators=(',', ': ')) + '\n')


def main():
    src, dst = sys.argv[1], sys.argv[2]
    if not os.path.isdir(dst):
        os.makedirs(dst)
    for gz, out in (('summ.json.gz', 'summ.json'), ('muts.json.gz', 'muts.json')):
        dump(json.load(gzip.GzipFile(os.path.join(src, gz))), os.path.join(dst, out))
    mutass = {}
    zf = zipfile.ZipFile(os.path.join(src, 'mutass.zip'))
    for name in zf.namelist():
        mutass[name[:-len('.json')]] = json.loads(zf.read(name))
    dump(mutass, os.path.join(dst, 'mutass.json'))


if __name__ == '__main__':
    main()
