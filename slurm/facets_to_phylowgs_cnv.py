#!/usr/bin/env python3
"""
facets_to_phylowgs_cnv.py — Convert FACETS hisens.cncf.txt to the CNV format
expected by create_phylowgs_inputs.py (morrislab/MSK version).

Expected output columns (tab-delimited):
    chromosome  start  end  major_cn  minor_cn  cellular_prevalence

Usage:
    python facets_to_phylowgs_cnv.py input.cncf.txt > cnv_for_phylowgs.txt
    python facets_to_phylowgs_cnv.py input.cncf.txt -o cnv_for_phylowgs.txt
"""

import argparse
import csv
import sys


def convert(input_path, output_fh):
    writer = csv.writer(output_fh, delimiter="\t")
    writer.writerow(["chromosome", "start", "end", "major_cn", "minor_cn", "cellular_prevalence"])

    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # FACETS columns: chrom, loc.start, loc.end, tcn.em, lcn.em, cf.em
            chrom = row["chrom"]
            start = row["loc.start"]
            end = row["loc.end"]

            tcn_raw = row.get("tcn.em", row.get("tcn", ""))
            lcn_raw = row.get("lcn.em", row.get("lcn", ""))
            cf_raw = row.get("cf.em", row.get("cf", ""))

            # Skip rows with missing CN data
            if tcn_raw in ("NA", "", "NaN") or lcn_raw in ("NA", "", "NaN"):
                continue

            tcn = int(float(tcn_raw))
            lcn = int(float(lcn_raw))
            major_cn = tcn - lcn

            # Skip diploid (normal CN)
            if major_cn == 1 and lcn == 1:
                continue

            # Cellular prevalence from cf.em (cellular fraction)
            if cf_raw in ("NA", "", "NaN"):
                cf = 1.0  # assume clonal if missing
            else:
                cf = float(cf_raw)

            writer.writerow([chrom, start, end, major_cn, lcn, cf])


def main():
    parser = argparse.ArgumentParser(description="Convert FACETS cncf.txt to PhyloWGS CNV format")
    parser.add_argument("input", help="FACETS hisens.cncf.txt file")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    args = parser.parse_args()

    if args.output:
        with open(args.output, "w", newline="") as out:
            convert(args.input, out)
    else:
        convert(args.input, sys.stdout)


if __name__ == "__main__":
    main()
