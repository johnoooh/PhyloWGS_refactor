#!/usr/bin/env python3
"""
facets_to_phylowgs_cnv.py — Convert FACETS gene-level or segment-level output
to the CNV format expected by create_phylowgs_inputs.py (morrislab/MSK version).

Supports two input formats:
  1. Gene-level (*.gene_level.txt) — deduplicates by segment
  2. Segment-level (*.cncf.txt) — direct conversion

Expected output columns (tab-delimited):
    chromosome  start  end  major_cn  minor_cn  cellular_prevalence

Usage:
    python facets_to_phylowgs_cnv.py input.gene_level.txt -o cnv.txt
    python facets_to_phylowgs_cnv.py input.cncf.txt -o cnv.txt
"""

import argparse
import csv
import sys


def detect_format(header):
    """Detect whether input is gene-level or segment-level FACETS output."""
    if "gene" in header and "seg_start" in header:
        return "gene_level"
    elif "loc.start" in header:
        return "cncf"
    else:
        raise ValueError(
            "Unrecognized FACETS format. Expected gene-level (gene, seg_start, seg_end) "
            "or cncf (loc.start, loc.end) columns. Got: %s" % ", ".join(header)
        )


def parse_cn_fields(row, fmt):
    """Extract chrom, start, end, tcn, lcn, cf from a row based on format."""
    if fmt == "gene_level":
        chrom = row["chrom"]
        start = row["seg_start"]
        end = row["seg_end"]
    else:  # cncf
        chrom = row["chrom"]
        start = row["loc.start"]
        end = row["loc.end"]

    tcn_raw = row.get("tcn.em", row.get("tcn", ""))
    lcn_raw = row.get("lcn.em", row.get("lcn", ""))
    cf_raw = row.get("cf.em", row.get("cf", ""))

    return chrom, start, end, tcn_raw, lcn_raw, cf_raw


def convert(input_path, output_fh):
    writer = csv.writer(output_fh, delimiter="\t")
    writer.writerow(["chromosome", "start", "end", "major_cn", "minor_cn", "cellular_prevalence"])

    seen_segments = set()  # (chrom, start, end) for dedup

    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fmt = detect_format(reader.fieldnames)

        for row in reader:
            # Gene-level files may have a filter column — skip filtered segments
            if fmt == "gene_level" and row.get("filter", "") not in ("", "PASS"):
                # Some FACETS gene-level files use empty string for PASS
                # Skip rows explicitly marked as filtered
                pass  # Don't skip — filter column semantics vary; include all for now

            chrom, start, end, tcn_raw, lcn_raw, cf_raw = parse_cn_fields(row, fmt)

            # Skip rows with missing CN data
            if tcn_raw in ("NA", "", "NaN") or lcn_raw in ("NA", "", "NaN"):
                continue

            tcn = int(float(tcn_raw))
            lcn = int(float(lcn_raw))
            major_cn = tcn - lcn

            # Sanity checks
            if major_cn < 0 or lcn < 0:
                continue

            # Deduplicate segments (gene-level has many genes per segment)
            seg_key = (chrom, start, end)
            if seg_key in seen_segments:
                continue
            seen_segments.add(seg_key)

            # Cellular prevalence from cf.em
            if cf_raw in ("NA", "", "NaN"):
                cf = 1.0  # assume clonal if missing
            else:
                cf = float(cf_raw)

            writer.writerow([chrom, start, end, major_cn, lcn, cf])


def main():
    parser = argparse.ArgumentParser(
        description="Convert FACETS gene-level or cncf output to PhyloWGS CNV format"
    )
    parser.add_argument("input", help="FACETS gene_level.txt or cncf.txt file")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    args = parser.parse_args()

    if args.output:
        with open(args.output, "w", newline="") as out:
            convert(args.input, out)
    else:
        convert(args.input, sys.stdout)


if __name__ == "__main__":
    main()
