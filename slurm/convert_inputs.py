#!/usr/bin/env python3
"""
convert_inputs.py — Convert MAF + FACETS .cncf to PhyloWGS ssm_data.txt + cnv_data.txt
Usage:
    python convert_inputs.py --sample-id SAMPLE --maf /path/to.maf \
        --facets /path/to/facets_hisens.cncf.txt --outdir /path/to/output/
"""

import argparse
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict

# ── Constants ─────────────────────────────────────────────────────────────────
EPSILON = 1e-6

# ── MAF → SSM ─────────────────────────────────────────────────────────────────

def maf_to_ssm(maf_path, sample_id, out_path):
    """Parse MAF and write ssm_data.txt for PhyloWGS."""
    ssms = []
    with open(maf_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('Hugo'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 50:
                continue
            # Standard MAF columns
            try:
                chrom = cols[4].lstrip('chr')
                pos = int(cols[5])
                ref = cols[10]
                alt = cols[12]
                variant_type = cols[9]

                # Skip indels and non-SNPs if desired
                if variant_type not in ('SNP', 'DNP', 'TNP'):
                    continue

                # t_ref_count, t_alt_count — columns vary by MAF flavor
                # Try common positions
                t_ref = None
                t_alt = None
                for header_guess in [('t_ref_count', 't_alt_count'), ('i_t_ref_count', 'i_t_alt_count')]:
                    pass

                # Fall back to column indices (common MSK MAF format)
                # t_ref_count = col 39, t_alt_count = col 40 in standard MAF
                t_ref = int(cols[39]) if cols[39].isdigit() else None
                t_alt = int(cols[40]) if cols[40].isdigit() else None

                if t_ref is None or t_alt is None or (t_ref + t_alt) == 0:
                    continue

                total = t_ref + t_alt
                gene = cols[0]
                ssm_id = f"s{len(ssms)}"
                ssms.append({
                    'id': ssm_id,
                    'gene': gene,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    't_ref': t_ref,
                    't_alt': t_alt,
                    'total': total,
                    'mu_r': 1e-3,   # error rate
                    'mu_v': 0.5,    # expected VAF if present
                })
            except (ValueError, IndexError):
                continue

    if not ssms:
        raise ValueError(f"No valid SSMs parsed from {maf_path}")

    with open(out_path, 'w') as f:
        f.write("id\tgene\ta\td\tmu_r\tmu_v\n")
        for s in ssms:
            # a = ref read count (alt reads NOT supporting variant)
            # d = total depth
            f.write(f"{s['id']}\t{s['gene']}\t{s['t_ref']}\t{s['total']}\t{s['mu_r']}\t{s['mu_v']}\n")

    print(f"[convert] {len(ssms)} SSMs written to {out_path}", file=sys.stderr)
    return ssms


# ── FACETS .cncf → CNV ────────────────────────────────────────────────────────

def facets_to_cnv(facets_path, ssms, out_path, purity=None):
    """
    Parse FACETS hisens.cncf.txt and write cnv_data.txt for PhyloWGS.
    
    FACETS columns: ID chrom loc.start loc.end seg cnlr.median mafR mafR.clust
                    cf.em tcn.em lcn.em
    cf.em = cellular fraction (tumor purity × local CCF)
    tcn.em = total copy number
    lcn.em = lesser (minor) copy number
    
    CNV read counts (a, d) are synthesized to match PhyloWGS conventions:
      d (total reads) = region_size × het_snp_rate × avg_read_depth
      a (ref reads)   = (1 - cellular_prevalence/2) × d
    This follows the original create_phylowgs_inputs.py logic.
    """
    # Estimate average read depth from SSMs
    ssm_depths = [s['total'] for s in ssms]
    avg_read_depth = sum(ssm_depths) / len(ssm_depths) if ssm_depths else 100
    het_snp_rate = 7.0 / 10000  # 7 het SNPs per 10 kb (PhyloWGS default)
    max_equiv_ssms = 3000  # cap per CNA region (PhyloWGS default)

    # Index SSMs by (chrom, pos) for overlap detection
    ssm_index = defaultdict(list)
    for s in ssms:
        ssm_index[s['chrom']].append(s)

    cnvs = []
    with open(facets_path) as f:
        header = f.readline().strip().split('\t')
        col = {h: i for i, h in enumerate(header)}

        for line in f:
            row = line.strip().split('\t')
            if len(row) < len(col):
                continue
            try:
                chrom = row[col['chrom']].lstrip('chr')
                start = int(row[col['loc.start']])
                end = int(row[col['loc.end']])
                tcn = row[col.get('tcn.em', col.get('tcn', -1))]
                lcn = row[col.get('lcn.em', col.get('lcn', -1))]
                cf = row[col.get('cf.em', col.get('cf', -1))]

                if tcn in ('NA', '', 'NaN') or lcn in ('NA', '', 'NaN'):
                    continue

                tcn = int(float(tcn))
                lcn = int(float(lcn))
                major = tcn - lcn

                if tcn == 2 and lcn == 1:
                    continue  # diploid, skip

                cf_val = float(cf) if cf not in ('NA', '', 'NaN') else 0.5
                cf_val = max(EPSILON, min(1.0 - EPSILON, cf_val))

                # Find overlapping SSMs — format: ssm_id,maternal_cn,paternal_cn
                # maternal = major allele copies, paternal = minor allele copies
                overlapping = []
                for s in ssm_index.get(chrom, []):
                    if start <= s['pos'] <= end:
                        overlapping.append(f"{s['id']},{major},{lcn}")

                # Synthesize read counts (matching PhyloWGS create_phylowgs_inputs.py)
                region_size = end - start
                equiv_ssms = region_size * het_snp_rate
                total_reads = int(round(equiv_ssms * avg_read_depth))
                total_reads = max(1, min(int(round(max_equiv_ssms * avg_read_depth)), total_reads))
                ref_reads = int((1 - cf_val / 2) * total_reads)

                cnvs.append({
                    'cnv_id': f"c{len(cnvs)}",
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'major_cn': major,
                    'minor_cn': lcn,
                    'cellular_prevalence': cf_val,
                    'ref_reads': ref_reads,
                    'total_reads': total_reads,
                    'ssms': overlapping,
                })
            except (ValueError, IndexError, KeyError):
                continue

    with open(out_path, 'w') as f:
        f.write("cnv\ta\td\tssms\tphysical_cnvs\n")
        for c in cnvs:
            ssm_str = ';'.join(c['ssms']) if c['ssms'] else ''
            # physical_cnvs: key=value pairs expected by parse_physical_cnvs
            # Multiple physical CNVs per logical CNV separated by ';'
            phys = (
                f"chrom={c['chrom']},"
                f"start={c['start']},"
                f"end={c['end']},"
                f"major_cn={c['major_cn']},"
                f"minor_cn={c['minor_cn']},"
                f"cell_prev={c['cellular_prevalence']:.4f}"
            )
            # a = ref reads, d = total reads (synthesized from region size + read depth)
            f.write(
                f"{c['cnv_id']}\t{c['ref_reads']}\t{c['total_reads']}\t"
                f"{ssm_str}\t{phys}\n"
            )

    print(f"[convert] {len(cnvs)} CNV segments written to {out_path}", file=sys.stderr)
    return cnvs


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Convert MAF + FACETS → PhyloWGS inputs")
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--maf', required=True)
    parser.add_argument('--facets', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ssm_path = os.path.join(args.outdir, 'ssm_data.txt')
    cnv_path = os.path.join(args.outdir, 'cnv_data.txt')

    ssms = maf_to_ssm(args.maf, args.sample_id, ssm_path)
    facets_to_cnv(args.facets, ssms, cnv_path)

    print(f"[convert] Done: {args.sample_id} → {args.outdir}")


if __name__ == '__main__':
    main()
