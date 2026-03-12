#!/usr/bin/env python3
"""
analyze_benchmark.py — Compare PhyloWGS implementation results across samples

Runs as the terminal SLURM job after all PhyloWGS runs complete.
Reads timing.json + summary.json from each results/<sample>/<impl>/ directory.
Outputs:
  - report.html     — visual comparison
  - summary.tsv     — machine-readable table
  - concordance.tsv — LLH + tree concordance per sample pair

Usage:
    python analyze_benchmark.py --results-dir /path/to/results \
        --outdir /path/to/analysis \
        --implementations optimized-python go-cpu go-cpu-opt go-gpu
"""

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# ── Data loading ───────────────────────────────────────────────────────────────

def load_timing(result_dir: Path) -> dict | None:
    p = result_dir / "timing.json"
    if not p.exists():
        return None
    with open(p) as f:
        return json.load(f)


def load_summary(result_dir: Path) -> dict | None:
    """Load PhyloWGS summary.json output."""
    for fname in ("summary.json", "results.json"):
        p = result_dir / fname
        if p.exists():
            with open(p) as f:
                return json.load(f)
    return None


def load_chain_llh(result_dir: Path) -> list[float]:
    """Extract per-sample LLH values from chain output files."""
    llhs = []
    for chain_file in sorted(result_dir.glob("chain_*_samples.txt")):
        with open(chain_file) as f:
            for line in f:
                try:
                    parts = line.strip().split('\t')
                    llh = float(parts[1])
                    llhs.append(llh)
                except (ValueError, IndexError):
                    continue
    # Also try Go JSON format
    if not llhs:
        for summary_file in result_dir.glob("*.json"):
            try:
                with open(summary_file) as f:
                    data = json.load(f)
                if isinstance(data, dict) and 'best_llh' in data:
                    llhs.append(float(data['best_llh']))
            except Exception:
                pass
    return llhs


def collect_results(results_dir: Path, implementations: list[str]) -> dict:
    """
    Returns:
        {
          sample_id: {
            impl: {
              timing: {...},
              best_llh: float,
              median_llh: float,
              llh_samples: [...],
              n_ssms: int,
              status: 'ok' | 'missing' | 'failed'
            }
          }
        }
    """
    data = defaultdict(dict)
    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sid = sample_dir.name
        for impl in implementations:
            impl_dir = sample_dir / impl
            if not impl_dir.exists():
                data[sid][impl] = {'status': 'missing'}
                continue

            timing = load_timing(impl_dir)
            if timing and timing.get('exit_code', 1) != 0:
                data[sid][impl] = {'status': 'failed', 'timing': timing}
                continue

            llhs = load_chain_llh(impl_dir)
            summary = load_summary(impl_dir)

            best_llh = max(llhs) if llhs else None
            median_llh = float(np.median(llhs)) if llhs else None

            # Try to get best LLH from summary
            if summary:
                best_llh = summary.get('best_llh', best_llh)
                median_llh = summary.get('median_llh', median_llh)

            data[sid][impl] = {
                'status': 'ok',
                'timing': timing,
                'wall_seconds': timing['wall_seconds'] if timing else None,
                'best_llh': best_llh,
                'median_llh': median_llh,
                'llh_samples': llhs,
                'n_llh_samples': len(llhs),
            }

    return dict(data)


# ── Analysis ───────────────────────────────────────────────────────────────────

def compute_speedup(data: dict, implementations: list[str], reference='optimized-python') -> dict:
    """Compute per-sample, per-implementation speedup vs reference."""
    speedups = defaultdict(list)
    for sid, impls in data.items():
        ref = impls.get(reference, {})
        ref_time = ref.get('wall_seconds')
        if not ref_time:
            continue
        for impl in implementations:
            if impl == reference:
                continue
            t = impls.get(impl, {}).get('wall_seconds')
            if t and t > 0:
                speedups[impl].append(ref_time / t)
    return {impl: np.array(v) for impl, v in speedups.items()}


def compute_llh_concordance(data: dict, implementations: list[str]) -> dict:
    """
    Check that LLH ranges overlap across implementations.
    Returns per-sample concordance scores.
    """
    concordance = {}
    for sid, impls in data.items():
        llh_ranges = {}
        for impl in implementations:
            rec = impls.get(impl, {})
            if rec.get('status') == 'ok' and rec.get('llh_samples'):
                llhs = np.array(rec['llh_samples'])
                llh_ranges[impl] = (np.min(llhs), np.max(llhs), np.median(llhs))

        # Check that best LLH values are within 5% relative error across impls
        best_llhs = {
            impl: impls[impl].get('best_llh')
            for impl in implementations
            if impls.get(impl, {}).get('status') == 'ok'
               and impls[impl].get('best_llh') is not None
        }

        if len(best_llhs) < 2:
            concordance[sid] = {'status': 'insufficient_data', 'ranges': llh_ranges}
            continue

        vals = list(best_llhs.values())
        ref_val = vals[0]
        max_rel_diff = max(abs(v - ref_val) / max(abs(ref_val), 1e-10) for v in vals)
        concordance[sid] = {
            'status': 'ok',
            'max_rel_diff': max_rel_diff,
            'concordant': max_rel_diff < 0.05,  # within 5%
            'best_llhs': best_llhs,
            'ranges': llh_ranges,
        }
    return concordance


# ── Reporting ─────────────────────────────────────────────────────────────────

def write_summary_tsv(data: dict, implementations: list[str], out_path: Path):
    with open(out_path, 'w') as f:
        header = ['sample_id', 'impl', 'status', 'wall_seconds', 'best_llh',
                  'median_llh', 'n_llh_samples']
        f.write('\t'.join(header) + '\n')
        for sid in sorted(data.keys()):
            for impl in implementations:
                rec = data[sid].get(impl, {'status': 'missing'})
                row = [
                    sid, impl,
                    rec.get('status', 'missing'),
                    str(rec.get('wall_seconds', 'NA')),
                    str(rec.get('best_llh', 'NA')),
                    str(rec.get('median_llh', 'NA')),
                    str(rec.get('n_llh_samples', 0)),
                ]
                f.write('\t'.join(row) + '\n')


def write_concordance_tsv(concordance: dict, out_path: Path):
    with open(out_path, 'w') as f:
        f.write('sample_id\tconcordant\tmax_rel_llh_diff\t' +
                '\t'.join(f'best_llh_{impl}' for impl in ['optimized-python', 'go-cpu', 'go-cpu-opt', 'go-gpu']) + '\n')
        for sid, rec in sorted(concordance.items()):
            if rec['status'] != 'ok':
                f.write(f"{sid}\tNA\tNA\n")
                continue
            llhs = rec.get('best_llhs', {})
            row = [
                sid,
                str(rec['concordant']),
                f"{rec['max_rel_diff']:.6f}",
            ] + [str(llhs.get(impl, 'NA')) for impl in
                 ['optimized-python', 'go-cpu', 'go-cpu-opt', 'go-gpu']]
            f.write('\t'.join(row) + '\n')


def write_html_report(data: dict, speedups: dict, concordance: dict,
                      implementations: list[str], out_path: Path):
    """Generate a minimal HTML report."""
    n_samples = len(data)
    n_ok = sum(
        1 for impls in data.values()
        if any(v.get('status') == 'ok' for v in impls.values())
    )
    n_concordant = sum(
        1 for c in concordance.values()
        if c.get('concordant', False)
    )

    speedup_rows = ""
    for impl, sups in speedups.items():
        if len(sups):
            speedup_rows += f"<tr><td>{impl}</td><td>{np.mean(sups):.2f}×</td>" \
                           f"<td>{np.median(sups):.2f}×</td>" \
                           f"<td>{np.min(sups):.2f}×</td>" \
                           f"<td>{np.max(sups):.2f}×</td></tr>\n"

    timing_rows = ""
    for sid in sorted(data.keys()):
        timing_rows += f"<tr><td>{sid}</td>"
        for impl in implementations:
            rec = data[sid].get(impl, {})
            t = rec.get('wall_seconds', 'NA')
            status = rec.get('status', 'missing')
            color = '#d4edda' if status == 'ok' else '#f8d7da'
            timing_rows += f"<td style='background:{color}'>{t}</td>"
        timing_rows += "</tr>\n"

    html = f"""<!DOCTYPE html>
<html>
<head>
<title>PhyloWGS Benchmark Report</title>
<style>
body {{ font-family: monospace; margin: 2em; }}
table {{ border-collapse: collapse; margin: 1em 0; }}
th, td {{ border: 1px solid #ccc; padding: 6px 12px; text-align: right; }}
th {{ background: #f0f0f0; }}
td:first-child {{ text-align: left; }}
h2 {{ margin-top: 2em; }}
.stat {{ font-size: 1.5em; font-weight: bold; }}
</style>
</head>
<body>
<h1>PhyloWGS Benchmark Report</h1>
<p>Samples: <span class="stat">{n_samples}</span> &nbsp;|&nbsp;
   Completed: <span class="stat">{n_ok}</span> &nbsp;|&nbsp;
   Concordant: <span class="stat">{n_concordant}/{n_samples}</span></p>

<h2>Speedup vs optimized-python</h2>
<table>
<tr><th>Implementation</th><th>Mean</th><th>Median</th><th>Min</th><th>Max</th></tr>
{speedup_rows}
</table>

<h2>Wall Time per Sample (seconds)</h2>
<table>
<tr><th>Sample</th>{''.join(f'<th>{impl}</th>' for impl in implementations)}</tr>
{timing_rows}
</table>

<h2>LLH Concordance</h2>
<p>Max relative LLH difference across implementations (threshold: 5%)</p>
<table>
<tr><th>Sample</th><th>Concordant</th><th>Max rel diff</th>
{''.join(f'<th>best_llh {impl}</th>' for impl in implementations)}</tr>
{''.join(
    f"<tr><td>{sid}</td>"
    f"<td>{'✅' if c.get('concordant') else '❌'}</td>"
    f"<td>{c.get('max_rel_diff', 'NA'):.4f}</td>"
    + ''.join(f"<td>{c.get('best_llhs', {{}}).get(impl, 'NA'):.1f}</td>" for impl in implementations)
    + "</tr>"
    for sid, c in sorted(concordance.items()) if c.get('status') == 'ok'
)}
</table>
</body>
</html>"""

    with open(out_path, 'w') as f:
        f.write(html)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results-dir', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--implementations', nargs='+',
                        default=['optimized-python', 'go-cpu', 'go-cpu-opt', 'go-gpu'])
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[analyze] Loading results from {args.results_dir}")
    data = collect_results(args.results_dir, args.implementations)
    print(f"[analyze] {len(data)} samples found")

    speedups = compute_speedup(data, args.implementations)
    concordance = compute_llh_concordance(data, args.implementations)

    # Print summary to stdout (captured in SLURM log)
    print("\n=== Speedup Summary (vs optimized-python) ===")
    for impl, sups in speedups.items():
        if len(sups):
            print(f"  {impl:25s}: {np.mean(sups):.2f}× mean, {np.median(sups):.2f}× median "
                  f"(n={len(sups)})")

    n_concordant = sum(1 for c in concordance.values() if c.get('concordant'))
    print(f"\n=== LLH Concordance: {n_concordant}/{len(concordance)} samples ===")
    for sid, c in sorted(concordance.items()):
        if c.get('status') == 'ok':
            status = "✓" if c['concordant'] else "✗"
            print(f"  {status} {sid}: max_rel_diff={c['max_rel_diff']:.4f}")

    # Write outputs
    write_summary_tsv(data, args.implementations, args.outdir / 'summary.tsv')
    write_concordance_tsv(concordance, args.outdir / 'concordance.tsv')
    write_html_report(data, speedups, concordance, args.implementations,
                      args.outdir / 'report.html')

    print(f"\n[analyze] Outputs written to {args.outdir}")
    print(f"  summary.tsv       — per-sample × per-impl timing + LLH")
    print(f"  concordance.tsv   — LLH agreement across implementations")
    print(f"  report.html       — visual summary")


if __name__ == '__main__':
    main()
