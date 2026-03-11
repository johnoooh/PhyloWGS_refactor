#!/usr/bin/env python3
"""Generate runtime plot for PhyloWGS benchmark series."""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

with open('benchmark_runs/results.json') as f:
    runs = json.load(f)

# Sort by total_iterations ascending
runs.sort(key=lambda r: r['total_iterations'])

labels = [r['label'] for r in runs]
x = [r['total_iterations'] for r in runs]
wall_times = [r['wall_time_s'] for r in runs]
sps = [r['samples_per_sec_per_chain'] for r in runs]
is_prod = [r['is_production'] for r in runs]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
fig.suptitle(
    'PhyloWGS Optimized — Runtime vs Iteration Count\n(11 SSMs, max 6 chains)',
    fontsize=14, fontweight='bold'
)

colors = ['#e74c3c' if p else '#2980b9' for p in is_prod]
sizes  = [160 if p else 80 for p in is_prod]
markers = ['*' if p else 'o' for p in is_prod]

# --- Panel 1: wall time ---
ax1.plot(x, wall_times, color='#7f8c8d', linewidth=1.2, zorder=1, linestyle='--')
for xi, yi, c, s, mk, lbl, prod in zip(x, wall_times, colors, sizes, markers, labels, is_prod):
    ax1.scatter(xi, yi, color=c, s=s, zorder=3, marker=mk)
    offset_y = 8 if not prod else 12
    ax1.annotate(lbl, (xi, yi), textcoords='offset points',
                 xytext=(5, offset_y), fontsize=8,
                 color='#c0392b' if prod else '#2c3e50', fontweight='bold' if prod else 'normal')

ax1.set_ylabel('Wall time (seconds)', fontsize=11)
ax1.set_yscale('log')
ax1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax1.grid(True, alpha=0.3, which='both')
ax1.set_title('Wall Clock Time', fontsize=11)

prod_patch = mpatches.Patch(color='#e74c3c', label='production-6c (★)')
other_patch = mpatches.Patch(color='#2980b9', label='other configs')
ax1.legend(handles=[prod_patch, other_patch], fontsize=9)

# --- Panel 2: samples/sec per chain ---
ax2.plot(x, sps, color='#7f8c8d', linewidth=1.2, zorder=1, linestyle='--')
for xi, yi, c, s, mk, lbl, prod in zip(x, sps, colors, sizes, markers, labels, is_prod):
    ax2.scatter(xi, yi, color=c, s=s, zorder=3, marker=mk)
    ax2.annotate(lbl, (xi, yi), textcoords='offset points',
                 xytext=(5, 6), fontsize=8,
                 color='#c0392b' if prod else '#2c3e50', fontweight='bold' if prod else 'normal')

ax2.set_xlabel('Total iterations = (burnin + mcmc_samples) × num_chains', fontsize=11)
ax2.set_ylabel('Samples / sec / chain', fontsize=11)
ax2.set_xscale('log')
ax2.grid(True, alpha=0.3, which='both')
ax2.set_title('Throughput per Chain', fontsize=11)
ax2.axhline(y=np.mean(sps), color='gray', linestyle=':', linewidth=1,
            label=f'mean = {np.mean(sps):.1f} samp/s/chain')
ax2.legend(fontsize=9)

plt.tight_layout()
plt.savefig('runtime_plot.png', dpi=150, bbox_inches='tight')
print("Saved runtime_plot.png")
