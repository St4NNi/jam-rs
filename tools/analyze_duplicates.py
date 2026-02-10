#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy<2",
#     "polars",
#     "pyarrow",
#     "matplotlib",
#     "seaborn",
#     "scipy",
# ]
# ///
"""
Analyze hash collision distribution from duplicate files.

Supports two formats:
- Old format: duplicates_{bindex}_{bin_id}.bin (hash only, u64)
- New format: duplicates_with_counts_{bindex}_{bin_id}.bin (hash u64 + count u32)
"""

import argparse
import struct
from pathlib import Path
from collections import defaultdict
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy import stats

# Publication-quality plot settings
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (8, 6),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})


def read_duplicates_old_format(filepath: Path) -> list[int]:
    """Read old format: just hash values (u64)."""
    hashes = []
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(8)
            if len(data) < 8:
                break
            hashes.append(struct.unpack('<Q', data)[0])
    return hashes


def read_duplicates_with_counts(filepath: Path) -> list[tuple[int, int]]:
    """Read new format: hash (u64) + count (u32) pairs."""
    duplicates = []
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(12)  # 8 bytes hash + 4 bytes count
            if len(data) < 12:
                break
            hash_val, count = struct.unpack('<QI', data)
            duplicates.append((hash_val, count))
    return duplicates


def load_all_duplicates(input_dir: Path, use_counts: bool = True) -> dict:
    """
    Load all duplicate files from directory.

    Returns dict with structure:
    {
        'by_bin': {(bindex, bin_id): [(hash, count), ...]},
        'all_counts': [count1, count2, ...],  # all multiplicities
        'all_hashes': [hash1, hash2, ...],    # all collision hashes
        'bin_stats': {(bindex, bin_id): {'total': N, 'unique': M}},
    }
    """
    result = {
        'by_bin': defaultdict(list),
        'all_counts': [],
        'all_hashes': [],
        'bin_stats': {},
    }

    # Try new format first, fall back to old
    if use_counts:
        pattern = "duplicates_with_counts_*_*.bin"
    else:
        pattern = "duplicates_*_*.bin"

    files = list(input_dir.glob(pattern))

    # If no new format files found, try old format
    if not files and use_counts:
        print("No duplicates_with_counts_*.bin files found, trying old format...")
        pattern = "duplicates_*_*.bin"
        files = list(input_dir.glob(pattern))
        use_counts = False

    if not files:
        raise FileNotFoundError(f"No duplicate files found in {input_dir}")

    print(f"Found {len(files)} duplicate files")

    for filepath in sorted(files):
        # Parse bindex and bin_id from filename
        parts = filepath.stem.split('_')
        if use_counts:
            # duplicates_with_counts_{bindex}_{bin_id}
            bindex = int(parts[-2])
            bin_id = int(parts[-1])
        else:
            # duplicates_{bindex}_{bin_id}
            bindex = int(parts[-2])
            bin_id = int(parts[-1])

        if use_counts:
            duplicates = read_duplicates_with_counts(filepath)
            for hash_val, count in duplicates:
                result['by_bin'][(bindex, bin_id)].append((hash_val, count))
                result['all_counts'].append(count)
                result['all_hashes'].append(hash_val)
        else:
            # Old format: assume count=2 (we only know it's a duplicate)
            hashes = read_duplicates_old_format(filepath)
            for hash_val in hashes:
                result['by_bin'][(bindex, bin_id)].append((hash_val, 2))
                result['all_counts'].append(2)
                result['all_hashes'].append(hash_val)

        result['bin_stats'][(bindex, bin_id)] = {
            'unique': len(result['by_bin'][(bindex, bin_id)]),
            'total': sum(c for _, c in result['by_bin'][(bindex, bin_id)])
        }

    return result


def plot_multiplicity_histogram(counts: list[int], output_path: Path, hash_name: str = "jamhash"):
    """
    Plot histogram of collision multiplicities.
    Figure for Section 1.1.2: Collision cluster sizes.
    """
    counts = np.array(counts)

    fig, ax = plt.subplots(figsize=(8, 5))

    # Use log bins for better visualization
    if counts.max() > 10:
        bins = np.logspace(np.log10(2), np.log10(counts.max() + 1), 50)
        ax.set_xscale('log')
    else:
        bins = np.arange(2, counts.max() + 2) - 0.5

    ax.hist(counts, bins=bins, edgecolor='black', alpha=0.7, color='#2E86AB')
    ax.set_yscale('log')

    ax.set_xlabel('Collision Multiplicity (k-mers per hash)')
    ax.set_ylabel('Frequency (log scale)')
    ax.set_title(f'{hash_name}: Distribution of Collision Multiplicities')

    # Add statistics text box
    stats_text = (
        f'N = {len(counts):,}\n'
        f'Mean = {counts.mean():.2f}\n'
        f'Median = {np.median(counts):.1f}\n'
        f'Max = {counts.max():,}'
    )
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_bin_distribution_heatmap(data: dict, output_path: Path, hash_name: str = "jamhash"):
    """
    Plot heatmap of collisions across hash space bins.
    Figure for Section 1.1.2: Spatial distribution across hash space.
    """
    # Create 16x128 matrix (bindex x bin_id)
    matrix = np.zeros((16, 128))

    for (bindex, bin_id), stats in data['bin_stats'].items():
        if bindex < 16 and bin_id < 128:
            matrix[bindex, bin_id] = stats['unique']

    fig, ax = plt.subplots(figsize=(14, 6))

    # Use log scale for better visualization if there's large variation
    if matrix.max() > 0:
        im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd',
                       norm=plt.matplotlib.colors.LogNorm(vmin=max(1, matrix[matrix > 0].min()),
                                                          vmax=matrix.max()))
    else:
        im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd')

    ax.set_xlabel('Bin ID (bits 53-59)')
    ax.set_ylabel('Bindex (bits 60-63)')
    ax.set_title(f'{hash_name}: Collision Distribution Across Hash Space')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Unique Collisions (log scale)')

    # Set tick marks
    ax.set_yticks(range(16))
    ax.set_xticks(range(0, 128, 16))

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_collision_rate_by_bin(data: dict, output_path: Path, hash_name: str = "jamhash"):
    """
    Plot collision rate per bin as a bar chart.
    Shows variation in collision density across hash space regions.
    """
    # Aggregate by bindex
    bindex_stats = defaultdict(lambda: {'unique': 0, 'bins_with_collisions': 0})

    for (bindex, bin_id), stats in data['bin_stats'].items():
        bindex_stats[bindex]['unique'] += stats['unique']
        if stats['unique'] > 0:
            bindex_stats[bindex]['bins_with_collisions'] += 1

    bindex_vals = sorted(bindex_stats.keys())
    collision_counts = [bindex_stats[b]['unique'] for b in bindex_vals]
    bins_with_collisions = [bindex_stats[b]['bins_with_collisions'] for b in bindex_vals]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Total collisions per bindex
    ax1.bar(bindex_vals, collision_counts, color='#2E86AB', edgecolor='black')
    ax1.set_xlabel('Bindex (hash bits 60-63)')
    ax1.set_ylabel('Total Unique Collisions')
    ax1.set_title(f'{hash_name}: Collisions by Hash Region')
    ax1.set_xticks(range(16))

    # Add mean line
    mean_collisions = np.mean(collision_counts)
    ax1.axhline(y=mean_collisions, color='red', linestyle='--', label=f'Mean: {mean_collisions:.0f}')
    ax1.legend()

    # Plot 2: Bins with collisions per bindex (out of 128)
    ax2.bar(bindex_vals, bins_with_collisions, color='#A23B72', edgecolor='black')
    ax2.set_xlabel('Bindex (hash bits 60-63)')
    ax2.set_ylabel('Bins with Collisions (out of 128)')
    ax2.set_title(f'{hash_name}: Coverage of Collision Bins')
    ax2.set_xticks(range(16))
    ax2.set_ylim(0, 128)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_seaborn_multiplicity_kde(counts: list[int], output_path: Path, hash_name: str = "jamhash"):
    """
    Plot KDE of collision multiplicities using seaborn.
    """
    counts = np.array(counts)

    # Set seaborn style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Histogram with KDE overlay
    ax1 = axes[0]
    sns.histplot(counts, kde=True, ax=ax1, color='#2E86AB', edgecolor='white',
                 stat='density', bins=min(50, len(np.unique(counts))))
    ax1.set_xlabel('Collision Multiplicity')
    ax1.set_ylabel('Density')
    ax1.set_title(f'{hash_name}: Multiplicity Distribution')

    # Plot 2: Box plot with swarm (if not too many points)
    ax2 = axes[1]
    if len(counts) < 1000:
        sns.boxplot(x=counts, ax=ax2, color='#2E86AB', width=0.3)
        sns.swarmplot(x=counts, ax=ax2, color='#A23B72', alpha=0.5, size=3)
    else:
        sns.boxplot(x=counts, ax=ax2, color='#2E86AB', width=0.5)
        sns.stripplot(x=counts, ax=ax2, color='#A23B72', alpha=0.3, size=2, jitter=True)
    ax2.set_xlabel('Collision Multiplicity')
    ax2.set_title(f'{hash_name}: Multiplicity Box Plot')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_seaborn_heatmap(data: dict, output_path: Path, hash_name: str = "jamhash"):
    """
    Plot seaborn heatmap of collisions across hash space bins.
    """
    # Create 16x128 matrix (bindex x bin_id)
    matrix = np.zeros((16, 128))

    for (bindex, bin_id), stats_val in data['bin_stats'].items():
        if bindex < 16 and bin_id < 128:
            matrix[bindex, bin_id] = stats_val['unique']

    # Apply log transform for better visualization (add 1 to handle zeros)
    log_matrix = np.log10(matrix + 1)

    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)

    fig, ax = plt.subplots(figsize=(16, 6))

    sns.heatmap(log_matrix, ax=ax, cmap='YlOrRd', cbar_kws={'label': 'log10(collisions + 1)'},
                xticklabels=16, yticklabels=True)

    ax.set_xlabel('Bin ID (bits 53-59)')
    ax.set_ylabel('Bindex (bits 60-63)')
    ax.set_title(f'{hash_name}: Collision Distribution Across Hash Space (log scale)')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_seaborn_violin(data: dict, output_path: Path, hash_name: str = "jamhash"):
    """
    Plot violin plot of collision counts per bindex region.
    """
    # Create DataFrame for seaborn
    records = []
    for (bindex, bin_id), duplicates in data['by_bin'].items():
        for hash_val, count in duplicates:
            records.append({'bindex': bindex, 'bin_id': bin_id, 'multiplicity': count})

    if not records:
        print("No data for violin plot")
        return

    df = pl.DataFrame(records).to_pandas()

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.1)

    fig, ax = plt.subplots(figsize=(14, 6))

    # Violin plot showing multiplicity distribution per bindex
    sns.violinplot(data=df, x='bindex', y='multiplicity', ax=ax,
                   hue='bindex', palette='viridis', inner='box', cut=0, legend=False)

    ax.set_xlabel('Bindex (hash bits 60-63)')
    ax.set_ylabel('Collision Multiplicity')
    ax.set_title(f'{hash_name}: Multiplicity Distribution by Hash Region')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_seaborn_ecdf(counts: list[int], output_path: Path, hash_name: str = "jamhash"):
    """
    Plot empirical cumulative distribution function of multiplicities.
    """
    counts = np.array(counts)

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)

    fig, ax = plt.subplots(figsize=(8, 6))

    sns.ecdfplot(counts, ax=ax, color='#2E86AB', linewidth=2)

    ax.set_xlabel('Collision Multiplicity')
    ax.set_ylabel('Cumulative Proportion')
    ax.set_title(f'{hash_name}: ECDF of Collision Multiplicities')

    # Add vertical lines for key percentiles
    for pct, color in [(50, 'green'), (90, 'orange'), (99, 'red')]:
        val = np.percentile(counts, pct)
        ax.axvline(x=val, color=color, linestyle='--', alpha=0.7,
                   label=f'{pct}th percentile: {val:.0f}')

    ax.legend(loc='lower right')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_seaborn_pairplot(data: dict, output_path: Path, hash_name: str = "jamhash"):
    """
    Create summary pair plot of bin-level statistics.
    """
    # Create DataFrame of bin statistics
    records = []
    for (bindex, bin_id), stats_val in data['bin_stats'].items():
        records.append({
            'bindex': bindex,
            'bin_id': bin_id,
            'unique_collisions': stats_val['unique'],
            'total_count': stats_val['total'],
            'mean_multiplicity': stats_val['total'] / stats_val['unique'] if stats_val['unique'] > 0 else 0
        })

    if len(records) < 5:
        print("Not enough data for pair plot")
        return

    df = pl.DataFrame(records).to_pandas()

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.0)

    # Create pair plot
    g = sns.pairplot(df[['unique_collisions', 'total_count', 'mean_multiplicity']],
                     diag_kind='kde', plot_kws={'alpha': 0.6, 's': 30},
                     diag_kws={'fill': True})
    g.fig.suptitle(f'{hash_name}: Bin-Level Collision Statistics', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_combined_summary(data: dict, counts: list[int], output_path: Path, hash_name: str = "jamhash"):
    """
    Create a combined 2x2 summary figure for the paper.
    """
    counts = np.array(counts)

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.0)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Histogram with KDE
    ax1 = axes[0, 0]
    sns.histplot(counts, kde=True, ax=ax1, color='#2E86AB', edgecolor='white',
                 stat='count', bins=min(50, len(np.unique(counts))))
    ax1.set_xlabel('Collision Multiplicity')
    ax1.set_ylabel('Count')
    ax1.set_title('A) Multiplicity Distribution')

    # 2. ECDF
    ax2 = axes[0, 1]
    sns.ecdfplot(counts, ax=ax2, color='#2E86AB', linewidth=2)
    ax2.set_xlabel('Collision Multiplicity')
    ax2.set_ylabel('Cumulative Proportion')
    ax2.set_title('B) Cumulative Distribution')
    # Add 90th and 99th percentile lines
    for pct, color in [(90, 'orange'), (99, 'red')]:
        val = np.percentile(counts, pct)
        ax2.axvline(x=val, color=color, linestyle='--', alpha=0.7,
                    label=f'P{pct}: {val:.0f}')
    ax2.legend(loc='lower right', fontsize=8)

    # 3. Heatmap
    ax3 = axes[1, 0]
    matrix = np.zeros((16, 128))
    for (bindex, bin_id), stats_val in data['bin_stats'].items():
        if bindex < 16 and bin_id < 128:
            matrix[bindex, bin_id] = stats_val['unique']
    log_matrix = np.log10(matrix + 1)
    sns.heatmap(log_matrix, ax=ax3, cmap='YlOrRd', cbar_kws={'label': 'log10(n+1)'},
                xticklabels=32, yticklabels=True)
    ax3.set_xlabel('Bin ID')
    ax3.set_ylabel('Bindex')
    ax3.set_title('C) Spatial Distribution')

    # 4. Bar chart by bindex
    ax4 = axes[1, 1]
    bindex_totals = defaultdict(int)
    for (bindex, bin_id), stats_val in data['bin_stats'].items():
        bindex_totals[bindex] += stats_val['unique']
    bindex_vals = sorted(bindex_totals.keys())
    collision_counts = [bindex_totals[b] for b in bindex_vals]
    sns.barplot(x=bindex_vals, y=collision_counts, ax=ax4, hue=bindex_vals, palette='viridis', legend=False)
    ax4.set_xlabel('Bindex (hash bits 60-63)')
    ax4.set_ylabel('Total Collisions')
    ax4.set_title('D) Collisions by Region')

    plt.suptitle(f'{hash_name}: Collision Analysis Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def compute_statistics(data: dict) -> dict:
    """Compute comprehensive statistics for the paper."""
    counts = np.array(data['all_counts'])
    hashes = np.array(data['all_hashes'])

    stats_dict = {
        'total_unique_collisions': len(counts),
        'total_collision_count': int(counts.sum()),
        'multiplicity': {
            'mean': float(counts.mean()),
            'median': float(np.median(counts)),
            'std': float(counts.std()),
            'min': int(counts.min()),
            'max': int(counts.max()),
            'percentile_95': float(np.percentile(counts, 95)),
            'percentile_99': float(np.percentile(counts, 99)),
        },
        'bins': {
            'total_possible': 16 * 128,
            'with_collisions': len(data['bin_stats']),
            'coverage_pct': len(data['bin_stats']) / (16 * 128) * 100,
        },
    }

    # Chi-squared test for uniformity across bins
    if len(data['bin_stats']) > 1:
        observed = [s['unique'] for s in data['bin_stats'].values()]
        expected = np.mean(observed)
        chi2, p_value = stats.chisquare(observed)
        stats_dict['uniformity'] = {
            'chi_squared': float(chi2),
            'p_value': float(p_value),
            'interpretation': 'uniform' if p_value > 0.05 else 'non-uniform'
        }

    # Coefficient of variation across bins with collisions
    bin_counts = [s['unique'] for s in data['bin_stats'].values()]
    if len(bin_counts) > 0:
        cv = np.std(bin_counts) / np.mean(bin_counts) if np.mean(bin_counts) > 0 else 0
        stats_dict['bins']['coefficient_of_variation'] = float(cv)

    return stats_dict


def print_statistics(stats_dict: dict):
    """Print formatted statistics summary."""
    print("\n" + "=" * 60)
    print("COLLISION ANALYSIS SUMMARY")
    print("=" * 60)

    print(f"\n{'Total unique collision hashes:':<40} {stats_dict['total_unique_collisions']:,}")
    print(f"{'Total collision count:':<40} {stats_dict['total_collision_count']:,}")

    print("\n--- Multiplicity Statistics ---")
    m = stats_dict['multiplicity']
    print(f"{'Mean multiplicity:':<40} {m['mean']:.2f}")
    print(f"{'Median multiplicity:':<40} {m['median']:.1f}")
    print(f"{'Std deviation:':<40} {m['std']:.2f}")
    print(f"{'Min multiplicity:':<40} {m['min']}")
    print(f"{'Max multiplicity:':<40} {m['max']}")
    print(f"{'95th percentile:':<40} {m['percentile_95']:.1f}")
    print(f"{'99th percentile:':<40} {m['percentile_99']:.1f}")

    print("\n--- Bin Distribution ---")
    b = stats_dict['bins']
    print(f"{'Total possible bins:':<40} {b['total_possible']}")
    print(f"{'Bins with collisions:':<40} {b['with_collisions']}")
    print(f"{'Coverage:':<40} {b['coverage_pct']:.2f}%")
    if 'coefficient_of_variation' in b:
        print(f"{'Coefficient of variation:':<40} {b['coefficient_of_variation']:.4f}")

    if 'uniformity' in stats_dict:
        print("\n--- Uniformity Test (Chi-squared) ---")
        u = stats_dict['uniformity']
        print(f"{'Chi-squared statistic:':<40} {u['chi_squared']:.2f}")
        print(f"{'P-value:':<40} {u['p_value']:.4e}")
        print(f"{'Interpretation:':<40} {u['interpretation']}")

    print("\n" + "=" * 60)


def save_statistics_csv(stats_dict: dict, output_path: Path):
    """Save statistics to CSV for inclusion in paper."""
    with open(output_path, 'w') as f:
        f.write("metric,value\n")
        f.write(f"total_unique_collisions,{stats_dict['total_unique_collisions']}\n")
        f.write(f"total_collision_count,{stats_dict['total_collision_count']}\n")

        m = stats_dict['multiplicity']
        f.write(f"multiplicity_mean,{m['mean']:.4f}\n")
        f.write(f"multiplicity_median,{m['median']:.1f}\n")
        f.write(f"multiplicity_std,{m['std']:.4f}\n")
        f.write(f"multiplicity_min,{m['min']}\n")
        f.write(f"multiplicity_max,{m['max']}\n")
        f.write(f"multiplicity_p95,{m['percentile_95']:.1f}\n")
        f.write(f"multiplicity_p99,{m['percentile_99']:.1f}\n")

        b = stats_dict['bins']
        f.write(f"bins_total,{b['total_possible']}\n")
        f.write(f"bins_with_collisions,{b['with_collisions']}\n")
        f.write(f"bins_coverage_pct,{b['coverage_pct']:.4f}\n")
        if 'coefficient_of_variation' in b:
            f.write(f"bins_cv,{b['coefficient_of_variation']:.6f}\n")

        if 'uniformity' in stats_dict:
            u = stats_dict['uniformity']
            f.write(f"chi_squared,{u['chi_squared']:.4f}\n")
            f.write(f"chi_squared_p_value,{u['p_value']:.4e}\n")

    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze hash collision distribution for research paper',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s /path/to/duplicates --output ./analysis
  %(prog)s /path/to/duplicates --hash-name murmur3 --no-counts
        """
    )
    parser.add_argument('input_dir', type=Path,
                        help='Directory containing duplicate files')
    parser.add_argument('--output', '-o', type=Path, default=Path('./collision_analysis'),
                        help='Output directory for plots and statistics')
    parser.add_argument('--hash-name', type=str, default='jamhash',
                        help='Name of hash function for plot titles')
    parser.add_argument('--no-counts', action='store_true',
                        help='Use old format files (without counts)')

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading duplicates from {args.input_dir}...")
    data = load_all_duplicates(args.input_dir, use_counts=not args.no_counts)

    if not data['all_counts']:
        print("No collision data found!")
        return

    # Compute statistics
    stats_dict = compute_statistics(data)
    print_statistics(stats_dict)

    # Save statistics
    save_statistics_csv(stats_dict, args.output / f'{args.hash_name}_collision_stats.csv')

    # Generate plots
    print("\nGenerating plots...")

    plot_multiplicity_histogram(
        data['all_counts'],
        args.output / f'{args.hash_name}_multiplicity_histogram.png',
        args.hash_name
    )

    plot_bin_distribution_heatmap(
        data,
        args.output / f'{args.hash_name}_bin_distribution_heatmap.png',
        args.hash_name
    )

    plot_collision_rate_by_bin(
        data,
        args.output / f'{args.hash_name}_collision_by_region.png',
        args.hash_name
    )

    # Seaborn plots
    print("\nGenerating seaborn plots...")

    plot_seaborn_multiplicity_kde(
        data['all_counts'],
        args.output / f'{args.hash_name}_seaborn_kde.png',
        args.hash_name
    )

    plot_seaborn_heatmap(
        data,
        args.output / f'{args.hash_name}_seaborn_heatmap.png',
        args.hash_name
    )

    plot_seaborn_violin(
        data,
        args.output / f'{args.hash_name}_seaborn_violin.png',
        args.hash_name
    )

    plot_seaborn_ecdf(
        data['all_counts'],
        args.output / f'{args.hash_name}_seaborn_ecdf.png',
        args.hash_name
    )

    plot_seaborn_pairplot(
        data,
        args.output / f'{args.hash_name}_seaborn_pairplot.png',
        args.hash_name
    )

    plot_combined_summary(
        data,
        data['all_counts'],
        args.output / f'{args.hash_name}_combined_summary.png',
        args.hash_name
    )

    print(f"\nAll outputs saved to: {args.output}")


if __name__ == '__main__':
    main()
