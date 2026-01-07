# /// script
# dependencies = [
#   "matplotlib>=3.10.3",
#   "numpy>=2.3.2",
#   "polars>=1.31.0"
# ]
# ///

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

def plot_avalanche_boxplots_from_separate_files(csv_files):
    """
    Create two boxplots: one for average bit flips and one for standard deviations
    Each CSV contains avg_bit_flips (column 1) and std_dev (column 2) with no headers
    """
    
    avg_data = []
    std_data = []
    labels = []
    
    for csv_file in csv_files:
        # Extract algorithm name from filename
        algo_name = csv_file.split('/')[-1].replace('.csv', '').replace('avalanche-', '')
        
        # Read the data (no headers)
        df = pl.read_csv(csv_file, has_header=False)
        
        # First column is avg_bit_flips, second is std_dev
        avg_bit_flips = df.get_column('column_1').to_list()
        std_devs = df.get_column('column_2').to_list()
        
        avg_data.append(avg_bit_flips)
        std_data.append(std_devs)
        labels.append(algo_name)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Colors for consistency between plots
    colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow', 'lightpink']
    
    # First boxplot: Average bit flips
    box_plot1 = ax1.boxplot(avg_data, labels=labels, patch_artist=True)
    for patch, color in zip(box_plot1['boxes'], colors[:len(labels)]):
        patch.set_facecolor(color)
    
    ax1.set_xlabel('Hash Algorithm')
    ax1.set_ylabel('Average Number of Bit Flips')
    ax1.set_title('Average Bit Flips Distribution\nAcross Hash Algorithms (1000 samples each)')
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='x', rotation=45)
    
    # Second boxplot: Standard deviations
    box_plot2 = ax2.boxplot(std_data, labels=labels, patch_artist=True)
    for patch, color in zip(box_plot2['boxes'], colors[:len(labels)]):
        patch.set_facecolor(color)
    
    ax2.set_xlabel('Hash Algorithm')
    ax2.set_ylabel('Standard Deviation of Bit Flips')
    ax2.set_title('Standard Deviation Distribution\nAcross Hash Algorithms (1000 samples each)')
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('out/avalanche_boxplots.png', dpi=300, bbox_inches='tight')
    plt.show()

# Main execution with the algorithms from your original script
csv_files = [
    'out/avalanche-jam.csv', 
    'out/avalanche-murmur.csv', 
    'out/avalanche-xxhash.csv'
]

plot_avalanche_boxplots_from_separate_files(csv_files)