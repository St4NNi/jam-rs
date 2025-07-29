# Copied (and modified) from FoldHash's avalanche test
# https://github.com/orlp/foldhash/blob/master/tools/plot-avalanche.py

# /// script
# dependencies = [
#   "matplotlib>=3.10.3",
#   "numpy>=2.3.2",
#   "polars>=1.31.0"
# ]
# ///


import matplotlib # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np
import polars as pl # type: ignore

def plot_avalanche(hashname):
    vals = pl.read_csv(f"out/avalanche-{hashname}.csv", has_header=False).to_numpy().reshape((64, 64))
        
    cm = matplotlib.colormaps["viridis"]
    plt.clf()
    plt.imshow(vals, cmap=cm, vmin=0, vmax=1, origin="lower")
    plt.colorbar()
    plt.xlabel("Input bit")
    plt.ylabel("Output bit")
    title = f"Worst-case avalanche diagram of {hashname}"
    plt.title(title)
    plt.savefig(f"out/avalanche-{hashname}.png")

plot_avalanche("jamhash")
plot_avalanche("murmur3")
plot_avalanche("xxhash3")