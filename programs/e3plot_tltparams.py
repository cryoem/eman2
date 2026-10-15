#!/usr/bin/env python
# Ported from /home5/pbaldwin/PlotEulersFromHDFfile.ipynb, 2026-08-07
#
# Plots tilt/shift parameters (tx, ty, ztilt, ytilt, xtilt) from a pipeline
# JSON file's "tlt_params" array. Works on both the SA-level JSON (e.g.
# info/<basename>_recon_tomo_final.json) and the full-frame JSON (e.g.
# info/<basename>_Full_Sh4.json) -- run once per file to get diagnostics
# for both the SA run and the large/full run.
#
# Usage:
#   e3plot_tltparams.py info/<basename>_recon_tomo_final.json
#   e3plot_tltparams.py info/<basename>_Full_Sh4.json --outdir diagnostics

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

TITLES = ['tx', 'ty', 'ztilt', 'ytilt', 'xtilt']
DELTA_TITLES = ['d' + t for t in TITLES]


def main():
    parser = argparse.ArgumentParser(
        usage="e3plot_tltparams.py <json_file> [options]",
        description="Plot tilt/shift parameters (tx,ty,ztilt,ytilt,xtilt) from a "
                     "pipeline JSON's tlt_params array. Run on both the SA-level "
                     "and full-frame JSON files to diagnose each stage.")
    parser.add_argument("json_file", type=str,
        help="Path to a JSON file with a 'tlt_params' key (e.g. "
             "info/<basename>_recon_tomo_final.json or info/<basename>_Full_Sh4.json).")
    parser.add_argument("--outdir", type=str, default=".",
        help="Directory to write output jpgs to (default: current directory).")
    parser.add_argument("--label", type=str, default=None,
        help="Override the dataset label used in filenames/titles "
             "(default: derived from json_file's basename).")
    options = parser.parse_args()

    with open(options.json_file) as f:
        data = json.load(f)

    if "tlt_params" not in data:
        raise SystemExit(f"Error: {options.json_file} has no 'tlt_params' key.")

    tlt_params = np.array(data["tlt_params"])
    n_tilts = tlt_params.shape[0]

    label = options.label or os.path.splitext(os.path.basename(options.json_file))[0]
    os.makedirs(options.outdir, exist_ok=True)

    print(f"{options.json_file}: {n_tilts} tilts, label='{label}'")

    # 1. Histogram panel
    plt.figure(figsize=(14, 3))
    for j, title in enumerate(TITLES):
        plt.subplot(1, 5, j + 1)
        plt.hist(tlt_params[:, j], 12)
        plt.title(title)
    plt.suptitle(f"tlt_params histograms — {label}  ({n_tilts} tilts)")
    plt.tight_layout()
    hist_out = os.path.join(options.outdir, f"{label}_tltParams_hist.jpg")
    plt.savefig(hist_out)
    plt.close()
    print(f"  wrote {hist_out}")

    # 2. Sequential (value vs tilt index) panel
    plt.figure(figsize=(14, 2.5))
    for j, title in enumerate(TITLES):
        ax = plt.subplot(1, 5, j + 1)
        plt.plot(tlt_params[:, j], '.')
        plt.xlabel("Tilt index")
        plt.title(title)
        if title == "ytilt":
            ax.yaxis.set_major_locator(MaxNLocator(nbins=15))
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
            ax.grid(True, alpha=0.3)
    plt.suptitle(f"tlt_params vs tilt index — {label}  ({n_tilts} tilts)")
    plt.tight_layout()
    seq_out = os.path.join(options.outdir, f"{label}_tltParams_sequential.jpg")
    plt.savefig(seq_out)
    plt.close()
    print(f"  wrote {seq_out}")

    # 3. Frame-to-frame delta panel
    if n_tilts > 1:
        diff_tlt_params = tlt_params[1:, :] - tlt_params[:-1, :]
        plt.figure(figsize=(14, 2.5))
        for j, title in enumerate(DELTA_TITLES):
            plt.subplot(1, 5, j + 1)
            plt.plot(diff_tlt_params[:, j], '.')
            plt.xlabel("Tilt index")
            plt.title(title)
        plt.suptitle(f"tlt_params frame-to-frame delta — {label}  ({n_tilts} tilts)")
        plt.tight_layout()
        delta_out = os.path.join(options.outdir, f"{label}_tltParams_delta.jpg")
        plt.savefig(delta_out)
        plt.close()
        print(f"  wrote {delta_out}")


if __name__ == "__main__":
    main()
