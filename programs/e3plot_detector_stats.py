#!/usr/bin/env python
# Ported from /home5/pbaldwin/GetStatsOn.ipynb, 2026-08-10
#
# Detector diagnostic: pixel-value histogram and "fill fraction" (fraction of
# pixels registering at least one electron count) for a single interior frame
# of a raw counting-mode movie. Designed to run silently as part of the import
# pipeline (see run_detector_stats() in e3ctomo_import.py) and write jpgs to
# diagnostics/ -- failures here must never abort the pipeline.
#
# The original notebook computed the discrete count histogram with a Python
# loop over every unique pixel value calling np.where() per value -- fine for
# a handful of small integer counts, but O(n * num_unique) and very slow if
# ever run on non-counting-mode data. Replaced with np.bincount(), which is
# the correct vectorized approach for small non-negative integer data and
# lets us just skip the plot (rather than hang) when the data doesn't look
# like raw counts.
#
# Usage:
#   e3plot_detector_stats.py <movie.hdf> [--frame N] [--outdir diagnostics]

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from EMAN3 import EMData, EMUtil, Region


def main():
    parser = argparse.ArgumentParser(
        usage="e3plot_detector_stats.py <movie.hdf> [options]",
        description="Detector pixel-value histogram and fill-fraction diagnostic "
                     "for a single interior frame of a raw counting-mode movie.")
    parser.add_argument("movie", type=str, help="Path to a movie HDF/MRC stack.")
    parser.add_argument("--frame", type=int, default=-1,
        help="Frame index to sample (default: an interior/middle frame).")
    parser.add_argument("--outdir", type=str, default=".",
        help="Directory to write output jpgs to (default: current directory).")
    parser.add_argument("--label", type=str, default=None,
        help="Override the dataset label used in filenames/titles "
             "(default: derived from the movie's basename).")
    options = parser.parse_args()

    nimg = EMUtil.get_image_count(options.movie)
    hdr = EMData(options.movie, 0, True)   # header-only, cheap regardless of file size
    nz = hdr["nz"]

    # Raw MRC movies are often read as a single 3D volume (nimg=1, nz=nframes)
    # rather than a stack of nimg separate 2D images -- a plain EMData(path, i)
    # read would try to allocate the *entire* volume in that case. Use a
    # Region to pull just one Z-slice instead.
    if nimg == 1 and nz > 1:
        nframes = nz
        nx, ny = hdr["nx"], hdr["ny"]
        frame = options.frame if options.frame >= 0 else nframes // 2
        frame = min(max(frame, 0), nframes - 1)
        img = EMData()
        img.read_image(options.movie, 0, False, Region(0, 0, frame, nx, ny, 1))
    else:
        nframes = nimg
        frame = options.frame if options.frame >= 0 else nframes // 2
        frame = min(max(frame, 0), nframes - 1)
        img = EMData(options.movie, frame)

    label = options.label or os.path.splitext(os.path.basename(options.movie))[0]
    os.makedirs(options.outdir, exist_ok=True)

    arr = img.numpy().flatten()
    n_pix = arr.size

    print(f"{options.movie}: frame {frame}/{nframes}, {n_pix} pixels, label='{label}'")

    # 1. General log-scale histogram of raw pixel values
    plt.figure(figsize=(6, 4))
    plt.hist(arr, 500, log=True)
    plt.xlabel("Pixel value")
    plt.ylabel("Number of pixels (log)")
    plt.title(f"Pixel value histogram — {label}\n(frame {frame}/{nframes})", fontsize=9)
    plt.tight_layout()
    hist_out = os.path.join(options.outdir, f"{label}_detectorHist_frame{frame}.jpg")
    plt.savefig(hist_out)
    plt.close()
    print(f"  wrote {hist_out}")

    # 2. Discrete electron-count histogram + fill fraction -- only meaningful
    #    (and only attempted) when the frame looks like raw integer counts.
    is_integer_like = np.allclose(arr, np.round(arr), atol=1e-3)
    vmin, vmax = float(arr.min()), float(arr.max())
    if is_integer_like and vmin >= 0 and vmax <= 64:
        counts = np.bincount(arr.astype(np.int64), minlength=int(vmax) + 1)
        values = np.arange(len(counts))
        fill_fraction = 1.0 - counts[0] / n_pix

        plt.figure(figsize=(10, 4))
        plt.bar(values, counts, log=True, width=1)
        plt.xlabel("Electron count")
        plt.ylabel("Number of pixels (log)")
        plt.title(f"Electron-count histogram — {label}  (frame {frame}/{nframes})  "
                  f"—  FillFraction = {fill_fraction:.4f}", fontsize=9)
        plt.tight_layout()
        counts_out = os.path.join(options.outdir, f"{label}_detectorCounts_frame{frame}.jpg")
        plt.savefig(counts_out)
        plt.close()
        print(f"  wrote {counts_out}")
        print(f"  FillFraction = {fill_fraction:.4f}")
    else:
        print(f"  Skipping electron-count histogram/fill-fraction (values not small "
              f"non-negative integers: min={vmin:.2f} max={vmax:.2f}) -- "
              f"frame doesn't look like raw counting-mode data.")


if __name__ == "__main__":
    main()
