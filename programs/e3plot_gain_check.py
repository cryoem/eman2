#!/usr/bin/env python
# Gain-correction efficacy diagnostic.
#
# Averages the SAME block of frames from the raw movie and from the gain-corrected
# stack, and shows both in real and Fourier space as a 2x2 panel. Fixed-pattern
# detector structure -- hot columns, the sensor grid, chip-boundary seams -- is
# multiplicative and does not average away, so it survives in the raw average and
# appears in its power spectrum as a sharp lattice of spikes and axis lines.
# Effective gain correction removes it: the real-space panel loses the pattern
# and the Fourier panel loses the spikes, leaving a smooth specimen spectrum.
#
# Companion to the RPS diagnostic in e3ctomo_import.py (_make_rps_diagnostic),
# which compares unshifted vs shift-corrected averages; this one isolates gain.
#
# Usage:
#   e3plot_gain_check.py <raw_movie.hdf> <gain_corrected.hdf> [--nframes 150]
#                        [--start N] [--outdir diagnostics] [--label NAME]

import argparse
import os
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from EMAN2 import EMData, EMUtil


def avg_block(path, s, e):
    """Mean of frames [s, e] read one at a time (blocks here are ~150 x 4096^2).

    EMData.numpy() returns a VIEW into the EMData's own buffer, not a copy. The
    obvious one-liner -- EMData(path,i).numpy().astype(float) -- silently yields
    zeros, because the temporary EMData is destroyed and its buffer freed before
    the copy is taken. Bind the EMData to a name so it outlives the copy.
    """
    acc = None
    for i in range(s, e + 1):
        em = EMData(path, i)
        a = np.array(em.numpy(), dtype=np.float64, copy=True)
        acc = a if acc is None else acc + a
        del em
    return acc / (e - s + 1)


def logps(a):
    """log10 power spectrum, windowed to keep edge wrap out of the spikes."""
    a = a - a.mean()
    n = min(a.shape)
    a = a[:n, :n] * np.outer(np.hanning(n), np.hanning(n))
    P = np.abs(np.fft.fftshift(np.fft.fft2(a))) ** 2
    return np.log10(P + 1e-12)


def show(ax, img, title, fourier=False):
    if fourier:
        lo, hi = np.percentile(img, (40, 99.9))
    else:
        lo, hi = np.percentile(img, (0.5, 99.5))
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, origin="lower")
    ax.set_title(title, fontsize=10.5)
    ax.set_xticks([]); ax.set_yticks([])


def main():
    p = argparse.ArgumentParser(
        usage="e3plot_gain_check.py <raw_movie.hdf> <gain_corrected.hdf> [options]",
        description="2x2 real/Fourier comparison of the same frames before and "
                    "after gain correction.")
    p.add_argument("raw", type=str, help="Raw (pre-gain) movie stack.")
    p.add_argument("corrected", type=str, help="Gain-corrected stack, same frame indexing.")
    p.add_argument("--nframes", type=int, default=150, help="Frames to average (default 150).")
    p.add_argument("--start", type=int, default=-1,
                   help="First frame; default centres the block in the movie.")
    p.add_argument("--outdir", type=str, default="diagnostics", help="Output directory.")
    p.add_argument("--label", type=str, default=None,
                   help="Override the label used in the filename and titles.")
    opt = p.parse_args()

    nr = EMUtil.get_image_count(opt.raw)
    nc = EMUtil.get_image_count(opt.corrected)
    if nr != nc:
        print(f"Warning: frame counts differ ({nr} vs {nc}); using the smaller.")
    n = min(nr, nc)

    start = opt.start if opt.start >= 0 else max(0, n // 2 - opt.nframes // 2)
    end = min(n - 1, start + opt.nframes - 1)

    label = opt.label or os.path.basename(opt.raw).replace(".hdf", "")
    apix = EMData(opt.raw, 0, True)["apix_x"]

    print(f"averaging frames {start}-{end} of {n}")
    raw = avg_block(opt.raw, start, end)
    cor = avg_block(opt.corrected, start, end)

    fig, ax = plt.subplots(2, 2, figsize=(10.4, 10.8))
    show(ax[0, 0], raw, "raw — real space")
    show(ax[0, 1], logps(raw), "raw — power spectrum (log)", fourier=True)
    show(ax[1, 0], cor, "gain corrected — real space")
    show(ax[1, 1], logps(cor), "gain corrected — power spectrum (log)", fourier=True)
    ax[0, 0].set_ylabel("before", fontsize=12, weight="bold")
    ax[1, 0].set_ylabel("after", fontsize=12, weight="bold")

    wrapped = "\n".join(textwrap.wrap(label, 62))
    fig.suptitle(f"Gain correction check\n{wrapped}\n"
                 f"mean of frames {start}–{end} ({end-start+1} frames), "
                 f"{raw.shape[1]}×{raw.shape[0]} at {apix:.3f} Å/pix",
                 fontsize=11, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    os.makedirs(opt.outdir, exist_ok=True)
    out = os.path.join(opt.outdir, f"GainCheck_{label}.jpg")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")

    # a couple of numbers to go with the pictures
    for tag, a in (("raw", raw), ("corrected", cor)):
        print(f"  {tag:9s} mean={a.mean():.4f}  sigma={a.std():.4f}  "
              f"sigma/mean={a.std()/abs(a.mean()):.4f}")


if __name__ == "__main__":
    main()
