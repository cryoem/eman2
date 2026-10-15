#!/usr/bin/env python
#
# Remove fiducial-bead shadows from a gain reference.
#
# A gain reference estimated from movies of a bead-bearing specimen carries the
# beads' own shadow: measured on GainFileFor_FiducialDec2025.hdf there are ~600
# compact dark blobs covering 0.34% of the detector, median depth -1.63% and worst
# -3.4%. Gain correction DIVIDES by this map, so every one of those spots makes
# the corrected image ~1.6% BRIGHT at a position fixed in detector coordinates --
# a ghost bead printed into every frame.
#
# The fix is to replace each bead pixel with a value drawn from the LOCAL mean and
# standard deviation. Three ways to fill were considered:
#   * mask them (gain -> 0)     : fixgain then substitutes, duplicating a
#                                 neighbour's noise. Measured to produce a visible
#                                 band and a diagonal X in the power spectrum.
#   * interpolate               : removes the real pixel-to-pixel gain variation,
#                                 which is NOT noise -- it is reproducible at
#                                 r=0.74 between independently measured references
#                                 and worth ~4% of genuine correction.
#   * draw from local mean+std  : preserves both the local level and the real
#                                 scatter. Used here.
#
# The chip-boundary mask is left exactly as found. Widening it was tried and made
# things materially worse.

import argparse
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, label, binary_dilation

from EMAN2 import EMData


def local_stats(G, good, sigma):
    """Local mean and std of G over `good` pixels only."""
    w = good.astype(np.float64)
    gw = gaussian_filter(np.where(good, G, 0.0), sigma)
    ww = gaussian_filter(w, sigma)
    ww = np.maximum(ww, 1e-6)
    mean = gw / ww
    g2w = gaussian_filter(np.where(good, G*G, 0.0), sigma)
    var = np.maximum(g2w/ww - mean*mean, 0.0)
    return mean, np.sqrt(var)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="e3gain_debead.py <gain_in.hdf> <gain_out.hdf> [options]")
    p.add_argument("gain_in")
    p.add_argument("gain_out")
    p.add_argument("--nsig", type=float, default=4.0,
                   help="Blob depth threshold in sigma of the bead-scale residual "
                        "(default 4).")
    p.add_argument("--minpix", type=int, default=6,
                   help="Smallest blob treated as a bead (default 6).")
    p.add_argument("--maxpix", type=int, default=800,
                   help="Largest blob treated as a bead (default 800); anything "
                        "bigger is structure, not a bead.")
    p.add_argument("--grow", type=int, default=2,
                   help="Dilate each blob by this many pixels to catch its "
                        "penumbra (default 2).")
    p.add_argument("--statsig", type=float, default=24.0,
                   help="Neighbourhood sigma for the local mean/std (default 24).")
    p.add_argument("--seed", type=int, default=0, help="RNG seed (default 0).")
    p.add_argument("--dryrun", action="store_true", help="Report only.")
    opt = p.parse_args()

    if not os.path.exists(opt.gain_in):
        sys.exit(f"ERROR: no such file: {opt.gain_in}")

    em = EMData(opt.gain_in, 0)
    G = np.array(em.numpy(), dtype=np.float64, copy=True)
    live = G > 0.3                                  # dead/masked pixels stay as they are
    med = np.median(G[live])
    Gn = np.where(live, G/med, 1.0)
    print(f"{opt.gain_in}: {G.shape}, median {med:.4f}, "
          f"{int((~live).sum())} px already masked", flush=True)

    band = gaussian_filter(Gn, 2.5) - gaussian_filter(Gn, 25)
    sd = band[200:-200, 200:-200].std()
    dark = band < -opt.nsig*sd
    lab, n = label(dark)
    sizes = np.bincount(lab.ravel())
    keep = np.zeros(sizes.shape[0], dtype=bool)
    keep[opt.minpix <= sizes] = True
    keep[sizes > opt.maxpix] = False
    keep[0] = False
    beads = keep[lab]
    if opt.grow > 0:
        k = 2*opt.grow + 1
        beads = binary_dilation(beads, np.ones((k, k), bool))
    beads &= live
    nb = int(keep.sum())
    print(f"  bead-scale scatter {100*sd:.3f}%; {nb} blobs -> {beads.sum()} px "
          f"({100*beads.mean():.3f}% of the detector)", flush=True)
    if nb:
        dep = band[beads].min()
        print(f"  deepest shadow {100*dep:+.2f}%  (these print as BRIGHT ghosts "
              f"of the same size after division)", flush=True)
    if opt.dryrun:
        print("  --dryrun: nothing written")
        return
    if nb == 0:
        sys.exit("ERROR: no bead shadows found; lower --nsig or check the input.")

    good = live & ~beads
    mean, std = local_stats(G, good, opt.statsig)
    rng = np.random.default_rng(opt.seed)
    draw = mean + rng.standard_normal(G.shape)*std
    GN = np.where(beads, draw, G)
    GN = np.where(live, GN, G)                      # never resurrect a masked pixel

    b = beads
    print(f"  replaced bead pixels: mean {G[b].mean():.4f} -> {GN[b].mean():.4f} "
          f"(local mean {mean[b].mean():.4f})", flush=True)
    print(f"  scatter there:        {G[b].std():.4f} -> {GN[b].std():.4f} "
          f"(local std {std[b].mean():.4f})", flush=True)

    out = EMData(em)
    out.numpy()[:] = GN.astype(np.float32)
    out.write_image(opt.gain_out, 0)
    del out, em
    chk = EMData(opt.gain_out, 0)
    a = np.array(chk.numpy(), dtype=np.float64, copy=True)
    hm, hs = float(chk["mean"]), float(chk["sigma"])
    del chk
    same_mask = int(((a <= 0.3) != (~live)).sum())
    print(f"wrote {opt.gain_out}: mean {hm:.4f} sigma {hs:.4f}")
    print(f"  mask unchanged: {same_mask == 0}")
    if same_mask:
        sys.exit("ERROR: the dead-pixel mask changed; do not use this file.")


if __name__ == "__main__":
    main()
