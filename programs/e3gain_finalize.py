#!/usr/bin/env python
#
# Condition a gain estimate into a production gain reference.
#
# THE TWO-STEP GAIN PRESCRIPTION
#   1.  e3movie.py --movies <...> --est_gain raw_gain.hdf
#   2.  e3gain_finalize.py raw_gain.hdf gain.hdf
#
# WHY STEP 2 EXISTS
# -----------------
# ct_functions.py applies the gain with the EMAN processor math.fixgain.counting
# (gainmin=2, gainmax=2). From libEM/processor.cpp it computes
#       gainmax = mean + 2*sigma,   gainmin = max(0.1, mean - 2*sigma)
# and a pixel whose gain falls OUTSIDE that window is not divided at all -- it is
# REPLACED by the first good pixel along a quadrant-dependent diagonal. Zero is
# how you tell the pipeline "substitute this pixel rather than correcting it".
#
# A raw estimate carries no such flags. The detector's chip-boundary pixels collect
# extra charge (+3.7% measured on Fiducial1), so they come out BRIGHT, sit inside
# the window, and get divided by a slightly wrong number -- printing a dark cross
# and rectangle edges into every corrected frame.
#
# Two real failures, both this mechanism:
#   * 08152025 -- a "repair" filled the dead seam lines with ~1.03, un-masking them.
#     Seams went from -0.34% to -3.66%; restoring the mask fixed it.
#   * Fiducials -- NO failure of this kind. Their references were already masked
#     correctly. A hypothesis that the faint 4x2 "panels" came from unmasked
#     shoulders beside each boundary was tested by widening the mask to +-8 px and
#     is DISPROVED: it made the images visibly worse (see --width). The panels are
#     a detector NOISE property -- the boundary columns have Fano (var/mean) 2.96
#     against 3.93 elsewhere, so they carry the right brightness and ~13% less
#     standard deviation, which reads as a faint smooth line. A gain reference is
#     one divisor per pixel, and under division by g the mean goes as 1/g while
#     Fano also goes as 1/g, so brightness and noise cannot both be corrected.
#     NO gain file can remove them. Do not try again through this tool.
#
# WHY THE GEOMETRY IS NOT AUTO-DETECTED
# -------------------------------------
# It is a fixed property of the sensor: a 4096 detector is 4 chips across and 2
# down, so the internal boundaries are at nx/4, nx/2, 3nx/4 and ny/2, plus the
# outer edges. Detecting them from the image was tried and proved fragile -- on an
# already-masked reference the width estimate ran away, and on an unmasked estimate
# the weaker boundaries were missed entirely. Known geometry is simpler and cannot
# silently half-fail. Override with --cols/--rows for a different sensor.

import argparse
import os
import sys

import numpy as np

from EMAN2 import EMData


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="e3gain_finalize.py <raw_gain.hdf> <out_gain.hdf> [options]")
    p.add_argument("raw_gain", help="Gain estimate, e.g. from e3movie.py --est_gain.")
    p.add_argument("out_gain", help="Production gain reference to write.")
    p.add_argument("--width", type=int, default=0,
                   help="EXTRA half-width masked either side of each chip line. "
                        "DEFAULT 0 -- mask the lines themselves and nothing more. "
                        "DO NOT RAISE THIS without looking at the result by eye. "
                        "Widening to 8 was measured on Fiducial1 and Fiducial3 and "
                        "it MADE THE IMAGES WORSE: every seam statistic improved "
                        "(tilt-series seams 1.56%% -> 0.53%%) while the faint 4x2 "
                        "panels became blatant dashed bands, and the GainCheck "
                        "power spectrum gained a diagonal X. fixgain does not "
                        "correct a masked pixel, it SUBSTITUTES a neighbour's "
                        "value, so a wide masked band replaces independent noise "
                        "with duplicated noise. The mean level improves and the "
                        "texture is destroyed.")
    p.add_argument("--cols", type=str, default="",
                   help="Comma-separated chip column boundaries. Default: the "
                        "quarter points of the image width, plus both edges.")
    p.add_argument("--rows", type=str, default="",
                   help="Comma-separated chip row boundaries. Default: the "
                        "midpoint of the image height, plus both edges.")
    p.add_argument("--keepstats", action="store_true",
                   help="Report only; write nothing.")
    opt = p.parse_args()

    if not os.path.exists(opt.raw_gain):
        sys.exit(f"ERROR: no such file: {opt.raw_gain}")

    em = EMData(opt.raw_gain, 0)
    G = np.array(em.numpy(), dtype=np.float64, copy=True)
    ny, nx = G.shape
    gm, gs = float(em["mean"]), float(em["sigma"])

    pos = G[G > 0]
    if pos.size == 0:
        sys.exit("ERROR: gain image is entirely <= 0")
    print(f"{opt.raw_gain}: {nx}x{ny}  median {np.median(pos):.4f}  "
          f"sigma {gs:.4f}  zeros {(G <= 0).sum()}")

    if opt.cols:
        cols = [int(v) for v in opt.cols.split(",") if v.strip() != ""]
    else:
        cols = [0, nx//4 - 1, nx//4, nx//2 - 1, nx//2,
                3*nx//4 - 1, 3*nx//4, nx - 1]
    if opt.rows:
        rows = [int(v) for v in opt.rows.split(",") if v.strip() != ""]
    else:
        rows = [0, ny//2 - 1, ny//2, ny - 1]
    cols = sorted({c for c in cols if 0 <= c < nx})
    rows = sorted({r for r in rows if 0 <= r < ny})
    print(f"  chip columns {cols}")
    print(f"  chip rows    {rows}")

    before = ((G < max(0.1, gm - 2*gs)) | (G > gm + 2*gs)).sum()
    seam_before = [float(G[:, c].mean()) for c in cols if 0 < c < nx-1]
    print(f"  seam columns currently read "
          + " ".join(f"{v:.4f}" for v in seam_before)
          + f"   (median elsewhere {np.median(pos):.4f})")

    W = opt.width
    mask = (G <= 0)
    for c in cols:
        mask[:, max(0, c-W):c+W+1] = True
    for r in rows:
        mask[max(0, r-W):r+W+1, :] = True
    GW = G.copy()
    GW[mask] = 0.0
    print(f"  masking +-{W} px around each: {mask.sum()} px ({100*mask.mean():.2f}%)"
          f"; the input would have flagged {before} ({100*before/G.size:.2f}%)")

    if opt.keepstats:
        print("  --keepstats: nothing written")
        return

    out = EMData(em)
    out.numpy()[:] = GW.astype(np.float32)
    out.write_image(opt.out_gain, 0)
    del out, em

    chk = EMData(opt.out_gain, 0)
    a = np.array(chk.numpy(), dtype=np.float64, copy=True)
    hm, hs = float(chk["mean"]), float(chk["sigma"])
    del chk
    lo = max(0.1, hm - 2*hs)
    bad_lines = [c for c in cols if not bool((a[:, c] < lo).all())] + \
                [f"row{r}" for r in rows if not bool((a[r, :] < lo).all())]
    print(f"wrote {opt.out_gain}: mean {hm:.4f} sigma {hs:.4f}")
    print(f"  fixgain window [{lo:.4f}, {hm+2*hs:.4f}]")
    if bad_lines:
        sys.exit(f"ERROR: these lines are still INSIDE the window and would be "
                 f"divided rather than substituted: {bad_lines}. Do not use this "
                 f"file -- it would reproduce the dark-cross artefact.")
    print(f"  all {len(cols)} columns and {len(rows)} rows will be SUBSTITUTED: OK")


if __name__ == "__main__":
    main()
