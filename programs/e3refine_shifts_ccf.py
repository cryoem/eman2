#!/usr/bin/env python
#
# Refine tilt-series shifts by cross-correlating CONSECUTIVE ALIGNED tilts.
#
# Why this exists
# ---------------
# e2tomogram fits 20-50 landmark patches. That leaves residual error the fit
# cannot see, because satisfying a few dozen patch positions does not guarantee
# the whole frame is registered. Measured on ct1 (08152025): the residual shift
# between consecutive aligned tilts is 0.00 px through the middle of the series
# but 1.4-1.7 px at both ends -- exactly where the stack looks jittery by eye.
# At bin 8 that is ~11 detector pixels.
#
# This does what the eye does: compare each tilt with its neighbours and nudge.
#
# The part that matters
# ---------------------
# A pairwise cross-correlation is LEAST reliable exactly where the residual is
# largest -- on ct1 the CCF peak SNR is 7.0 at the ends against 15.5 in the
# middle, because high-tilt frames are dimmer and lower contrast. Trusting each
# pair individually would inject noise precisely where we are trying to fix
# something. So the corrections are solved globally:
#
#   minimise  sum_i w_i * ||(c_i - c_{i-1}) - r_i||^2        (match the data)
#           + lam * sum_i ||c_{i+1} - 2 c_i + c_{i-1}||^2    (stay smooth)
#           + gauge * sum_i ||c_i||^2                        (pin the mean)
#
# with w_i from the CCF peak SNR. A low-confidence tilt is then carried by its
# well-measured neighbours instead of being believed on its own -- which is the
# continuity an eye uses and a pairwise method lacks.
#
# Smoothing does NOT span the negative/positive junction: the shift genuinely
# jumps there (the stage returns and the specimen settles elsewhere) and
# smoothing across it erodes a real offset.
#
# The tool re-measures after applying and reports before/after. A sign error
# would show up immediately as residuals getting worse.

import argparse
import json
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

from EMAN2 import EMData, EMUtil


def prep(path, i, hp, crop):
    e = EMData(path, i)
    a = np.array(e.numpy(), dtype=np.float64, copy=True)
    del e
    if crop > 0:
        a = a[crop:-crop, crop:-crop]
    a = a - gaussian_filter(a, hp)
    s = a.std()
    if s > 0:
        a = a / s
    return a * np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))


def ccf(a, b, maxs):
    F = np.fft.fft2(a) * np.conj(np.fft.fft2(b))
    c = np.fft.fftshift(np.real(np.fft.ifft2(F)))
    cy, cx = np.array(c.shape) // 2
    w = c[cy-maxs:cy+maxs+1, cx-maxs:cx+maxs+1]
    j = np.unravel_index(np.argmax(w), w.shape)
    pk = w[j]
    # parabolic refinement to sub-pixel
    def sub(v0, v1, v2):
        d = v0 - 2*v1 + v2
        return 0.0 if abs(d) < 1e-12 else 0.5*(v0 - v2)/d
    oy = sub(w[max(j[0]-1,0), j[1]], pk, w[min(j[0]+1, w.shape[0]-1), j[1]])
    ox = sub(w[j[0], max(j[1]-1,0)], pk, w[j[0], min(j[1]+1, w.shape[1]-1)])
    bg = np.std(w)
    return (j[1]-maxs+ox), (j[0]-maxs+oy), (pk/bg if bg > 0 else 0.0)


def measure(path, hp, crop, maxs, verbose):
    n = EMUtil.get_image_count(path)
    prev = prep(path, 0, hp, crop)
    r = np.zeros((n, 2))
    w = np.zeros(n)
    for i in range(1, n):
        cur = prep(path, i, hp, crop)
        dx, dy, snr = ccf(prev, cur, maxs)
        r[i] = (dx, dy)
        w[i] = snr
        prev = cur
        if verbose and i % 20 == 0:
            print(f"    tilt {i}: residual ({dx:+.2f},{dy:+.2f}) snr {snr:.1f}", flush=True)
    return r, w, n


def solve(r, w, n, lam, gauge, junction):
    """Least squares for per-tilt corrections, one axis at a time."""
    out = np.zeros((n, 2))
    for ax in (0, 1):
        rows, cols, vals, rhs = [], [], [], []
        k = 0
        for i in range(1, n):                       # data: c_i - c_{i-1} = r_i
            ww = np.sqrt(max(w[i], 0.0))
            if ww <= 0:
                continue
            rows += [k, k]; cols += [i, i-1]; vals += [ww, -ww]
            rhs.append(ww * r[i, ax]); k += 1
        sl = np.sqrt(lam)
        for i in range(1, n-1):                     # smoothness, not across the junction
            if junction >= 0 and (i-1 <= junction <= i+1):
                continue
            rows += [k, k, k]; cols += [i-1, i, i+1]; vals += [sl, -2*sl, sl]
            rhs.append(0.0); k += 1
        sg = np.sqrt(gauge)
        for i in range(n):                          # gauge
            rows.append(k); cols.append(i); vals.append(sg); rhs.append(0.0); k += 1
        A = np.zeros((k, n))
        A[rows, cols] = vals
        # accumulate duplicates correctly
        A[:] = 0.0
        for rr, cc, vv in zip(rows, cols, vals):
            A[rr, cc] += vv
        sol, *_ = np.linalg.lstsq(A, np.array(rhs), rcond=None)
        out[:, ax] = sol
    return out


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="e3refine_shifts_ccf.py <recon_tomo_final.json> <measured_stack.hdf> [options]")
    p.add_argument("json")
    p.add_argument("stack", help="the ALIGNED stack, e.g. diagnostics/FidStack_measured_<bn>_SA.hdf")
    p.add_argument("--out", default="", help="output JSON (default <in>_ccf.json)")
    p.add_argument("--lam", type=float, default=10.0,
                   help="smoothness weight (default 10). Higher trusts the "
                        "neighbours more and the individual CCF less.")
    p.add_argument("--gauge", type=float, default=1e-3, help="weight pinning the mean correction to zero")
    p.add_argument("--maxshift", type=int, default=20, help="largest residual searched, px (default 20)")
    p.add_argument("--hp", type=float, default=20.0, help="high-pass sigma before correlating")
    p.add_argument("--crop", type=int, default=40, help="border to drop before correlating")
    p.add_argument("--junction", type=int, default=-2,
                   help="last tilt index of the negative acquisition; smoothing "
                        "will not span it. -2 = auto from info/sa_chunks_*, "
                        "-1 = disable")
    p.add_argument("--dryrun", action="store_true", help="measure and report, write nothing")
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args()

    for f in (a.json, a.stack):
        if not os.path.exists(f):
            sys.exit(f"ERROR: no such file: {f}")
    rec = json.load(open(a.json))
    tp = np.array(rec["tlt_params"], dtype=float)
    em = EMData(a.stack, 0, True); sz = em["nx"]; del em
    scale = 1024.0 / sz                     # tlt_params are in 1024-px units
    print(f"{a.stack}: {sz} px, tlt_params scale x{scale:g}")

    junction = a.junction
    if junction == -2:
        junction = -1
        bn = os.path.basename(a.json).replace("_recon_tomo_final.json", "")
        cj = os.path.join(os.path.dirname(a.json), f"sa_chunks_{bn}.json")
        if os.path.exists(cj):
            d = json.load(open(cj))
            ch = d.get("chunks") or []
            ns = d.get("n_sel_neg")
            if ch and ns is not None:
                for i, (s, e) in enumerate(ch):
                    if e >= ns - 1:
                        junction = i - 1
                        break
                print(f"  junction after tilt {junction} (from {os.path.basename(cj)})")
        if junction < 0:
            print("  WARNING: could not determine the junction; smoothing will span it")

    print("  measuring residuals between consecutive aligned tilts ...", flush=True)
    r, w, n = measure(a.stack, a.hp, a.crop, a.maxshift, a.verbose)
    d0 = np.hypot(r[1:, 0], r[1:, 1])
    print(f"  BEFORE: median |residual| {np.median(d0):.2f} px   mean {d0.mean():.2f}   max {d0.max():.2f}")
    print(f"          CCF snr median {np.median(w[1:]):.1f}  (weights come from this)")

    c = solve(r, w, n, a.lam, a.gauge, junction)
    pred = np.hypot(*(np.diff(c, axis=0) - r[1:]).T)
    print(f"  AFTER (predicted): median |residual| {np.median(pred):.2f} px   "
          f"mean {pred.mean():.2f}   max {pred.max():.2f}")
    print(f"  correction applied: median {np.median(np.hypot(c[:,0],c[:,1])):.2f} px, "
          f"max {np.hypot(c[:,0],c[:,1]).max():.2f} px (stack units)")
    if np.median(pred) > np.median(d0):
        print("  WARNING: predicted residual is WORSE than before -- check the sign convention")

    if a.dryrun:
        print("  --dryrun: nothing written")
        return
    tp2 = tp.copy()
    tp2[:, 0] -= c[:, 0] * scale
    tp2[:, 1] -= c[:, 1] * scale
    rec2 = dict(rec); rec2["tlt_params"] = tp2.tolist()
    out = a.out or a.json.replace(".json", "_ccf.json")
    json.dump(rec2, open(out, "w"), indent=1)
    j0 = np.sqrt(np.diff(tp[:,0],2)**2 + np.diff(tp[:,1],2)**2).std()
    j1 = np.sqrt(np.diff(tp2[:,0],2)**2 + np.diff(tp2[:,1],2)**2).std()
    print(f"wrote {out}")
    print(f"  jitter (rms 2nd difference of the shifts): {j0:.2f} -> {j1:.2f}")


if __name__ == "__main__":
    main()
