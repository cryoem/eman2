#!/usr/bin/env python
# Fiducial placement error by reprojection.
#
# Reprojects a tomogram at the Euler angles of each tilt and compares the result
# with the measured tilt, and -- independently -- measures where each fiducial
# bead actually sits in every tilt relative to where the alignment model says it
# should be. The second measurement is the quantitative one: it needs no
# tomogram, only the tilt images and the refined tilt parameters, so it carries
# no reconstruction or interpolation error of its own.
#
# Conventions follow e2tomogram.py exactly:
#   * get_xf_pos(tpm, pk) (e2tomogram.py:1502) maps a 3D point to its 2D
#     position in a tilt, relative to the image centre.
#   * reconstruct() (e2tomogram.py:2117) aligns each tilt before insertion by
#     translating -int(tx), -int(ty) and rotating -ztilt, then inserts it with
#     Transform({"type":"xyz","ytilt":..,"xtilt":..}). A reprojection at those
#     same two angles is therefore comparable to the ALIGNED tilt, not the raw
#     one -- comparing against the raw tilt would fold the alignment shifts back
#     in and inflate every residual.
#
# Usage:
#   e3fid_reproject.py --tomo tomograms/<name>.hdf --tilts tiltseries/<name>.hdf \
#                      --ttparams <name>_ttparams.npy --pks <name>_pks.npy \
#                      [--bin 2] [--box 48] [--outdir diagnostics] [--label NAME]
#
#   --skip_reproject   bead residuals only (seconds; no tomogram read)
#   --skip_beads       reprojection/CCC only

import argparse
import hashlib
import json
import os
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from EMAN2 import EMData, EMUtil, Region, Transform

try:
    from scipy.ndimage import maximum_filter
    _HAS_SCIPY = True
except ImportError:                                  # nearest-peak needs scipy;
    _HAS_SCIPY = False                               # fall back to strongest-peak


# ---------------------------------------------------------------- geometry ---

def get_xf_pos(tpm, pk):
    """2D position of 3D point pk in the tilt with parameters tpm.

    Verbatim from e2tomogram.py:1502 so the prediction matches the model that
    was actually refined. tpm is [tx, ty, ztilt, ytilt, xtilt].
    """
    xf0 = Transform({"type": "xyz", "xtilt": tpm[4], "ytilt": tpm[3],
                     "ztilt": tpm[2], "tx": tpm[0], "ty": tpm[1]})
    p1 = xf0.transform([float(pk[0]), float(pk[1]), float(pk[2])])
    return np.array([p1[0], p1[1]])


def align_tilt(img, tpm, pad=None):
    """Reproduce the aligned tilt that went into the reconstruction."""
    m = img.copy()
    m.process_inplace("filter.ramp")
    m.process_inplace("normalize")
    nx, ny = m["nx"], m["ny"]
    if pad is None:
        pad = nx
    m.process_inplace("mask.decayedge2d", {"width": int(pad // 40)})
    p2 = m.get_clip(Region(nx // 2 - pad // 2, ny // 2 - pad // 2, pad, pad), fill=0)
    tx, ty = get_xf_pos(tpm, [0, 0, 0])
    p2.translate(-int(tx), -int(ty), 0)
    p2.rotate(-float(tpm[2]), 0, 0)
    return p2


def to_np(em):
    """EMData.numpy() is a VIEW into the object's buffer -- copy explicitly."""
    return np.array(em.numpy(), dtype=np.float64, copy=True)


def xcorr_shift(a, b, maxsh=10):
    """Offset of b relative to a, by cross-correlation with sub-pixel refinement.

    Returns (dx, dy) or None if the peak sits on the edge of the search range,
    which means the true shift exceeds maxsh and the fitted value would be
    meaningless.
    """
    a = a - a.mean()
    b = b - b.mean()
    if a.std() < 1e-9 or b.std() < 1e-9:
        return None
    F = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b)))))
    c = F.shape[0] // 2
    w = F[c - maxsh:c + maxsh + 1, c - maxsh:c + maxsh + 1]
    if w.shape[0] < 3 or w.shape[1] < 3:
        return None
    iy, ix = np.unravel_index(np.argmax(w), w.shape)
    if iy == 0 or ix == 0 or iy == w.shape[0] - 1 or ix == w.shape[1] - 1:
        return None
    dx = 0.5 * (w[iy, ix - 1] - w[iy, ix + 1]) / \
        (w[iy, ix - 1] - 2 * w[iy, ix] + w[iy, ix + 1] + 1e-12)
    dy = 0.5 * (w[iy - 1, ix] - w[iy + 1, ix]) / \
        (w[iy - 1, ix] - 2 * w[iy, ix] + w[iy + 1, ix] + 1e-12)
    return (ix + dx - maxsh, iy + dy - maxsh)


# ------------------------------------------------------------ bead finding ---

def refine_peak(box, polarity, rad=3, minsnr=3.0):
    """Sub-pixel bead centre within box, as an offset from the box centre.

    polarity is +1 if beads are bright, -1 if dark. Returns (dx, dy, strength)
    or None if no acceptable peak is present.

    Takes the candidate peak NEAREST the box centre, not the strongest. The
    search box can contain several beads, and a neighbouring bead is often the
    brighter one; picking by strength then reports the distance to the wrong
    bead, which manifests as a residual comparable to the box half-width. That
    is a mis-association, not a placement error, and it silently inflates the
    tail of the distribution.
    """
    b = (box - box.mean()) * polarity
    s = b.std()
    if s < 1e-9:
        return None
    n = b.shape[0]
    c = n // 2

    if _HAS_SCIPY:
        mx = maximum_filter(b, size=2 * rad + 1)
        cand = np.argwhere((b == mx) & (b > minsnr * s))
        if not len(cand):
            return None
        # nearest to the predicted (box-centre) position
        k = cand[np.argmin(np.hypot(cand[:, 0] - c, cand[:, 1] - c))]
        iy, ix = int(k[0]), int(k[1])
    else:
        iy, ix = np.unravel_index(np.argmax(b), b.shape)
        if b[iy, ix] <= minsnr * s:
            return None
    # reject peaks pinned to the border: the true centre is outside the box
    if ix < rad or iy < rad or ix >= n - rad or iy >= n - rad:
        return None
    # centre of mass of a small window about that peak (sub-pixel, noise robust)
    w = b[iy - rad:iy + rad + 1, ix - rad:ix + rad + 1].copy()
    w -= w.min()
    tot = w.sum()
    if tot <= 0:
        return None
    yy, xx = np.mgrid[-rad:rad + 1, -rad:rad + 1]
    cy = (w * yy).sum() / tot
    cx = (w * xx).sum() / tot
    return (ix + cx - c, iy + cy - c, b[iy, ix] / s)


def measure_beads(tilts, ttparams, pks, box, minsnr, verbose):
    """Residual (measured - predicted) bead position for every bead and tilt.

    Returns res[ntilt, nbead, 2] with NaN where the bead was not found.
    """
    ntilt, nbead = len(ttparams), len(pks)
    im0 = EMData(tilts, 0, True)
    nx, ny = im0["nx"], im0["ny"]
    res = np.full((ntilt, nbead, 2), np.nan)
    snr = np.full((ntilt, nbead), np.nan)

    # polarity: look at the zero-tilt image, at the predicted bead centres
    izero = int(np.argmin(np.abs(ttparams[:, 3])))
    ref = to_np(EMData(tilts, izero))
    dev = []
    for p in pks:
        c = get_xf_pos(ttparams[izero], p) + np.array([nx // 2, ny // 2])
        x, y = int(round(c[0])), int(round(c[1]))
        if box <= x < nx - box and box <= y < ny - box:
            sub = ref[y - box // 2:y + box // 2, x - box // 2:x + box // 2]
            dev.append(sub[box // 2 - 2:box // 2 + 3,
                           box // 2 - 2:box // 2 + 3].mean() - sub.mean())
    polarity = 1.0 if (len(dev) and np.mean(dev) > 0) else -1.0
    print(f"  bead polarity: {'bright' if polarity > 0 else 'dark'} "
          f"(mean central deviation {np.mean(dev) if dev else float('nan'):+.4f})")

    for i in range(ntilt):
        a = to_np(EMData(tilts, i))
        for j, p in enumerate(pks):
            c = get_xf_pos(ttparams[i], p) + np.array([nx // 2, ny // 2])
            x, y = int(round(c[0])), int(round(c[1]))
            if not (box // 2 <= x < nx - box // 2 and box // 2 <= y < ny - box // 2):
                continue
            sub = a[y - box // 2:y + box // 2, x - box // 2:x + box // 2]
            r = refine_peak(sub, polarity)
            if r is None or r[2] < minsnr:
                continue
            # the integer part of the predicted position was rounded away above,
            # so add the rounding back to keep the residual sub-pixel exact
            res[i, j, 0] = r[0] + (x - c[0])
            res[i, j, 1] = r[1] + (y - c[1])
            snr[i, j] = r[2]
        if verbose and i % 10 == 0:
            n = np.isfinite(res[i, :, 0]).sum()
            print(f"    tilt {i:3d}  ytilt {ttparams[i,3]:+7.2f}  {n}/{nbead} beads")
    return res, snr


# ------------------------------------------------------------ reprojection ---

def reproject(tomo_path, tilts, ttparams, pks, binning, lowpass, highpass,
              cccrad, obox, verbose, stackdir=None, label=None, stackbase=None,
              pscale=1.0, ttparams_path=""):
    """Per-tilt CCC between the reprojected tomogram and the aligned tilt.

    Two corrections are essential for the number to mean anything:

    * Matched band. The tomogram was built from ramp-filtered slices and then
      low-pass filtered (40 A here), and it carries slow density gradients the
      tilt does not. Correlating it against an unfiltered tilt compares a
      band-limited image with a broadband one, and the mismatched frequencies
      enter the CCC as pure disagreement. Both are band-passed identically.

    * Field of view. At tilt theta the measured image sees specimen out to
      x/cos(theta) in tomogram coordinates, but the volume only extends to its
      own half width, so the periphery of the tilt has no counterpart in the
      reprojection at all. The CCC is therefore reported over a central disc
      (--cccrad) that stays inside the volume at every tilt; the full-frame
      value is returned alongside it for reference.

    Also returns the per-bead OFFSET between the reprojection and the measured
    tilt, by local cross-correlation in a box about each predicted bead. This
    is the quantity that actually answers "is the fiducial in the right place":
    it asks whether the bead in the reconstructed volume lands where the bead
    in the data lands. The CCC does not answer that -- a fiducial is only a few
    pixels across, so a 2-3 pixel offset destroys its correlation while still
    looking well placed by eye, and the CCC then reads as catastrophic when the
    geometry is nearly right.

    Note this differs from the residual measured by measure_beads(): that one
    compares the alignment MODEL against the data and needs no tomogram, while
    this compares the RECONSTRUCTION against the data. They can disagree, and
    where they do, the difference is attributable to the reconstruction rather
    than to the alignment.
    """
    print(f"  reading {tomo_path}")
    tomo = EMData(tomo_path)
    if binning > 1:
        tomo.process_inplace("math.meanshrink", {"n": binning})
        print(f"    shrunk by {binning} -> "
              f"{tomo['nx']}x{tomo['ny']}x{tomo['nz']}")

    apix_out = tomo["apix_x"]

    def bandpass(em):
        if lowpass > 0:
            em.process_inplace("filter.lowpass.gauss", {"cutoff_freq": 1.0 / lowpass})
        if highpass > 0:
            em.process_inplace("filter.highpass.gauss", {"cutoff_freq": 1.0 / highpass})
        em.process_inplace("normalize")
        return em

    ntilt = len(ttparams)
    ccc = np.full(ntilt, np.nan)          # central disc
    ccc_full = np.full(ntilt, np.nan)     # whole frame
    off = np.full((ntilt, len(pks), 2), np.nan)   # reprojection -> measured
    panels = {}

    # Optional stacks for side-by-side inspection in e2display. Identical
    # indexing, so tilt i is image i in each: flipping between them at a fixed
    # frame shows the reprojection against its own measured tilt.
    #
    # The 'measured' stack is the aligned tilt series. It depends only on the
    # tilts, the tilt parameters, the binning and the band-pass -- never on the
    # tomogram -- so analysing several reconstructions built from ONE alignment
    # (SA, stage 2, ...) reproduces the same images.
    #
    # It used to be named for the dataset alone, without the run label, so those
    # runs could share one file. But runs of DIFFERENT alignments share a dataset
    # too, and the second silently overwrote the first while 'reproj' and 'diff'
    # kept their labels: Fiducial2/diagnostics/ ended up holding a fixytilt
    # reproj beside a no-fixytilt measured, 29 minutes apart, which compare as a
    # spurious mismatch. It now carries the same label as the other two, so a
    # run folder always holds three files describing one alignment.
    stacks = {}
    if stackdir:
        os.makedirs(stackdir, exist_ok=True)
        for tag in ("reproj", "diff"):
            pth = os.path.join(stackdir, f"FidStack_{tag}_{label}.hdf")
            if os.path.isfile(pth):
                os.remove(pth)
            stacks[tag] = pth
        mpth = os.path.join(stackdir, f"FidStack_measured_{label or stackbase}.hdf")
        # Reuse is only valid if the stack was written with the SAME parameters.
        # The aligned tilts depend on the tilt parameters (including any pixel
        # rescaling), the binning and the band-pass -- not on the tomogram. An
        # earlier version compared only the frame count, so a stack written from
        # wrongly-scaled parameters was silently reused by every later run: the
        # numbers stayed correct (meas is recomputed in memory) while the file on
        # disk kept translations 4x too large, showing enormous apparent motion.
        # The signature must capture everything the aligned tilts depend on.
        # Parameters alone are not enough: rebuilding a dataset (e.g. after the
        # per-half subaveraging fix) leaves scale/bin/filters/count unchanged
        # while the tilt series and the refined parameters are both different,
        # so a content digest of those two inputs is included as well.
        def _digest(path):
            try:
                h = hashlib.md5()
                h.update(str(os.path.getsize(path)).encode())
                h.update(str(int(os.path.getmtime(path))).encode())
                return h.hexdigest()[:12]
            except Exception:
                return "?"
        sig = (f"scale={pscale:.6f} bin={binning} lp={lowpass} hp={highpass} "
               f"n={ntilt} tilts={_digest(tilts)} tpm={_digest(ttparams_path)}")
        reuse = False
        if os.path.isfile(mpth):
            try:
                if (EMUtil.get_image_count(mpth) == ntilt and
                        EMData(mpth, 0, True).get_attr_default("fid_sig", "") == sig):
                    reuse = True
            except Exception:
                reuse = False
        if reuse:
            print(f"  measured stack already present and matches, reusing {mpth}")
        else:
            if os.path.isfile(mpth):
                print(f"  measured stack present but parameters differ, regenerating")
                os.remove(mpth)
            stacks["measured"] = mpth
            stack_sig = sig

    for i in range(ntilt):
        xf = Transform({"type": "xyz", "ytilt": float(ttparams[i, 3]),
                        "xtilt": float(ttparams[i, 4])})
        pj = bandpass(tomo.project("standard", xf))

        meas = align_tilt(EMData(tilts, i), ttparams[i])
        if binning > 1:
            meas.process_inplace("math.meanshrink", {"n": binning})
        meas = bandpass(meas)

        a, b = to_np(pj), to_np(meas)
        ny, nx = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
        a, b = a[:ny, :nx], b[:ny, :nx]
        yy, xx = np.mgrid[0:ny, 0:nx]
        rr = np.hypot(yy - ny // 2, xx - nx // 2)
        # ignore the zero border the high-tilt mask leaves behind
        valid = (np.abs(a) > 1e-9) & (np.abs(b) > 1e-9)

        def _c(m):
            if m.sum() < 100:
                return np.nan
            av, bv = a[m] - a[m].mean(), b[m] - b[m].mean()
            return float((av * bv).sum() / np.sqrt((av**2).sum() * (bv**2).sum()))

        ccc[i] = _c(valid & (rr < cccrad * min(ny, nx) // 2))
        ccc_full[i] = _c(valid)
        panels[i] = (b, a, b - a)

        # per-bead offset, reprojection -> measured, in the binned frame
        H = obox // 2
        for j, pk in enumerate(pks):
            q = get_xf_pos(ttparams[i], pk) / float(binning) + \
                np.array([nx // 2, ny // 2])
            x, y = int(round(q[0])), int(round(q[1]))
            if not (H <= x < nx - H and H <= y < ny - H):
                continue
            s = xcorr_shift(a[y - H:y + H, x - H:x + H],
                            b[y - H:y + H, x - H:x + H])
            if s is not None:
                off[i, j] = s

        if stacks:
            dif = meas.copy()
            dif.sub(pj)
            for tag, em in (("reproj", pj), ("measured", meas), ("diff", dif)):
                if tag not in stacks:          # measured reused from a prior run
                    continue
                em["apix_x"] = em["apix_y"] = em["apix_z"] = apix_out
                em["ytilt"] = float(ttparams[i, 3])
                em["xtilt"] = float(ttparams[i, 4])
                em["tilt_idx"] = i
                if tag != "measured":   # CCC belongs to a reconstruction, and
                    em["ccc_central"] = float(ccc[i])   # measured is shared
                else:
                    em["fid_sig"] = stack_sig
                em.write_image(stacks[tag], i)
        if verbose and i % 10 == 0:
            print(f"    tilt {i:3d}  ytilt {ttparams[i,3]:+7.2f}  "
                  f"CCC {ccc[i]:.4f} (full frame {ccc_full[i]:.4f})")
    # Verify what actually landed on disk. Images are written one at a time, so
    # an interrupted run leaves a short-but-valid HDF rather than a corrupt one:
    # it opens, displays and looks fine, and only the frame count reveals it.
    for tag in sorted(stacks):
        try:
            n = EMUtil.get_image_count(stacks[tag])
        except Exception as e:
            print(f"  WARNING: cannot read back {stacks[tag]}: {e}")
            continue
        if n == ntilt:
            print(f"  wrote {stacks[tag]}  ({n} images)")
        else:
            print(f"  WARNING: {stacks[tag]} has {n} images, expected {ntilt} "
                  f"-- truncated, rerun to regenerate")
    return ccc, ccc_full, off, panels


# ------------------------------------------------------------------ plots ----

def plot_montage(panels, ttparams, ccc, label, out, ntile=5):
    idx = np.linspace(0, len(ttparams) - 1, ntile).astype(int)
    fig, ax = plt.subplots(3, ntile, figsize=(3.0 * ntile, 9.4))
    for c, i in enumerate(idx):
        meas, rep, dif = panels[i]
        for r, (img, nm) in enumerate(((meas, "measured"), (rep, "reprojected"),
                                       (dif, "difference"))):
            lo, hi = np.percentile(img, (1, 99))
            ax[r, c].imshow(img, cmap="gray", vmin=lo, vmax=hi, origin="lower")
            ax[r, c].set_xticks([]); ax[r, c].set_yticks([])
            if c == 0:
                ax[r, c].set_ylabel(nm, fontsize=11, weight="bold")
        ax[0, c].set_title(f"{ttparams[i,3]:+.1f}°\nCCC {ccc[i]:.3f}", fontsize=10)
    fig.suptitle(f"Tomogram reprojected at the tilt Euler angles\n{label}",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_offsets(off, ttparams, pks, apix_bin, label, out, zeroid=None):
    """Reprojection-vs-measured bead offsets: the placement error itself."""
    ang = ttparams[:, 3]
    mag = np.hypot(off[:, :, 0], off[:, :, 1])
    med = np.nanmedian(mag, axis=1)

    fig = plt.figure(figsize=(13.6, 8.6))
    gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.30)

    # (a) |offset| vs tilt angle
    a0 = fig.add_subplot(gs[0, 0])
    for j in range(off.shape[1]):
        a0.plot(ang, mag[:, j] * apix_bin, color="#8899aa", lw=0.5, alpha=0.45)
    a0.plot(ang, med * apix_bin, color="#c0392b", lw=2.2, label="median over beads")
    a0.set_xlabel("tilt angle (deg)"); a0.set_ylabel("|offset| (Å)")
    a0.set_title("(a) reprojection vs measured", fontsize=11, weight="bold")
    a0.grid(alpha=0.25); a0.legend(fontsize=8)

    # (b) signed components -- separates drift along the tilt axis from across it
    a1 = fig.add_subplot(gs[0, 1])
    a1.plot(ang, np.nanmedian(off[:, :, 0], axis=1) * apix_bin, color="#2f4a6d",
            lw=1.8, label="dx  (⊥ tilt axis)")
    a1.plot(ang, np.nanmedian(off[:, :, 1], axis=1) * apix_bin, color="#1a8f4c",
            lw=1.8, label="dy  (∥ tilt axis)")
    a1.axhline(0, color="#666", lw=0.8)
    if zeroid is not None:
        a1.axvline(ang[zeroid], color="#c0392b", ls=":", lw=1.3)
    a1.set_xlabel("tilt angle (deg)"); a1.set_ylabel("median offset (Å)")
    a1.set_title("(b) signed components", fontsize=11, weight="bold")
    a1.grid(alpha=0.25); a1.legend(fontsize=8)

    # (c) negative vs positive half
    a2 = fig.add_subplot(gs[0, 2])
    zi = zeroid if zeroid is not None else int(np.argmin(np.abs(ang)))
    negm = mag[:zi + 1][np.isfinite(mag[:zi + 1])] * apix_bin
    posm = mag[zi + 1:][np.isfinite(mag[zi + 1:])] * apix_bin
    bins = np.linspace(0, np.nanpercentile(np.r_[negm, posm], 97), 34)
    a2.hist(negm, bins=bins, alpha=0.62, color="#2f4a6d",
            label=f"negative half (med {np.median(negm):.0f} Å)")
    a2.hist(posm, bins=bins, alpha=0.62, color="#c0392b",
            label=f"positive half (med {np.median(posm):.0f} Å)")
    a2.set_xlabel("|offset| (Å)"); a2.set_ylabel("count")
    a2.set_title("(c) by half-sweep", fontsize=11, weight="bold")
    a2.legend(fontsize=8); a2.grid(alpha=0.25)

    # (d) offset vs distance from the tilt axis
    a3 = fig.add_subplot(gs[1, 0])
    pb = np.nanmedian(mag, axis=0)
    rad = np.abs(pks[:, 0])
    a3.scatter(rad * apix_bin / 1e4, pb * apix_bin, s=40, c="#2f4a6d", zorder=3)
    g = np.isfinite(pb)
    if g.sum() > 2:
        k = np.polyfit(rad[g], pb[g], 1)
        xs = np.linspace(rad.min(), rad.max(), 40)
        a3.plot(xs * apix_bin / 1e4, np.polyval(k, xs) * apix_bin, "--",
                color="#c0392b", lw=1.6)
    a3.set_xlabel("distance from tilt axis (µm)")
    a3.set_ylabel("median |offset| (Å)")
    a3.set_title("(d) vs radius", fontsize=11, weight="bold")
    a3.grid(alpha=0.25)

    # (e) offset vs bead height -- probes a z-dependent geometry error
    a4 = fig.add_subplot(gs[1, 1])
    a4.scatter(pks[:, 2] * apix_bin / 1e4, pb * apix_bin, s=40, c="#1a8f4c", zorder=3)
    a4.set_xlabel("bead z (µm)"); a4.set_ylabel("median |offset| (Å)")
    a4.set_title("(e) vs height in slab", fontsize=11, weight="bold")
    a4.grid(alpha=0.25)

    # (f) mean offset vector per bead
    a5 = fig.add_subplot(gs[1, 2])
    mx = np.nanmedian(off[:, :, 0], axis=0)
    my = np.nanmedian(off[:, :, 1], axis=0)
    a5.quiver(pks[:, 0] * apix_bin / 1e4, pks[:, 1] * apix_bin / 1e4, mx, my,
              angles="xy", color="#8e44ad", width=0.006)
    a5.scatter(pks[:, 0] * apix_bin / 1e4, pks[:, 1] * apix_bin / 1e4, s=9, c="#333")
    a5.axvline(0, color="#c0392b", ls=":", lw=1.2)
    a5.set_xlabel("x (µm)"); a5.set_ylabel("y (µm)")
    a5.set_title("(f) median offset vector", fontsize=11, weight="bold")
    a5.grid(alpha=0.25); a5.set_aspect("equal", adjustable="datalim")

    wrapped = "\n".join(textwrap.wrap(label, 78))
    fig.suptitle("Fiducial placement: reconstruction vs data\n" + wrapped,
                 fontsize=12.5, weight="bold")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_residuals(res, snr, ttparams, pks, apix, label, out):
    ang = ttparams[:, 3]
    mag = np.hypot(res[:, :, 0], res[:, :, 1])
    per_tilt = np.nanmean(mag, axis=1)
    per_bead = np.nanmean(mag, axis=0)
    # distance of each bead from the tilt axis (the y axis through the centre)
    rad = np.abs(pks[:, 0])

    fig = plt.figure(figsize=(13.6, 9.0))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.28)

    # (a) residual vs tilt angle
    a0 = fig.add_subplot(gs[0, 0])
    for j in range(res.shape[1]):
        a0.plot(ang, mag[:, j] * apix, color="#8899aa", lw=0.6, alpha=0.55)
    a0.plot(ang, per_tilt * apix, color="#c0392b", lw=2.2, label="mean over beads")
    a0.set_xlabel("tilt angle (deg)"); a0.set_ylabel("|residual| (Å)")
    a0.set_title("(a) placement error vs tilt angle", fontsize=11, weight="bold")
    a0.grid(alpha=0.25); a0.legend(fontsize=8)

    # (b) residual vs distance from tilt axis
    a1 = fig.add_subplot(gs[0, 1])
    a1.scatter(rad * apix / 1e4, per_bead * apix, s=42, c="#2f4a6d", zorder=3)
    if np.isfinite(per_bead).sum() > 2:
        g = np.isfinite(per_bead)
        k = np.polyfit(rad[g], per_bead[g], 1)
        xs = np.linspace(rad.min(), rad.max(), 50)
        a1.plot(xs * apix / 1e4, np.polyval(k, xs) * apix, "--", color="#c0392b",
                lw=1.6, label=f"slope {k[0]*apix:.3f} Å/px")
        a1.legend(fontsize=8)
    a1.set_xlabel("distance from tilt axis (µm)")
    a1.set_ylabel("mean |residual| (Å)")
    a1.set_title("(b) placement error vs radius", fontsize=11, weight="bold")
    a1.grid(alpha=0.25)

    # (c) residual vector field, averaged per bead
    a2 = fig.add_subplot(gs[0, 2])
    mx = np.nanmean(res[:, :, 0], axis=0)
    my = np.nanmean(res[:, :, 1], axis=0)
    a2.quiver(pks[:, 0] * apix / 1e4, pks[:, 1] * apix / 1e4, mx, my,
              angles="xy", color="#1a8f4c", width=0.006)
    a2.scatter(pks[:, 0] * apix / 1e4, pks[:, 1] * apix / 1e4, s=9, c="#333")
    a2.axvline(0, color="#c0392b", ls=":", lw=1.2)
    a2.set_xlabel("x (µm)"); a2.set_ylabel("y (µm)")
    a2.set_title("(c) mean residual vector per bead\n(red = tilt axis)",
                 fontsize=11, weight="bold")
    a2.grid(alpha=0.25); a2.set_aspect("equal", adjustable="datalim")

    # (d) x and y components vs angle -- separates the two error modes
    a3 = fig.add_subplot(gs[1, 0])
    a3.plot(ang, np.nanmean(res[:, :, 0], axis=1) * apix, color="#2f4a6d",
            lw=1.8, label="x (⊥ tilt axis)")
    a3.plot(ang, np.nanmean(res[:, :, 1], axis=1) * apix, color="#1a8f4c",
            lw=1.8, label="y (∥ tilt axis)")
    a3.axhline(0, color="#666", lw=0.8)
    a3.set_xlabel("tilt angle (deg)"); a3.set_ylabel("mean residual (Å)")
    a3.set_title("(d) signed components vs angle", fontsize=11, weight="bold")
    a3.grid(alpha=0.25); a3.legend(fontsize=8)

    # (e) scatter of all residuals
    a4 = fig.add_subplot(gs[1, 1])
    fx = res[:, :, 0][np.isfinite(res[:, :, 0])] * apix
    fy = res[:, :, 1][np.isfinite(res[:, :, 1])] * apix
    a4.scatter(fx, fy, s=5, alpha=0.35, c="#2f4a6d")
    a4.axhline(0, color="#666", lw=0.8); a4.axvline(0, color="#666", lw=0.8)
    lim = np.nanpercentile(np.abs(np.r_[fx, fy]), 99) * 1.15
    a4.set_xlim(-lim, lim); a4.set_ylim(-lim, lim)
    a4.set_aspect("equal")
    a4.set_xlabel("x residual (Å)"); a4.set_ylabel("y residual (Å)")
    a4.set_title(f"(e) all residuals  (n={len(fx)})", fontsize=11, weight="bold")
    a4.grid(alpha=0.25)

    # (f) how many beads were located in each tilt
    a5 = fig.add_subplot(gs[1, 2])
    a5.plot(ang, np.isfinite(mag).sum(axis=1), color="#c0392b", lw=1.8)
    a5.set_ylim(0, res.shape[1] + 1)
    a5.set_xlabel("tilt angle (deg)"); a5.set_ylabel("beads located")
    a5.set_title("(f) detection completeness", fontsize=11, weight="bold")
    a5.grid(alpha=0.25)

    wrapped = "\n".join(textwrap.wrap(label, 78))
    fig.suptitle(f"Fiducial placement error from the alignment model\n{wrapped}",
                 fontsize=12.5, weight="bold")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ------------------------------------------------------------------- main ----

def main():
    p = argparse.ArgumentParser(
        description="Reproject a tomogram at the tilt Euler angles and measure "
                    "fiducial placement errors.")
    p.add_argument("--tomo", required=True, help="Tomogram (bin4/SA is fine).")
    p.add_argument("--tilts", required=True, help="Tilt series, same apix as the tomogram.")
    p.add_argument("--ttparams", required=True, help="<name>_ttparams.npy")
    p.add_argument("--pks", required=True, help="<name>_pks.npy")
    p.add_argument("--bin", type=int, default=2,
                   help="Shrink the tomogram before reprojecting (default 2). "
                        "Affects only the pictures, never the bead residuals.")
    p.add_argument("--box", type=int, default=48, help="Bead search box (default 48).")
    p.add_argument("--lowpass", type=float, default=40.0,
                   help="Low-pass both images to this resolution in Å before "
                        "correlating (default 40, matching the tomogram). 0 disables.")
    p.add_argument("--highpass", type=float, default=2000.0,
                   help="High-pass both images at this resolution in Å (default "
                        "2000) to remove slow gradients. 0 disables.")
    p.add_argument("--pscale", type=float, default=0.0,
                   help="Divide ttparams tx/ty and pks by this factor to bring "
                        "them into the tilt series' pixel frame. 0 (default) "
                        "auto-detects from info/<base>_info.json; 1 disables.")
    p.add_argument("--obox", type=int, default=32,
                   help="Box for the reprojection-vs-measured offset "
                        "cross-correlation, in binned pixels (default 32).")
    p.add_argument("--cccrad", type=float, default=0.5,
                   help="Report the CCC over a central disc of this fraction of "
                        "the image half-width (default 0.5); see reproject().")
    p.add_argument("--minsnr", type=float, default=3.0,
                   help="Reject bead peaks weaker than this many sigma (default 3).")
    p.add_argument("--ntile", type=int, default=5, help="Tilts shown in the montage.")
    p.add_argument("--outdir", type=str, default="diagnostics")
    p.add_argument("--label", type=str, default=None)
    p.add_argument("--nostack", action="store_true",
                   help="Do not write the reproj/measured/diff HDF stacks.")
    p.add_argument("--skip_reproject", action="store_true")
    p.add_argument("--skip_beads", action="store_true")
    p.add_argument("--verbose", action="store_true", default=True)
    opt = p.parse_args()

    ttparams = np.load(opt.ttparams)
    pks = np.load(opt.pks)
    ntilt = EMUtil.get_image_count(opt.tilts)

    # --- coordinate frame -------------------------------------------------
    # The refined parameters are stored in the pixel frame of whatever tilt
    # series the alignment ran on. With --patchtrack that is the FULL
    # RESOLUTION series (e.g. 4096 at 3.714 A/pix), while the tomogram and the
    # _Sh4 tilt series are bin4 (1024 at 14.856). Feeding unscaled parameters
    # then puts every predicted bead a factor of 4 out -- beads land outside
    # the volume, almost none are found, and the CCC reads ~0, which looks
    # exactly like a diverged alignment rather than a unit mismatch.
    #
    # e2tomogram_stage2.py does the same rescaling (tlt_params[:,:2] /= scale).
    # Here both the translations and the 3D peak coordinates need it, since
    # get_xf_pos() rotates the peak and then adds the translation.
    scale = opt.pscale
    if scale == 0.0:
        scale = 1.0
        base = os.path.basename(opt.ttparams).replace("_ttparams.npy", "")
        ij = os.path.join("info", base + "_info.json")
        if os.path.isfile(ij):
            try:
                with open(ij) as fh:
                    unb = json.load(fh).get("apix_unbin")
                if unb:
                    scale = EMData(opt.tilts, 0, True)["apix_x"] / float(unb)
            except Exception as e:
                print(f"  could not read {ij} ({e}); assuming pscale=1")
    if abs(scale - 1.0) > 1e-6:
        print(f"  rescaling tilt parameters and peaks by 1/{scale:.4f} "
              f"(alignment ran at a finer pixel size than this tilt series)")
        ttparams = ttparams.copy()
        ttparams[:, :2] /= scale
        pks = pks / scale
    if len(ttparams) != ntilt:
        print(f"Error: {len(ttparams)} tilt parameter rows but {ntilt} tilt images.")
        return 1

    h = EMData(opt.tilts, 0, True)
    apix = h["apix_x"]
    label = opt.label or os.path.basename(opt.tomo).replace(".hdf", "")
    # Dataset identity for the shared 'measured' stack: taken from the tilt
    # parameter file, which names the dataset and not the reconstruction.
    stackbase = os.path.basename(opt.ttparams).replace("_ttparams.npy", "")
    os.makedirs(opt.outdir, exist_ok=True)

    print(f"{ntilt} tilts, {len(pks)} beads, {h['nx']}x{h['ny']} at {apix:.3f} Å/pix")
    print(f"  tilt range {ttparams[:,3].min():+.2f} .. {ttparams[:,3].max():+.2f} deg")

    if not opt.skip_beads:
        print("\nMeasuring bead positions...")
        res, snr = measure_beads(opt.tilts, ttparams, pks, opt.box,
                                 opt.minsnr, opt.verbose)
        mag = np.hypot(res[:, :, 0], res[:, :, 1])
        found = np.isfinite(mag)
        print(f"\n  located {found.sum()}/{found.size} bead-in-tilt instances "
              f"({100.0*found.sum()/found.size:.1f}%)")
        if found.sum():
            print(f"  mean   |residual| = {np.nanmean(mag)*apix:8.2f} Å "
                  f"({np.nanmean(mag):.3f} px)")
            print(f"  median |residual| = {np.nanmedian(mag)*apix:8.2f} Å")
            print(f"  90th pct          = {np.nanpercentile(mag,90)*apix:8.2f} Å")
            print(f"  RMS               = {np.sqrt(np.nanmean(mag**2))*apix:8.2f} Å")
            out = os.path.join(opt.outdir, f"FidResid_{label}.png")
            plot_residuals(res, snr, ttparams, pks, apix, label, out)
            txt = os.path.join(opt.outdir, f"FidResid_{label}.txt")
            with open(txt, "w") as fh:
                fh.write("# tilt  ytilt_deg  bead  dx_px  dy_px  |d|_A  snr\n")
                for i in range(ntilt):
                    for j in range(len(pks)):
                        if np.isfinite(mag[i, j]):
                            fh.write(f"{i:5d} {ttparams[i,3]:10.3f} {j:5d} "
                                     f"{res[i,j,0]:8.3f} {res[i,j,1]:8.3f} "
                                     f"{mag[i,j]*apix:9.2f} {snr[i,j]:7.2f}\n")
            print(f"  wrote {txt}")

    if not opt.skip_reproject:
        print("\nReprojecting...")
        ccc, ccc_full, off, panels = reproject(
            opt.tomo, opt.tilts, ttparams, pks, opt.bin, opt.lowpass,
            opt.highpass, opt.cccrad, opt.obox, opt.verbose,
            None if opt.nostack else opt.outdir, label, stackbase, scale,
            opt.ttparams)
        print(f"  CCC (central r<{opt.cccrad}): mean {np.nanmean(ccc):.4f}  "
              f"min {np.nanmin(ccc):.4f}  max {np.nanmax(ccc):.4f}")
        print(f"  CCC (full frame)           : mean {np.nanmean(ccc_full):.4f}")
        plot_montage(panels, ttparams, ccc, label,
                     os.path.join(opt.outdir, f"FidReproj_{label}.jpg"), opt.ntile)
        fig, ax = plt.subplots(figsize=(8.8, 5.0))
        ax.plot(ttparams[:, 3], ccc, "o-", color="#2f4a6d", lw=1.8, ms=4,
                label=f"central disc (r < {opt.cccrad}·half-width)")
        ax.plot(ttparams[:, 3], ccc_full, "s--", color="#999", lw=1.2, ms=3,
                label="full frame")
        ax.axhline(0, color="#666", lw=0.8)
        ax.set_xlabel("tilt angle (deg)")
        ax.set_ylabel("CCC (reprojection vs measured)")
        ax.set_title(f"Reprojection agreement per tilt\n{label}",
                     fontsize=11.5, weight="bold")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8.5)
        out = os.path.join(opt.outdir, f"FidReprojCCC_{label}.png")
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out}")
        np.savetxt(os.path.join(opt.outdir, f"FidReprojCCC_{label}.txt"),
                   np.c_[ttparams[:, 3], ccc, ccc_full], fmt="%10.4f",
                   header="ytilt_deg   ccc_central   ccc_full")

        # --- placement error: reconstruction vs data -------------------------
        apix_bin = apix * opt.bin
        mag = np.hypot(off[:, :, 0], off[:, :, 1])
        good = np.isfinite(mag)
        zi = int(np.argmin(np.abs(ttparams[:, 3])))
        print(f"\n  Placement offset (reprojection vs measured), "
              f"{good.sum()}/{good.size} beads:")
        print(f"    median  {np.nanmedian(mag)*apix_bin:8.1f} Å "
              f"({np.nanmedian(mag):.2f} px at {apix_bin:.2f} Å/pix)")
        print(f"    90th    {np.nanpercentile(mag[good],90)*apix_bin:8.1f} Å")
        for nm, sl in (("negative half", slice(0, zi + 1)),
                       ("positive half", slice(zi + 1, len(ttparams)))):
            s = mag[sl][np.isfinite(mag[sl])]
            dy = np.nanmedian(off[sl, :, 1])
            if len(s):
                print(f"    {nm}: median {np.median(s)*apix_bin:7.1f} Å   "
                      f"median dy {dy*apix_bin:+7.1f} Å")
        plot_offsets(off, ttparams, pks, apix_bin, label,
                     os.path.join(opt.outdir, f"FidOffset_{label}.png"), zi)
        with open(os.path.join(opt.outdir, f"FidOffset_{label}.txt"), "w") as fh:
            fh.write(f"# reprojection -> measured bead offsets, "
                     f"{apix_bin:.4f} A/pix (bin {opt.bin})\n")
            fh.write("# tilt  ytilt_deg  bead   dx_px    dy_px    |d|_A\n")
            for i in range(len(ttparams)):
                for j in range(len(pks)):
                    if np.isfinite(mag[i, j]):
                        fh.write(f"{i:5d} {ttparams[i,3]:10.3f} {j:5d} "
                                 f"{off[i,j,0]:8.3f} {off[i,j,1]:8.3f} "
                                 f"{mag[i,j]*apix_bin:9.2f}\n")
        print(f"  wrote {os.path.join(opt.outdir, f'FidOffset_{label}.txt')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
