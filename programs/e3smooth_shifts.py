#!/usr/bin/env python
"""Re-solve a tilt-series alignment with a preference for smooth image shifts.

THE PROBLEM THIS SOLVES
-----------------------
At high tilt the sideways shift tx and the landmarks' depth z stop being
distinguishable: moving the picture sideways and moving the beads deeper produce
almost the same predicted bead positions, and the ambiguity grows with sin(tilt).
The consequence is not a single correct alignment but a family of them, all
fitting the observed beads about equally well -- some smooth, some kinked.

e2tomogram.py refines each tilt independently with nothing tying it to its
neighbours, so it can and does switch between members of that family part-way
through a series. Measured on Fiducial4 (08419/08417): the fit follows one
branch to tilt 117, hops to another at tilt 118 (a 39 px step in tx while ty
stays perfectly smooth), and falls back by tilt 123. The beads still fit -- 19
of 20 found, residuals among the best in the series -- because the hopped-to
solution is genuinely valid. It simply looks wrong, and it is wrong in the sense
that the specimen did not move that way.

e2tomogram's own --xdrift avoids the hop by walking outward from the centre and
allowing only +/-5 whole pixels of change per tilt. That works, but the coarse
integer search adds shimmer everywhere: on Fiducial4 the frame-to-frame jitter
rose from 2.68 px to 3.66 px, degrading the low tilts that were already correct.

WHAT THIS DOES INSTEAD
----------------------
Starts from the alignment e2tomogram already produced and re-solves the whole
series at once, minimising

    (bead reprojection error)  +  lambda * (second difference of tx)

continuously, with no integer grid. Where the data pins tx down -- near zero
tilt -- the data term is stiff and the penalty barely moves it. Where the data
is ambiguous -- high tilt -- the data term is flat along the degenerate
direction, so the penalty picks the smooth member of the family. The correction
therefore concentrates itself exactly where the problem is, without touching
tilts that were already right.

Nothing here modifies e2tomogram.py. It reads that program's own output JSON,
writes a corrected one, and the tomogram is rebuilt with e2tomogram_stage2.py,
which reconstructs from supplied parameters and performs no fitting.

GAUGE
-----
The ambiguity means the solution is not unique even after smoothing: shifting
every landmark's z by d and every tx by -d*sin(tilt) leaves both terms nearly
unchanged. The mean landmark position is therefore pinned to its input value so
the optimiser cannot drift along that direction.
"""

import os, sys, json, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e3fid_reproject import get_xf_pos, measure_beads          # noqa: E402


def projection_matrices(tlt_params):
    """2x3 matrix R_i with  predicted = R_i @ pk + [tx_i, ty_i].

    get_xf_pos applies the rotation then adds (tx,ty), so the map is exactly
    linear in the landmark coordinate. Verified against get_xf_pos over 1500
    random cases to 1.1e-4 px, i.e. to EMAN's single precision.
    """
    R = np.zeros((len(tlt_params), 2, 3))
    for i, t in enumerate(tlt_params):
        rot = [0.0, 0.0, t[2], t[3], t[4]]
        for k, e in enumerate(np.eye(3)):
            R[i, :, k] = get_xf_pos(rot, e)
    return R


def pack(tx, ty, pks):
    return np.concatenate([tx, ty, pks.ravel()])


def unpack(v, n, m):
    return v[:n], v[n:2 * n], v[2 * n:].reshape(m, 3)


def build_problem(n, m, R, meas, mask, lam, lam_y, w, pks0_mean, gauge_w):
    """Residual function and its exact Jacobian.

    Every residual is LINEAR in the unknowns (tx, ty, landmark positions):
    the projection is  R_i @ pk + [tx_i, ty_i], the smoothness term is a second
    difference, and the gauge is a mean. The Jacobian is therefore constant, so
    it is built once as a sparse matrix. Supplying it turns each optimiser step
    from 316 numerical probes into one sparse solve, which is what makes a
    lambda sweep practical.
    """
    from scipy.sparse import coo_matrix
    ii, jj = np.nonzero(mask)
    nd = len(ii)
    nrow = 2 * nd + 2 * (n - 2) + 3
    ncol = 2 * n + 3 * m

    rows, cols, vals = [], [], []
    for k in range(nd):
        i_, j_ = ii[k], jj[k]
        rx, ry = 2 * k, 2 * k + 1
        rows += [rx, ry];               cols += [i_, n + i_];        vals += [1.0, 1.0]
        for c in range(3):
            rows += [rx, ry]
            cols += [2 * n + 3 * j_ + c, 2 * n + 3 * j_ + c]
            vals += [R[i_, 0, c], R[i_, 1, c]]
    base = 2 * nd
    for r in range(n - 2):
        cf = lam * w[r + 1]
        rows += [base + r] * 3; cols += [r, r + 1, r + 2]; vals += [cf, -2 * cf, cf]
        cf = lam_y * w[r + 1]
        rows += [base + (n - 2) + r] * 3
        cols += [n + r, n + r + 1, n + r + 2]; vals += [cf, -2 * cf, cf]
    base = 2 * nd + 2 * (n - 2)
    for g in range(3):
        for j_ in range(m):
            rows.append(base + g); cols.append(2 * n + 3 * j_ + g)
            vals.append(gauge_w / m)
    J = coo_matrix((vals, (rows, cols)), shape=(nrow, ncol)).tocsr()

    def resid(v):
        tx, ty, pks = unpack(v, n, m)
        pred = np.einsum('ikl,jl->ijk', R, pks)
        pred[:, :, 0] += tx[:, None]
        pred[:, :, 1] += ty[:, None]
        data = (pred[ii, jj] - meas[ii, jj]).ravel()
        d2x = tx[2:] - 2 * tx[1:-1] + tx[:-2]
        d2y = ty[2:] - 2 * ty[1:-1] + ty[:-2]
        sm = np.concatenate([lam * w[1:-1] * d2x, lam_y * w[1:-1] * d2y])
        gauge = gauge_w * (pks.mean(axis=0) - pks0_mean)
        return np.concatenate([data, sm, gauge])

    return resid, (lambda v: J)


def report(tag, tx, ty, R, pks, meas, mask):
    ii, jj = np.nonzero(mask)
    pred = np.einsum('ikl,jl->ijk', R, pks)
    pred[:, :, 0] += tx[:, None]
    pred[:, :, 1] += ty[:, None]
    d = np.linalg.norm(pred[ii, jj] - meas[ii, jj], axis=1)
    step = np.hypot(np.diff(tx), np.diff(ty))
    k = int(np.argmax(step))
    print(f"  {tag:22s} bead err median {np.median(d):6.3f} px   "
          f"mean {d.mean():6.3f} px |  step median {np.median(step):5.2f} px   "
          f"max {step.max():6.2f} px at {k}->{k+1}")
    return np.median(d), np.median(step), step.max()


def main():
    ap = argparse.ArgumentParser(
        description="Re-solve a tilt alignment preferring smooth image shifts.")
    ap.add_argument("json", help="e2tomogram --stage2prep record, "
                                 "info/<base>_recon_tomo_final.json")
    ap.add_argument("--out", default="", help="output JSON (default <in>_smooth.json)")
    ap.add_argument("--lam", type=float, default=3.0,
                    help="smoothness weight for tx (default 3.0)")
    ap.add_argument("--lamy", type=float, default=-1.0,
                    help="smoothness weight for ty (default: same as --lam)")
    ap.add_argument("--anglepower", type=float, default=0.0,
                    help="weight smoothing by |sin(tilt)|**p. 0 = uniform "
                         "(default); the data term already stiffens at low tilt.")
    ap.add_argument("--box", type=int, default=48, help="bead search box (default 48)")
    ap.add_argument("--minsnr", type=float, default=3.0, help="bead SNR floor")
    ap.add_argument("--gauge", type=float, default=100.0,
                    help="weight pinning the mean landmark position")
    ap.add_argument("--tilts", default="", help="override the tilt series path")
    ap.add_argument("--meas", default="",
                    help="cache file for the measured bead positions. Measuring "
                         "is the slow step and is identical for every smoothing "
                         "weight, so a lambda sweep should reuse one cache.")
    ap.add_argument("--iters", type=int, default=3,
                    help="outer iterations of measure->optimise (default 3). "
                         "One pass cannot escape a bad alignment: bead positions "
                         "are measured by searching near where the CURRENT "
                         "alignment predicts them, so the data handed to the "
                         "optimiser already assumes the answer. Alternating lets "
                         "the solution migrate.")
    ap.add_argument("--junction", type=int, default=-2,
                    help="index of the last tilt of the negative acquisition. The "
                         "shift genuinely jumps across the neg/pos boundary -- the "
                         "stage returns and the specimen settles elsewhere -- so "
                         "smoothing must NOT span it or it erodes a real offset "
                         "(seen on Fiducial4: 48 px shrinking to 32 px over "
                         "iterations). -2 = auto-detect from the tilt angles, "
                         "-1 = disable (treat the series as continuous).")
    ap.add_argument("--stoptol", type=float, default=1.25,
                    help="stop when the bead fit gets worse than this multiple "
                         "of the best seen (default 1.25), and return the best "
                         "iterate. The alternating scheme is not stable "
                         "indefinitely: once the alignment drifts far enough the "
                         "bead search latches onto the wrong features and the fit "
                         "collapses abruptly. Measured on Fiducial4 at lam=30, the "
                         "median bead error ran 1.16 -> 1.13 px over iterations "
                         "6-16 and then jumped to 5.63 px at iteration 17.")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    with open(a.json) as f:
        rec = json.load(f)
    tlt = np.array(rec["tlt_params"], dtype=float)
    pks = np.array(rec["pks"], dtype=float)
    tilts = a.tilts or rec["inputname"]
    n, m = len(tlt), len(pks)
    print(f"=== {os.path.basename(a.json)} ===")
    print(f"  {n} tilts, {m} landmarks, tilt series {tilts}")
    if not os.path.isfile(tilts):
        sys.exit(f"Error: tilt series not found: {tilts}")

    R = projection_matrices(tlt)          # rotations are held fixed throughout
    lam_y = a.lam if a.lamy < 0 else a.lamy
    # locate the neg/pos junction and exempt it from the smoothness penalty
    zi = a.junction
    if zi == -2:
        # Prefer the chunk record: it states how many tilts came from the
        # negative acquisition, which is the actual boundary. Do NOT infer it
        # from a sign change in ytilt -- the zero-reference tilt is written as
        # exactly 0.00, so the sign flips one tilt early (Fiducial4: sign flips
        # after 64, the true boundary is after 65).
        zi = None
        bn = rec.get("basename", "")
        cj = os.path.join("info", f"sa_chunks_{bn}.json")
        if bn and os.path.isfile(cj):
            try:
                ch = json.load(open(cj))
                nneg = sum(1 for c in ch["chunks"] if c[0] < ch["n_sel_neg"])
                zi = nneg - 1
                print(f"  neg/pos junction after tilt {zi} (from {cj}: "
                      f"{nneg} negative tilts)")
            except Exception as e:
                print(f"  could not read {cj} ({e})")
        if zi is None:
            sgn = np.sign(tlt[:, 3])
            flips = np.nonzero(np.diff(sgn) != 0)[0]
            zi = int(flips[0]) if len(flips) else int(np.argmin(np.abs(tlt[:, 3])))
            print(f"  WARNING: no chunk record; guessing junction after tilt {zi} "
                  f"from the angle sign change. Verify this.")
    elif zi >= 0:
        print(f"  neg/pos junction given as after tilt {zi}")
    else:
        print("  no junction: series treated as continuous")

    w = np.ones(n)
    if a.anglepower > 0:
        w = np.abs(np.sin(np.deg2rad(tlt[:, 3]))) ** a.anglepower
        w /= max(w.max(), 1e-9)
    if zi >= 0:
        # a second difference centred at i spans tilts i-1,i,i+1; zero the two
        # that straddle the junction so the real offset survives untouched
        for c in (zi, zi + 1):
            if 1 <= c <= n - 2:
                w[c] = 0.0

    from scipy.optimize import least_squares
    cur = tlt.copy()
    pks_cur = pks.copy()
    pks0_mean = pks.mean(axis=0)

    best = {"err": np.inf, "tlt": cur.copy(), "pks": pks_cur.copy(), "it": 0}
    for it in range(a.iters):
        print(f"\n--- iteration {it+1}/{a.iters} ---")
        use_cache = (it == 0 and a.meas and os.path.isfile(a.meas))
        res = None
        if use_cache:
            z = np.load(a.meas)
            if (z["res"].shape == (n, m, 2) and np.allclose(z["tlt"], cur, atol=1e-9)
                    and np.allclose(z["pks"], pks_cur, atol=1e-9)
                    and float(z["box"]) == a.box and float(z["minsnr"]) == a.minsnr):
                res = z["res"]
                print(f"  reusing measured positions from {a.meas}")
        if res is None:
            print("  measuring bead positions with the current alignment...")
            res, _snr = measure_beads(tilts, cur, pks_cur, a.box, a.minsnr, a.verbose)
            if it == 0 and a.meas:
                np.savez(a.meas, res=res, snr=_snr, tlt=cur, pks=pks_cur,
                         box=a.box, minsnr=a.minsnr)

        mask = np.isfinite(res[:, :, 0]) & np.isfinite(res[:, :, 1])
        print(f"  located {mask.sum()}/{n*m} bead-in-tilt instances "
              f"({100*mask.sum()/(n*m):.1f}%)")

        pred = np.einsum('ikl,jl->ijk', R, pks_cur)
        pred[:, :, 0] += cur[:, 0][:, None]
        pred[:, :, 1] += cur[:, 1][:, None]
        meas = pred + res

        report("  before", cur[:, 0], cur[:, 1], R, pks_cur, meas, mask)
        fun, jac = build_problem(n, m, R, meas, mask, a.lam, lam_y, w,
                                 pks0_mean, a.gauge)
        v0 = pack(cur[:, 0].copy(), cur[:, 1].copy(), pks_cur.copy())
        sol = least_squares(fun, v0, jac=jac, method="trf", loss="soft_l1",
                            f_scale=2.0, tr_solver="lsmr", max_nfev=400, verbose=0)
        tx, ty, pks_cur = unpack(sol.x, n, m)
        cur[:, 0], cur[:, 1] = tx, ty
        e_after, _, _ = report("  after ", tx, ty, R, pks_cur, meas, mask)
        d = np.abs(np.diff(tx))
        print(f"  worst step (excl. the neg/pos seam): "
              f"{np.delete(d, min(int(np.argmax(d)), len(d)-1)).max():.1f} px"
              f"    median step {np.median(d):.2f} px")

        # the fit against FRESH measurements is the honest score for this iterate
        if e_after < best["err"]:
            best = {"err": e_after, "tlt": cur.copy(), "pks": pks_cur.copy(),
                    "it": it + 1}
            print(f"  best so far (bead err {e_after:.3f} px)")
        elif e_after > best["err"] * a.stoptol:
            print(f"  bead fit degraded to {e_after:.3f} px, over {a.stoptol:g}x the "
                  f"best ({best['err']:.3f} px at iteration {best['it']}). Stopping.")
            break

    if best["it"] and (best["it"] != a.iters):
        print(f"\n  returning iteration {best['it']} (bead err {best['err']:.3f} px)")
    cur = best["tlt"]; pks_cur = best["pks"]

    tlt_out = cur
    tx, ty = cur[:, 0], cur[:, 1]
    pks_new = pks_cur
    dd = np.abs(tx - tlt[:, 0])
    print(f"\n  tx changed from input by: median {np.median(dd):.2f} px, "
          f"max {dd.max():.2f} px at tilt {int(np.argmax(dd))}")
    print(f"  landmark mean moved {np.linalg.norm(pks_new.mean(0)-pks0_mean):.3f} px")

    out = a.out or a.json.replace(".json", "_smooth.json")
    rec2 = dict(rec)
    rec2["tlt_params"] = tlt_out.tolist()
    rec2["pks"] = pks_new.tolist()
    rec2["smoothed_by"] = {"tool": "e3smooth_shifts.py", "lam": a.lam,
                           "lamy": lam_y, "anglepower": a.anglepower,
                           "iters": a.iters, "source": os.path.abspath(a.json)}
    with open(out, "w") as f:
        json.dump(rec2, f, indent=1)
    print(f"\n  wrote {out}")
    print(f"  rebuild with:  e2tomogram_stage2.py {out} <label> --filterres=40.0")


if __name__ == "__main__":
    main()
