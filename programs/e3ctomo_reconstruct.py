#!/usr/bin/env python
# e3ctomo_reconstruct.py  --  Continuous-tilt reconstruction launcher.
from past.utils import old_div
from builtins import range
from EMAN2 import *
import sys, os, json, subprocess
import numpy as np

# Pre-load pick_ranges.json so parser defaults can reference pick_ranges_data
# (also injected by e2projectmanager.py's eval namespace)
pick_ranges_data = {}
try:
    _pr = os.path.join(os.getcwd(), "info", "pick_ranges.json")
    if os.path.isfile(_pr):
        with open(_pr) as _f:
            pick_ranges_data = json.load(_f)
except Exception:
    pass


# ---------------------------------------------------------------------------
# Per-SA-tilt angles, derived from the *selected* frame ranges of each half.
#
# The stage pauses at the start of each movie while frames are still being
# acquired; e3pick_ranges.py trims those away, so each half's retained range
# spans exactly that half's nominal sweep -- the angle ramp maps onto the
# SELECTED range, not the raw movie. The two halves are separate acquisitions
# and may sweep at different rates, which a single uniform --tltstep cannot
# express (forcing one rate onto a two-rate scan inflates the apparent angular
# range), so we hand e2tomogram.py an explicit --rawtlt angle list instead.
# ---------------------------------------------------------------------------

def _compute_sa_tilt_angles(sa_path, neg_start, neg_end, pos_start, pos_end, avgseq_sa,
                            neg_start_angle, neg_end_angle, pos_start_angle, pos_end_angle):
    """Return (angles, zeroid) for every tilt of the SA tiltseries.

    Mirrors e3pick_ranges.py's combined-LST convention: the neg half is stored
    reversed, so combined index neg_start holds the *last*-acquired neg frame
    (neg_end_angle) and neg_end the first retained one (neg_start_angle); the
    pos half is direct, pos_start->pos_start_angle and pos_end->pos_end_angle.
    The final LST is [selected neg]+[selected pos], and the SA tiltseries
    groups every avgseq_sa consecutive LST entries into one output tilt.

    Note the two halves can overlap in angle (e.g. neg sweeping +2 -> -60 and
    pos 0 -> +60 both cover 0..+2), so the returned list is not necessarily
    monotonic; e2tomogram.py argsorts it for coarse alignment.
    """
    n_sa = EMUtil.get_image_count(sa_path)
    n_sel_neg = neg_end - neg_start + 1
    n_sel_pos = pos_end - pos_start + 1
    total = n_sel_neg + n_sel_pos

    def _angle_of_lst(j):
        if j < n_sel_neg:
            f = j / float(neg_end - neg_start)
            return neg_end_angle + f * (neg_start_angle - neg_end_angle)
        f = (j - n_sel_neg) / float(pos_end - pos_start)
        return pos_start_angle + f * (pos_end_angle - pos_start_angle)

    # Chunk boundaries. e3pick_ranges.py subaverages each half separately so no
    # chunk straddles the neg/pos boundary, and records the exact LST row ranges.
    # Use them when present; otherwise fall back to uniform chunking over the
    # combined list, which is how older series were built.
    chunks = None
    _bn = os.path.splitext(os.path.basename(sa_path))[0]
    for _cand in (os.path.join("info", f"sa_chunks_{_bn}.json"),
                  os.path.join("info", f"sa_chunks_{_bn.split('__')[0]}.json")):
        if os.path.isfile(_cand):
            try:
                with open(_cand) as _f:
                    _d = json.load(_f)
                if _d.get("avgseq_sa") == avgseq_sa and len(_d.get("chunks", [])) == n_sa:
                    chunks = _d["chunks"]
                    print(f"  using chunk boundaries from {_cand} "
                          f"(scheme={_d.get('scheme')}, no chunk crosses the seam)")
                else:
                    print(f"  {_cand} does not match this series "
                          f"({len(_d.get('chunks', []))} chunks vs {n_sa} images, "
                          f"avgseq {_d.get('avgseq_sa')} vs {avgseq_sa}); ignoring")
            except Exception as _e:
                print(f"  could not read {_cand} ({_e}); using uniform chunking")
            break
    if chunks is None:
        print("  no chunk record found; assuming uniform chunking over the "
              "combined list (a chunk may straddle the neg/pos seam)")
        chunks = [[min(k * avgseq_sa, total - 1),
                   min((k + 1) * avgseq_sa - 1, total - 1)] for k in range(n_sa)]

    angles = np.zeros(n_sa)
    for k in range(n_sa):
        j0, j1 = chunks[k]
        angles[k] = 0.5 * (_angle_of_lst(j0) + _angle_of_lst(j1))

    zeroid = int(np.argmin(np.abs(angles)))
    negrate = abs(neg_start_angle - neg_end_angle) / (neg_end - neg_start) * avgseq_sa
    posrate = abs(pos_end_angle - pos_start_angle) / (pos_end - pos_start) * avgseq_sa
    print(f"  tilt angles: {n_sa} SA tilts, {angles.min():+.2f} to {angles.max():+.2f} deg "
          f"(neg {negrate:.5f}, pos {posrate:.5f} deg/SAtilt) -> zeroid={zeroid} "
          f"(angle {angles[zeroid]:+.3f})")
    if max(negrate, posrate) > 0 and abs(negrate - posrate) / max(negrate, posrate) > 0.05:
        print(f"  WARNING: neg/pos sweep rates differ by "
              f"{abs(negrate - posrate) / max(negrate, posrate) * 100:.1f}% -- check the angle "
              f"inputs and the Pick Ranges selection for each half.")
    return angles, zeroid


def _write_rawtlt(angles, sa_path):
    """Write the per-tilt angle list that e2tomogram.py --rawtlt expects."""
    out = os.path.join("info", os.path.splitext(os.path.basename(sa_path))[0] + ".rawtlt")
    os.makedirs("info", exist_ok=True)
    np.savetxt(out, angles, fmt="%.6f")
    print(f"  wrote {out}")
    return out


# ---------------------------------------------------------------------------
# Path helpers  (all accept optional basename; fall back to auto-derived names)
# ---------------------------------------------------------------------------

def _fullsh4_path(lst_path, basename=None, meanshrink=4):
    tag = f"Sh{meanshrink}"
    if basename:
        return os.path.join("sets", f"{basename}_Full_{tag}.hdf")
    stem = os.path.splitext(os.path.basename(lst_path))[0]
    return os.path.join("sets", f"{stem}_{tag}.hdf")

def _tomogram_json(sa_path):
    """Path that e2tomogram.py --stage2prep writes: info/<sa_base>_recon_tomo_final.json."""
    stem = os.path.splitext(os.path.basename(sa_path))[0]
    first = stem.split("__")[0]
    return os.path.join("info", f"{first}_recon_tomo_final.json")

def _stage2_label(basename=None, interp_json=None):
    """Label passed to e2tomogram_stage2.py as baseNameForOutputs."""
    if basename:
        return basename
    if interp_json:
        import re
        stem = os.path.splitext(os.path.basename(interp_json))[0]
        stem = re.sub(r"_Full_Sh\d+$", "", stem)
        return stem
    return "stage2"

def _interp_json(fullsh4_path, basename=None, meanshrink=4):
    tag = f"Sh{meanshrink}"
    if basename:
        return os.path.join("info", f"{basename}_Full_{tag}.json")
    stem = os.path.splitext(os.path.basename(fullsh4_path))[0]
    return os.path.join("info", f"{stem}.json")

def _interp_log(fullsh4_path, basename=None, meanshrink=4):
    tag = f"Sh{meanshrink}"
    if basename:
        return os.path.join("info", f"{basename}_Full_{tag}.txt")
    stem = os.path.splitext(os.path.basename(fullsh4_path))[0]
    return os.path.join("info", f"{stem}.txt")


# ---------------------------------------------------------------------------
# Pipeline steps  (all paths passed in; no internal path derivation)
# ---------------------------------------------------------------------------

def _run_logged(cmd_str, log_path):
    """Run cmd_str, streaming stdout+stderr live to log_path (unbuffered)."""
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    with open(log_path, "w") as lf:
        proc = subprocess.Popen(cmd_str, shell=True, stdout=lf,
                                stderr=subprocess.STDOUT, env=env)
        return proc.wait()


def run_tomogram(sa_path, zeroid, rawtlt_path, threads, clipz, tltkeep, extrapad=True, moretile=False, keeptmp=False, patchtrack=-1, bxsz=32, filterres=40.0, fixytilt=False, niter="2,1,1,1", tltax=0.5, npk=20, compressbits=6, correctrot=False, posz=False, xdrift=False, autoclipxy=False):
    """Run e2tomogram.py --stage2prep: reconstruct and write rich JSON."""
    n_sa = EMUtil.get_image_count(sa_path)
    cmd = (
        f"e2tomogram.py {sa_path}"
        f" --stage2prep"
        f" --outsuffix=_SA"
        f" --tltax={tltax}"
        f" --zeroid={zeroid}"
        f" --rawtlt={rawtlt_path}"
        f" --npk={npk}"
        f" --tltkeep={tltkeep}"
        f" --outsize=1k"
        f" --niter={niter}"
        f" --bytile"
        f" --pkkeep=0.9"
        f" --compressbits={compressbits}"
        f" --clipz={clipz}"
        f" --bxsz={bxsz}"
        f" --filterres={filterres}"
        f" --rmbeadthr=-1.0"
        f" --threads={threads}"
        f" --patchtrack={patchtrack}"
    )
    if not keeptmp:
        cmd += " --notmp"
    if extrapad:
        cmd += " --extrapad"
    if moretile:
        cmd += " --moretile"
    if fixytilt:
        cmd += " --fixytilt"
    # Steve Ludtke's platelet recommendation: patch tracking alone
    # (--niter=0 --patchtrack=2 --correctrot --posz --xdrift --autoclipxy).
    # Platelets lack the strong localized features landmark alignment needs, so
    # refinement has nothing to lock onto. Note fixytilt is inert when niter=0,
    # since there is no refinement to constrain.
    for _flag, _on in (("correctrot", correctrot), ("posz", posz),
                       ("xdrift", xdrift), ("autoclipxy", autoclipxy)):
        if _on:
            cmd += f" --{_flag}"
    print(f"\n=== Step 3 — tomogram reconstruction + stage2prep ({n_sa} SA frames, zeroid={zeroid}, rawtlt={rawtlt_path}) ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Error: e2tomogram.py exited with code {ret}")
        sys.exit(ret)


def _copy_tlt_keys_to_sa_info(tomo_json, basename):
    """Move tlt_file/tlt_params/ali_loss to the info file read by the SA tomogram.

    e2tomogram.py always writes these keys to info/<basename>_info.json (derived from
    the input tiltseries, which splits on '__').

    Two output cases based on input resolution:
      - 1K input  → output: tomograms/<basename>_SA.hdf
                    info_name: info/<basename>_SA_info.json  (different from src)
                    → MOVE keys there; stage-2 tomogram (also reads _info.json) loses AliTilts
      - 4K input (patchtrack) → output: tomograms/<basename>__bin4_SA.hdf
                    info_name: info/<basename>_info.json  (same as src — splits on __)
                    → do nothing; keys are already in the right file
    """
    src_path      = os.path.join("info", f"{basename}_info.json")
    bin4_tomo     = os.path.join("tomograms", f"{basename}__bin4_SA.hdf")

    if os.path.isfile(bin4_tomo):
        # Patchtrack / 4K case: output info file == src; keys already in the right place.
        print(f"  Patchtrack output {bin4_tomo}: tlt keys remain in {src_path}")
        return

    # Normal 1K case: move keys → _SA_info.json so stage-2 tomogram loses AliTilts.
    sa_path = os.path.join("info", f"{basename}_SA_info.json")
    if not os.path.isfile(src_path) or not os.path.isfile(sa_path):
        return
    try:
        with open(src_path) as f:
            src = json.load(f)
        with open(sa_path) as f:
            sa = json.load(f)
        moved = []
        for key in ("tlt_file", "tlt_params", "ali_loss"):
            if key in src:
                sa[key] = src.pop(key)
                moved.append(key)
        with open(sa_path, "w") as f:
            json.dump(sa, f, indent=2)
        with open(src_path, "w") as f:
            json.dump(src, f, indent=2)
        print(f"  Moved {moved} → {sa_path}")
    except Exception as e:
        print(f"  Warning: could not update SA info: {e}")


def _read_sa_chunks(basename, avgseq, n_sel_neg, n_tot):
    """Chunk boundaries recorded by e3pick_ranges.py, or None.

    The shrunk SA series is grouped by _even_chunks(): each half is split at the
    acquisition boundary and its leftover distributed evenly, so every tilt holds
    avgseq +/- 1 frames. Recomputing that here would be a second implementation
    of the same rule, free to drift from the first. Reading the record instead
    makes the full-resolution series match the shrunk one chunk for chunk.

    Without this the full-res series is built with a fixed stride and the two
    disagree on the number of tilts -- measured on 4 of the 5 RatNeuron sets
    (ct1 51 vs 52, ct2 57 vs 58, ct3 58 vs 60, ct5 57 vs 59).
    """
    if not basename:
        return None
    cand = os.path.join("info", f"sa_chunks_{basename}.json")
    if not os.path.isfile(cand):
        return None
    try:
        with open(cand) as f:
            d = json.load(f)
    except Exception as e:
        print(f"  could not read {cand} ({e})")
        return None
    if d.get("avgseq_sa") != avgseq:
        print(f"  {cand}: avgseq {d.get('avgseq_sa')} != {avgseq}; ignoring")
        return None
    if d.get("n_sel_neg") != n_sel_neg:
        print(f"  {cand}: n_sel_neg {d.get('n_sel_neg')} != {n_sel_neg}; ignoring")
        return None
    ch = d.get("chunks") or []
    if not ch or ch[-1][1] != n_tot - 1:
        print(f"  {cand}: chunks end at {ch[-1][1] if ch else '?'}, expected "
              f"{n_tot-1}; ignoring")
        return None
    print(f"  using chunk boundaries from {cand} ({len(ch)} tilts, scheme="
          f"{d.get('scheme')})")
    return ch


def make_fullres_sa(lst_path, sa_fullres_path, avgseq, compressbits, n_sel_neg=0,
                    basename=None):
    """Build a full-resolution (meanshrink=1) SA tilt series for patchtracking.

    Unlike the standard SA series (built in Pick Ranges with spatial shrink),
    this keeps the full camera resolution so patchtracking has enough pixels.
    Output: tiltseries/<basename>__SA{avgseq}.hdf

    MUST match the shrunk series tilt for tilt, since patchtracking refines an
    alignment derived from that series. The grouping therefore comes from
    info/sa_chunks_<basename>.json, written by e3pick_ranges.py, rather than
    being recomputed here. See _read_sa_chunks().

    Falls back to the older per-half fixed-stride behaviour when no record is
    available, and says so -- that path can disagree with the shrunk series on
    the number of tilts, so it is a fallback, not an equivalent.
    """
    if os.path.isfile(sa_fullres_path):
        os.remove(sa_fullres_path)
        print(f"  Removed existing: {sa_fullres_path}")
    os.makedirs("tiltseries", exist_ok=True)
    hdr = EMData(lst_path, 0, True)
    apix = hdr["apix_x"]
    print(f"\n=== Step 3a — full-resolution SA tilt series for patchtracking ===")

    n_tot = EMUtil.get_image_count(lst_path)
    chunks = _read_sa_chunks(basename, avgseq, n_sel_neg, n_tot)

    if chunks is not None:
        with open(lst_path) as f:
            rows = [ln for ln in f if not ln.startswith("#")]
        idx = 0
        for a, b in chunks:
            acc = None
            for row in rows[a:b + 1]:
                fld = row.split()
                im = EMData(fld[1], int(fld[0]))
                if acc is None:
                    acc = im
                else:
                    acc.add(im)
            if acc is None:
                print(f"  WARNING: empty chunk at rows {a}-{b}")
                continue
            acc.mult(1.0 / float(b - a + 1))
            acc["apix_x"] = acc["apix_y"] = acc["apix_z"] = apix
            acc.write_image(f"{sa_fullres_path}:{compressbits}", idx)
            idx += 1
        sz = [c[1] - c[0] + 1 for c in chunks]
        n_neg_ch = sum(1 for c in chunks if c[0] < n_sel_neg)
        print(f"  Subaveraged from the recorded boundaries: {idx} tilts "
              f"({n_neg_ch} neg + {idx - n_neg_ch} pos), chunk sizes "
              f"{min(sz)}-{max(sz)} (target {avgseq})")
        return

    print("  no usable chunk record; falling back to fixed-stride chunking.")
    print("  WARNING: this can produce a different tilt count from the shrunk "
          "SA series. Verify before trusting the patchtracked result.")
    if n_sel_neg <= 0 or n_sel_neg >= n_tot:
        print("  (no neg/pos split available; single pass over the combined list)")
        cmd = (f"e2proc2d.py {lst_path} {sa_fullres_path}:{compressbits}"
               f" --avgseq={avgseq} --apix={apix:.4f}")
        print(cmd)
        ret = os.system(cmd)
        if ret != 0:
            print(f"Error: e2proc2d.py (fullres SA) exited with code {ret}")
            sys.exit(ret)
        return

    with open(lst_path) as _f:
        rows = [ln for ln in _f if not ln.startswith("#")]
    tmp_n, tmp_p = sa_fullres_path + ".negtmp.hdf", sa_fullres_path + ".postmp.hdf"
    lst_n, lst_p = sa_fullres_path + ".neg.lst", sa_fullres_path + ".pos.lst"
    with open(lst_n, "w") as _f:
        _f.write("#LST\n"); _f.writelines(rows[:n_sel_neg])
    with open(lst_p, "w") as _f:
        _f.write("#LST\n"); _f.writelines(rows[n_sel_neg:])
    for _l, _o in ((lst_n, tmp_n), (lst_p, tmp_p)):
        if os.path.isfile(_o):
            os.remove(_o)
        cmd = (f"e2proc2d.py {_l} {_o}:{compressbits}"
               f" --avgseq={avgseq} --apix={apix:.4f}")
        print(cmd)
        if os.system(cmd) != 0:
            print(f"Error: e2proc2d.py (fullres SA) failed on {_l}")
            sys.exit(1)
    idx = 0
    for _o in (tmp_n, tmp_p):
        for _i in range(EMUtil.get_image_count(_o)):
            EMData(_o, _i).write_image(f"{sa_fullres_path}:{compressbits}", idx)
            idx += 1
        os.remove(_o)
    for _l in (lst_n, lst_p):
        os.remove(_l)
    n_neg_ch = -(-n_sel_neg // avgseq)
    n_pos_ch = -(-(n_tot - n_sel_neg) // avgseq)
    print(f"  Subaveraged each half separately: {idx} tilts "
          f"({n_neg_ch} neg + {n_pos_ch} pos), no chunk crosses the boundary")


def make_fullframe_sh4(lst_path, fullsh4_path, compressbits, meanshrink=4):
    """Shrink the combined LST (no sub-averaging) → sets/."""
    if os.path.isfile(fullsh4_path):
        os.remove(fullsh4_path)
        print(f"  Removed existing: {fullsh4_path}")
    hdr = EMData(lst_path, 0, True)
    output_apix = hdr["apix_x"] * meanshrink
    cmd = f"e2proc2d.py {lst_path} {fullsh4_path}:{compressbits} --meanshrink={meanshrink} --apix={output_apix:.4f}"
    print(f"\n=== Step 4 — full-frame shrink ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Error: e2proc2d.py exited with code {ret}")
        sys.exit(ret)


def run_plot_tltparams(json_path, outdir="diagnostics"):
    """Best-effort tlt_params diagnostic plots; a failure here must not abort the pipeline."""
    cmd = f"e3plot_tltparams.py {json_path} --outdir {outdir}"
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e3plot_tltparams.py exited with code {ret} (non-fatal)")


def run_measure(basename, sa_tiltseries, outdir="diagnostics", binning=2, ntile=6):
    """Best-effort aligned-tilt stack + bead diagnostics; must not abort the pipeline.

    Writes diagnostics/FidStack_{measured,reproj,diff}_<basename>_SA.hdf. The
    measured stack is the aligned tilt series -- the thing an operator actually
    looks at to judge whether a run worked -- so it is produced automatically
    rather than left as a second command to type.

    Doing it here also removes two mistakes that are easy to make by hand:

    * the tomogram's name depends on whether patchtracking ran. A 1K input gives
      tomograms/<base>_SA.hdf; a 4K input (patchtrack>0) gives
      tomograms/<base>__bin4_SA.hdf. Asking for the wrong one fails outright.
    * the ttparams and the tilt series must be the pair the alignment actually
      produced. Mixing a fresh alignment with a stale series gives a plausible
      stack describing neither.

    Non-fatal by design: on RatNeuron ct3 the reconstruction succeeded and only
    the measurement failed (59 parameter rows against 58 tilt images). Aborting
    there would have discarded a good 34-minute run.

    sa_tiltseries must be the SHRUNK series (<base>__SA<n>_Sh<m>.hdf), NOT the
    full-resolution one e2tomogram aligns when patchtracking. Either would give a
    correct answer -- e3fid_reproject rescales the parameters from apix_unbin --
    but the full-res version produces a needlessly large stack and takes far
    longer for a picture of the same thing.
    """
    bin4 = os.path.join("tomograms", f"{basename}__bin4_SA.hdf")
    plain = os.path.join("tomograms", f"{basename}_SA.hdf")
    tomo = bin4 if os.path.isfile(bin4) else plain
    if not os.path.isfile(tomo):
        print(f"Warning: no SA tomogram at {bin4} or {plain}; skipping measurement")
        return
    tpm = f"{basename}_ttparams.npy"
    pks = f"{basename}_pks.npy"
    for f in (tpm, pks, sa_tiltseries):
        if not os.path.isfile(f):
            print(f"Warning: {f} missing; skipping measurement")
            return
    cmd = (f"e3fid_reproject.py --tomo {tomo} --tilts {sa_tiltseries}"
           f" --ttparams {tpm} --pks {pks}"
           f" --bin {binning} --ntile {ntile} --outdir {outdir}"
           f" --label {basename}_SA")
    print(f"\n=== Step 1b — aligned tilt stack + diagnostics ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e3fid_reproject.py exited with code {ret} (non-fatal)")
    else:
        print(f"  wrote {outdir}/FidStack_measured_{basename}_SA.hdf")


def run_interpolate_eulers(tomo_json, fullsh4_path, interp_json, interp_log):
    """Interpolate SA-level Euler angles to every full frame."""
    cmd = (f"e2InterpolateJsonEulersAmongFramesV2.py"
           f" --JsonIn {tomo_json}"
           f" --finalstackName {fullsh4_path}"
           f" --JsonOut {interp_json}")
    print(f"\n=== Step 5 — interpolate Euler angles to full frames ===")
    print(cmd + f" > {interp_log}")
    ret = _run_logged(cmd, interp_log)
    if ret != 0:
        print(f"Error: e2InterpolateJsonEulersAmongFramesV2.py exited with code {ret}")
        sys.exit(ret)


def run_stage2(interp_json, label, log_path, filterres=40.0):
    """Run e2tomogram_stage2.py → tomograms/<label>.hdf."""
    cmd = (f"e2tomogram_stage2.py"
           f" {interp_json} {label}"
           f" --filterres={filterres}")
    print(f"\n=== Step 6 — stage-2 reconstruction ===")
    print(cmd + f" > {log_path}")
    ret = _run_logged(cmd, log_path)
    if ret != 0:
        print(f"Error: e2tomogram_stage2.py exited with code {ret}")
        sys.exit(ret)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    progname = os.path.basename(sys.argv[0])
    usage = """e3ctomo_reconstruct.py [options]

Four-step continuous-tilt reconstruction pipeline:
  1. e2tomogram.py --stage2prep  (reconstruct SA tilt series + write rich JSON)
  2. Full-frame shrink (no SA)  ->  sets/
  3. Interpolate SA Euler angles to all full frames
  4. e2tomogram_stage2.py  ->  tomograms/
Steps 1 and 2 of the original pipeline (build LST, sub-average) are now
performed automatically by e3pick_ranges.py after range selection.
"""
    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    parser.add_argument("--basename",     type=str,  default="",   help="Short label for output files (e.g. CT07). Auto-loaded from Pick Ranges.", guitype='strbox', row=0, col=0, rowspan=1, colspan=3)
    parser.add_argument("--avgseq_sa",     type=int,  default=60,        help="Frames averaged per tilt in SA series (must match Pick Ranges).", guitype='intbox', row=1, col=0, rowspan=1, colspan=1)
    parser.add_argument("--meanshrink_sa", type=int,  default=4,        help="Spatial shrink used when building SA tiltseries in Pick Ranges (4K→1K: 4; 2K→1K: 2).", guitype='intbox', row=1, col=1, rowspan=1, colspan=1)
    parser.add_argument("--negStartAngle", type=float, default=0.0, help="Nominal start angle (deg) of the FULL negative-tilt acquisition sweep (before Pick Ranges selection). Required, along with negEndAngle/posStartAngle/posEndAngle, to compute zeroid and tltstep.", guitype='floatbox', row=2, col=0, rowspan=1, colspan=1)
    parser.add_argument("--negEndAngle",   type=float, default=0.0, help="Nominal end angle (deg) of the FULL negative-tilt acquisition sweep.", guitype='floatbox', row=2, col=1, rowspan=1, colspan=1)
    parser.add_argument("--posStartAngle", type=float, default=0.0, help="Nominal start angle (deg) of the FULL positive-tilt acquisition sweep.", guitype='floatbox', row=2, col=2, rowspan=1, colspan=1)
    parser.add_argument("--posEndAngle",   type=float, default=0.0, help="Nominal end angle (deg) of the FULL positive-tilt acquisition sweep.", guitype='floatbox', row=3, col=0, rowspan=1, colspan=1)
    parser.add_argument("--fixytilt", default=False, action="store_true", help="Hold each tilt's tilt_y angle at the computed rawtlt value instead of refining it. Recommended for continuous tilt: neighbouring tilts are a fraction of a degree apart, so per-tilt tilt_y refinement over-fits noise and scrambles the known stage angles, which the frame-level interpolation then amplifies into large angular sampling-density swings.", guitype='boolbox', row=3, col=1, rowspan=1, colspan=1)
    parser.add_argument("--compressbits", type=int,  default=6,    help="Compression bits for output.", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
    parser.add_argument("--meanshrink",   type=int,  default=4,    help="Spatial shrink for full-frame stack (2 for 2K cameras, 4 for 4K cameras).", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
    parser.add_argument("--threads",      type=int,  default=1,    help="CPU threads for e2tomogram.py --stage2prep.", guitype='intbox', row=4, col=2, rowspan=1, colspan=1)
    parser.add_argument("--clipz",        type=int,  default=768,  help="Z clip size (pixels).", guitype='intbox', row=5, col=0, rowspan=1, colspan=1)
    parser.add_argument("--tltkeep",      type=float, default=1.0, help="Fraction of tilt images to keep (1.0 = keep all).", guitype='floatbox', row=5, col=1, rowspan=1, colspan=1)
    parser.add_argument("--patchtrack",   type=int,  default=-1,   help="Patch tracking iterations for e2tomogram.py (-1=off, 2=recommended when used).", guitype='intbox', row=5, col=2, rowspan=1, colspan=1)
    parser.add_argument("--extrapad",   default=False,  action="store_true", help="Pass --extrapad to e2tomogram.py --stage2prep (recommended, but roughly triples the per-tile memory footprint of the bin4-tiling reconstruction step; uncheck to cut memory use).", guitype='boolbox', row=6, col=0, rowspan=1, colspan=1)
    parser.add_argument("--moretile",   default=False, action="store_true", help="Pass --moretile to e2tomogram.py --stage2prep (more tiles, slower).", guitype='boolbox', row=6, col=1, rowspan=1, colspan=1)
    parser.add_argument("--keeptmp",    default=False, action="store_true", help="Keep tomorecon/ temp directory (contains aligned tilt series for inspection).", guitype='boolbox', row=6, col=2, rowspan=1, colspan=1)
    parser.add_argument("--skip_tomogram",  default=False, action="store_true", help="Skip step 1: e2tomogram.py --stage2prep.", guitype='boolbox', row=7, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_fullsh4",   default=False, action="store_true", help="Skip step 2: full-frame shrink stack.", guitype='boolbox', row=7, col=1, rowspan=1, colspan=1)
    parser.add_argument("--skip_interp",    default=False, action="store_true", help="Skip step 3: Euler interpolation.", guitype='boolbox', row=7, col=2, rowspan=1, colspan=1)
    parser.add_argument("--skip_patchsa",   default=False, action="store_true", help="Skip Step 3a: rebuilding the full-res SA series for patchtracking. Safe to check once tiltseries/<basename>__SA<avgseq_sa>.hdf already exists and you're only iterating on later params (threads/extrapad/clipz/bxsz/etc.), which don't affect this file.", guitype='boolbox', row=9, col=0, rowspan=1, colspan=1)
    parser.add_argument("--niter",        type=str,  default="2,1,1,1", help="Refinement iterations for e2tomogram.py. '2,1,1,1' is the landmark-based default. Set to '0' to disable landmark refinement and rely on patch tracking alone (with --patchtrack=2): recommended for platelets and other thick cellular specimens, which lack the strong localized features landmark alignment needs.", guitype='strbox', row=9, col=1, rowspan=1, colspan=1)
    parser.add_argument("--tltax",        type=float, default=0.5,  help="Initial tilt axis angle in degrees (default 0.5). Check against the images before changing; a wrong value rotates the whole geometry.", guitype='floatbox', row=9, col=2, rowspan=1, colspan=1)
    parser.add_argument("--npk",          type=int,  default=20,   help="Number of landmarks/peaks for e2tomogram.py (default 20).", guitype='intbox', row=10, col=0, rowspan=1, colspan=1)
    parser.add_argument("--correctrot",   default=False, action="store_true", help="Pass --correctrot to e2tomogram.py (rotate the tomogram so the tilt axis is aligned).", guitype='boolbox', row=10, col=1, rowspan=1, colspan=1)
    parser.add_argument("--posz",         default=False, action="store_true", help="Pass --posz to e2tomogram.py (auto positioning along z). Recommended for platelets.", guitype='boolbox', row=10, col=2, rowspan=1, colspan=1)
    parser.add_argument("--xdrift",       default=False, action="store_true", help="Pass --xdrift to e2tomogram.py (extra correction for drift along x). Recommended for platelets.", guitype='boolbox', row=11, col=0, rowspan=1, colspan=1)
    parser.add_argument("--autoclipxy",   default=False, action="store_true", help="Pass --autoclipxy to e2tomogram.py (clip x/y to the useful region).", guitype='boolbox', row=11, col=1, rowspan=1, colspan=1)
    parser.add_argument("--bxsz",         type=int,  default=32,   help="Box size for feature tracking in e2tomogram.py (default 32; increase for patchtracking on large fields).", guitype='intbox', row=8, col=1, rowspan=1, colspan=1)
    parser.add_argument("--skip_measure", default=False, action="store_true", help="Skip the aligned-tilt stack and bead diagnostics (step 1b).", guitype='boolbox', row=12, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_stage2",    default=False, action="store_true", help="Skip step 4: stage-2 reconstruction.", guitype='boolbox', row=8, col=0, rowspan=1, colspan=1)
    parser.add_argument("--filterres",    type=float, default=40.0, help="Target lowpass resolution (Angstroms) for both the SA and stage-2 full reconstructions. Each stage recomputes its own Nyquist cutoff from this against its own volume's apix, so the two are directly comparable at equal resolution.", guitype='floatbox', row=8, col=2, rowspan=1, colspan=1)

    parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross-platform PPID.")

    (options, args) = parser.parse_args()

    logid = E2init(sys.argv, options.ppid)

    # Auto-load basename from pick_ranges.json if not supplied
    if not options.basename:
        _pr = os.path.join("info", "pick_ranges.json")
        if os.path.isfile(_pr):
            try:
                with open(_pr) as _f:
                    options.basename = json.load(_f).get("basename", "")
            except Exception:
                pass

    if not options.basename:
        print("Error: --basename is required (or run Pick Ranges first).")
        E2end(logid)
        return

    # Compute all output paths upfront from basename + avgseq
    bn           = options.basename
    lst_path     = os.path.join("sets", f"{bn}_combined.lst")
    _sa_sh_tag   = f"_Sh{options.meanshrink_sa}" if options.meanshrink_sa > 1 else ""
    sa_path      = os.path.join("tiltseries", f"{bn}__SA{options.avgseq_sa}{_sa_sh_tag}.hdf")
    # Full-res SA (meanshrink=1) used as input to e2tomogram.py when patchtrack > 0
    sa_fullres_path = os.path.join("tiltseries", f"{bn}__SA{options.avgseq_sa}.hdf")
    # e2tomogram.py reads tomo_json from the path it was given, so derive from whichever SA we pass
    tomo_sa_path = sa_fullres_path if options.patchtrack > 0 else sa_path
    fullsh4_path = _fullsh4_path(lst_path, bn, options.meanshrink)
    tomo_json    = _tomogram_json(tomo_sa_path)
    interp_json  = _interp_json(fullsh4_path, bn, options.meanshrink)
    interp_log   = _interp_log(fullsh4_path, bn, options.meanshrink)
    s2_label     = _stage2_label(bn, interp_json)
    s2_log       = os.path.join("tomograms", f"{s2_label}_stage2.txt")

    print(f"\nOutput paths:")
    if options.patchtrack > 0:
        print(f"  SA tiltseries (full-res for patchtrack): {sa_fullres_path}")
    print(f"  SA tiltseries: {sa_path}")
    print(f"  Tomogram JSON: {tomo_json}")
    print(f"  Full Sh4     : {fullsh4_path}")
    print(f"  Interp JSON  : {interp_json}")
    print(f"  Stage-2 tomo : tomograms/{s2_label}.hdf")

    # Step 1
    if options.skip_tomogram:
        print(f"\n--- Step 1 skipped (--skip_tomogram)  [{tomo_sa_path}]")
    else:
        if (options.negStartAngle == 0.0 and options.negEndAngle == 0.0 and
                options.posStartAngle == 0.0 and options.posEndAngle == 0.0):
            print("Error: --negStartAngle/--negEndAngle/--posStartAngle/--posEndAngle are "
                  "required (the nominal full-scan angles for each half, from your "
                  "acquisition scheme) -- needed to compute zeroid and tltstep.")
            E2end(logid)
            return
        if options.negStartAngle == options.negEndAngle:
            print("Error: --negStartAngle and --negEndAngle must differ.")
            E2end(logid)
            return
        if options.posStartAngle == options.posEndAngle:
            print("Error: --posStartAngle and --posEndAngle must differ.")
            E2end(logid)
            return

        # Frame ranges must come from THIS dataset's own pick_ranges file. A project can
        # hold several datasets, and info/pick_ranges.json only reflects whichever was
        # picked last -- reading it blindly builds the rawtlt from another dataset's
        # ranges, which fails silently and produces a plausible-looking wrong tomogram.
        _prb = os.path.join("info", f"pick_ranges_{bn}.json")
        _pr_used = _prb if os.path.isfile(_prb) else os.path.join("info", "pick_ranges.json")
        try:
            with open(_pr_used) as _f:
                _pr_data = json.load(_f)
        except Exception as e:
            print(f"Error: could not read {_pr_used} ({e}) -- run Pick Ranges first.")
            E2end(logid)
            return
        _pr_bn = _pr_data.get("basename")
        if _pr_bn and _pr_bn != bn:
            print(f"Error: {_pr_used} holds frame ranges for '{_pr_bn}', not '{bn}'. "
                  f"Re-run Pick Ranges for {bn}, or restore info/pick_ranges_{bn}.json.")
            E2end(logid)
            return
        neg_start = _pr_data.get("neg_start")
        neg_end   = _pr_data.get("neg_end")
        pos_start = _pr_data.get("pos_start")
        pos_end   = _pr_data.get("pos_end")
        if None in (neg_start, neg_end, pos_start, pos_end):
            print(f"Error: missing neg_start/neg_end/pos_start/pos_end in {_pr_used} -- "
                  f"run Pick Ranges first.")
            E2end(logid)
            return
        print(f"  frame ranges from {_pr_used}: neg [{neg_start},{neg_end}] pos [{pos_start},{pos_end}]")
        # Number of negative-half rows in the combined LST: the split point the
        # full-res SA series must respect so no chunk straddles the acquisition
        # boundary (see SeamTilt_20260823.tex).
        _n_sel_neg = neg_end - neg_start + 1

        if options.patchtrack > 0:
            if options.skip_patchsa:
                print(f"\n--- Step 3a skipped (--skip_patchsa)  [{sa_fullres_path}]")
            else:
                print(f"\n--- Step 3a: building full-res SA series for patchtracking ---")
                make_fullres_sa(lst_path, sa_fullres_path, options.avgseq_sa,
                                options.compressbits, n_sel_neg=_n_sel_neg,
                                basename=bn)

        sa_angles, zeroid = _compute_sa_tilt_angles(
            tomo_sa_path, neg_start, neg_end, pos_start, pos_end,
            options.avgseq_sa, options.negStartAngle, options.negEndAngle,
            options.posStartAngle, options.posEndAngle)
        rawtlt_path = _write_rawtlt(sa_angles, tomo_sa_path)

        run_tomogram(tomo_sa_path, zeroid, rawtlt_path, options.threads, options.clipz, options.tltkeep,
                     extrapad=options.extrapad, moretile=options.moretile, keeptmp=options.keeptmp,
                     niter=options.niter, tltax=options.tltax, npk=options.npk,
                     compressbits=options.compressbits, correctrot=options.correctrot,
                     posz=options.posz, xdrift=options.xdrift, autoclipxy=options.autoclipxy,
                     patchtrack=options.patchtrack, bxsz=options.bxsz, filterres=options.filterres,
                     fixytilt=options.fixytilt)
        _copy_tlt_keys_to_sa_info(tomo_json, bn)
        run_plot_tltparams(tomo_json)
        if options.skip_measure:
            print("\n--- Step 1b skipped (--skip_measure)")
        else:
            run_measure(bn, sa_path)

    # Step 2
    if options.skip_fullsh4:
        print(f"\n--- Step 2 skipped (--skip_fullsh4)  [{fullsh4_path}]")
    else:
        make_fullframe_sh4(lst_path, fullsh4_path, options.compressbits, options.meanshrink)

    # Step 3
    if options.skip_interp:
        print(f"\n--- Step 3 skipped (--skip_interp)  [{interp_json}]")
    else:
        run_interpolate_eulers(tomo_json, fullsh4_path, interp_json, interp_log)
        run_plot_tltparams(interp_json)

    # Step 4
    if options.skip_stage2:
        print(f"\n--- Step 4 skipped (--skip_stage2)  [tomograms/{s2_label}.hdf]")
    else:
        run_stage2(interp_json, s2_label, s2_log, filterres=options.filterres)

    E2end(logid)

if __name__ == "__main__":
    main()
