#!/usr/bin/env python
# e3ctomo_reconstruct.py  --  Continuous-tilt reconstruction launcher.
from past.utils import old_div
from builtins import range
from EMAN2 import *
import sys, os, json


# ---------------------------------------------------------------------------
# Path helpers  (all accept optional basename; fall back to auto-derived names)
# ---------------------------------------------------------------------------

def _fullsh4_path(lst_path, basename=None):
    if basename:
        return os.path.join("sets", f"{basename}_Full_Sh4.hdf")
    stem = os.path.splitext(os.path.basename(lst_path))[0]
    return os.path.join("sets", f"{stem}_Sh4.hdf")

def _tomogram_json(sa_path):
    """Path that e2tomogram.py --stage2prep writes: info/<sa_base>_recon_tomo_final.json."""
    stem = os.path.splitext(os.path.basename(sa_path))[0]
    first = stem.split("__")[0]
    return os.path.join("info", f"{first}_recon_tomo_final.json")

def _stage2_label(basename=None, interp_json=None):
    """Label passed to e2tomogram_stage2_Oct7th2025v2.py as baseNameForOutputs."""
    if basename:
        return basename
    if interp_json:
        stem = os.path.splitext(os.path.basename(interp_json))[0]
        if stem.endswith("_Full_Sh4"):
            stem = stem[:-len("_Full_Sh4")]
        return stem
    return "stage2"

def _interp_json(fullsh4_path, basename=None):
    if basename:
        return os.path.join("info", f"{basename}_Full_Sh4.json")
    stem = os.path.splitext(os.path.basename(fullsh4_path))[0]
    return os.path.join("info", f"{stem}.json")

def _interp_log(fullsh4_path, basename=None):
    if basename:
        return os.path.join("info", f"{basename}_Full_Sh4.txt")
    stem = os.path.splitext(os.path.basename(fullsh4_path))[0]
    return os.path.join("info", f"{stem}.txt")


# ---------------------------------------------------------------------------
# Pipeline steps  (all paths passed in; no internal path derivation)
# ---------------------------------------------------------------------------


def run_tomogram(sa_path, tilt_range, threads, clipz, tltkeep):
    """Run /home3/Python/e2tomogram.py --stage2prep: reconstruct and write rich JSON."""
    n_sa = EMUtil.get_image_count(sa_path)
    tltstep = tilt_range / n_sa
    cmd = (
        f"/home3/Python/e2tomogram.py {sa_path}"
        f" --stage2prep"
        f" --outsuffix=_SA"
        f" --tltax=0.5"
        f" --zeroid=-1"
        f" --tltstep={tltstep:.6f}"
        f" --npk=20"
        f" --tltkeep={tltkeep}"
        f" --outsize=1k"
        f" --niter=2,1,1,1"
        f" --bytile"
        f" --notmp"
        f" --pkkeep=0.9"
        f" --compressbits=6"
        f" --clipz={clipz}"
        f" --bxsz=32"
        f" --filterres=40.0"
        f" --extrapad"
        f" --moretile"
        f" --rmbeadthr=-1.0"
        f" --threads={threads}"
        f" --patchtrack=-1"
    )
    print(f"\n=== Step 3 — tomogram reconstruction + stage2prep ({n_sa} SA frames, tltstep={tltstep:.4f}°) ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e2tomogram.py exited with code {ret}")


def make_fullframe_sh4(lst_path, fullsh4_path, compressbits):
    """Shrink the combined LST (no sub-averaging) → sets/."""
    cmd = f"e2proc2d.py {lst_path} {fullsh4_path}:{compressbits} --meanshrink=4"
    print(f"\n=== Step 4 — full-frame shrink ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e2proc2d.py exited with code {ret}")


def run_interpolate_eulers(tomo_json, fullsh4_path, interp_json, interp_log):
    """Interpolate SA-level Euler angles to every full frame."""
    cmd = (f"e2InterpolateJsonEulersAmongFramesV2.py"
           f" --JsonIn {tomo_json}"
           f" --finalstackName {fullsh4_path}"
           f" --JsonOut {interp_json}"
           f" > {interp_log}")
    print(f"\n=== Step 5 — interpolate Euler angles to full frames ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e2InterpolateJsonEulersAmongFramesV2.py exited with code {ret}")


def run_stage2(interp_json, label, log_path):
    """Run e2tomogram_stage2_Oct7th2025v2.py → tomograms/<label>.hdf."""
    cmd = (f"e2tomogram_stage2_Oct7th2025v2.py"
           f" {interp_json} {label}"
           f" > {log_path}")
    print(f"\n=== Step 6 — stage-2 reconstruction ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e2tomogram_stage2_Oct7th2025v2.py exited with code {ret}")


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
  4. e2tomogram_stage2_Oct7th2025v2.py  ->  tomograms/
Steps 1 and 2 of the original pipeline (build LST, sub-average) are now
performed automatically by e3pick_ranges.py after range selection.
"""
    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    parser.add_argument("--basename",     type=str,  default=pick_ranges_data.get("basename",""),   help="Short label for output files (e.g. CT07). Auto-loaded from Pick Ranges.", guitype='strbox', row=0, col=0, rowspan=1, colspan=3)
    parser.add_argument("--avgseq_sa",     type=int,  default=60,   help="Frames averaged per tilt in SA series (must match Pick Ranges).", guitype='intbox', row=1, col=0, rowspan=1, colspan=1)
    parser.add_argument("--tiltRange",    type=int,  default=120,  help="Total tilt range in degrees.", guitype='intbox', row=1, col=1, rowspan=1, colspan=1)
    parser.add_argument("--compressbits", type=int,  default=6,    help="Compression bits for output.", guitype='intbox', row=1, col=2, rowspan=1, colspan=1)
    parser.add_argument("--threads",      type=int,  default=2,    help="CPU threads for e2tomogram.py --stage2prep.", guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
    parser.add_argument("--clipz",        type=int,  default=768,  help="Z clip size (pixels).", guitype='intbox', row=2, col=1, rowspan=1, colspan=1)
    parser.add_argument("--tltkeep",      type=float, default=1.0, help="Fraction of tilt images to keep (1.0 = keep all).", guitype='floatbox', row=2, col=2, rowspan=1, colspan=1)
    parser.add_argument("--skip_tomogram",  default=False, action="store_true", help="Skip step 1: e2tomogram.py --stage2prep.", guitype='boolbox', row=3, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_fullsh4",   default=False, action="store_true", help="Skip step 2: full-frame Sh4 stack.", guitype='boolbox', row=3, col=1, rowspan=1, colspan=1)
    parser.add_argument("--skip_interp",    default=False, action="store_true", help="Skip step 3: Euler interpolation.", guitype='boolbox', row=3, col=2, rowspan=1, colspan=1)
    parser.add_argument("--skip_stage2",    default=False, action="store_true", help="Skip step 4: stage-2 reconstruction.", guitype='boolbox', row=4, col=0, rowspan=1, colspan=1)

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
    sa_path      = os.path.join("tiltseries", f"{bn}__SA{options.avgseq_sa}_Sh4.hdf")
    fullsh4_path = _fullsh4_path(lst_path, bn)
    tomo_json    = _tomogram_json(sa_path)
    interp_json  = _interp_json(fullsh4_path, bn)
    interp_log   = _interp_log(fullsh4_path, bn)
    s2_label     = _stage2_label(bn, interp_json)
    s2_log       = os.path.join("tomograms", f"{s2_label}_stage2.txt")

    print(f"\nOutput paths:")
    print(f"  SA tiltseries: {sa_path}")
    print(f"  Tomogram JSON: {tomo_json}")
    print(f"  Full Sh4     : {fullsh4_path}")
    print(f"  Interp JSON  : {interp_json}")
    print(f"  Stage-2 tomo : tomograms/{s2_label}.hdf")

    # Step 1
    if options.skip_tomogram:
        print(f"\n--- Step 1 skipped (--skip_tomogram)  [{sa_path}]")
    else:
        run_tomogram(sa_path, options.tiltRange, options.threads, options.clipz, options.tltkeep)

    # Step 2
    if options.skip_fullsh4:
        print(f"\n--- Step 2 skipped (--skip_fullsh4)  [{fullsh4_path}]")
    else:
        make_fullframe_sh4(lst_path, fullsh4_path, options.compressbits)

    # Step 3
    if options.skip_interp:
        print(f"\n--- Step 3 skipped (--skip_interp)  [{interp_json}]")
    else:
        run_interpolate_eulers(tomo_json, fullsh4_path, interp_json, interp_log)

    # Step 4
    if options.skip_stage2:
        print(f"\n--- Step 4 skipped (--skip_stage2)  [tomograms/{s2_label}.hdf]")
    else:
        run_stage2(interp_json, s2_label, s2_log)

    E2end(logid)

if __name__ == "__main__":
    main()
