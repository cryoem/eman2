#!/usr/bin/env python
# e3ctomo_reconstruct.py  --  Continuous-tilt reconstruction launcher.
from past.utils import old_div
from builtins import range
from EMAN2 import *
import sys, os, json


# ---------------------------------------------------------------------------
# Path helpers  (all accept optional basename; fall back to auto-derived names)
# ---------------------------------------------------------------------------

def _lst_path(movie_neg, movie_pos, basename=None):
    if basename:
        return os.path.join("sets", f"{basename}_combined.lst")
    neg_base = os.path.splitext(os.path.basename(movie_neg))[0]
    pos_base = os.path.splitext(os.path.basename(movie_pos))[0]
    return os.path.join("sets", f"{neg_base}__{pos_base}__combined.lst")

def _sa_path(lst_path, avgseq, basename=None):
    if basename:
        # double-underscore before SA tag so EMAN2 base_name() returns just basename
        return os.path.join("tiltseries", f"{basename}__SA{avgseq}_Sh4.hdf")
    stem = os.path.splitext(os.path.basename(lst_path))[0]
    return os.path.join("tiltseries", f"{stem}__SA{avgseq}_Sh4.hdf")

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

def make_combined_lst(movie_neg, movie_pos, neg_start, neg_end, pos_start, pos_end, lst_path):
    """
    Write the combined LST.

    Combined structure (mirrors the thumbnail movie):
      positions 0 .. N_neg-1  : neg frames, reversed
                                 combined index i  ->  original neg frame  N_neg-1-i
      positions N_neg .. end   : pos frames, normal order
                                 combined index i  ->  pos frame  i - N_neg

    neg_start, neg_end  : combined indices bounding the neg section to keep
    pos_start, pos_end  : combined indices bounding the pos section to keep
                          (pos_start >= N_neg by construction from pick_ranges)
    """
    n_neg = EMUtil.get_image_count(movie_neg)

    os.makedirs("sets", exist_ok=True)

    neg_start = max(0, neg_start)
    neg_end   = min(n_neg - 1, neg_end)
    pos_end_max = n_neg + EMUtil.get_image_count(movie_pos) - 1
    pos_start = max(n_neg, pos_start)
    pos_end   = min(pos_end_max, pos_end)

    with open(lst_path, "w") as f:
        f.write("#LST\n")
        for i in range(neg_start, neg_end + 1):
            f.write(f"{n_neg - 1 - i}\t{movie_neg}\n")
        for i in range(pos_start, pos_end + 1):
            f.write(f"{i - n_neg}\t{movie_pos}\n")

    n_neg_written = neg_end - neg_start + 1
    n_pos_written = pos_end - pos_start + 1
    print(f"Wrote combined LST: {lst_path}")
    print(f"  Neg frames kept : {n_neg_written}  "
          f"(combined {neg_start}–{neg_end}  →  orig neg {n_neg-1-neg_start}–{n_neg-1-neg_end})")
    print(f"  Pos frames kept : {n_pos_written}  "
          f"(pos movie {pos_start-n_neg}–{pos_end-n_neg})")


def make_tiltseries_sa(lst_path, sa_path, avgseq, compressbits):
    """Sub-average and shrink the combined LST → tiltseries/."""
    os.makedirs("tiltseries", exist_ok=True)
    cmd = (f"e2proc2d.py {lst_path} {sa_path}:{compressbits} "
           f"--avgseq={avgseq} --meanshrink=4")
    print(f"\n=== Step 2 — sub-average + shrink ===")
    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print(f"Warning: e2proc2d.py exited with code {ret}")


def run_tomogram(sa_path, tilt_range, threads, clipz, tltkeep):
    """Run /home3/Python/e2tomogram.py --stage2prep: reconstruct and write rich JSON."""
    n_sa = EMUtil.get_image_count(sa_path)
    tltstep = tilt_range / n_sa
    cmd = (
        f"/home3/Python/e2tomogram.py {sa_path}"
        f" --stage2prep"
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

Six-step continuous-tilt reconstruction pipeline:
  1. Build combined LST  (neg reversed + pos, trimmed to picked ranges)
  2. Sub-average + shrink  ->  tiltseries/
  3. e2tomogram.py --stage2prep  (reconstruct + write rich JSON)
  4. Full-frame shrink (no SA)  ->  sets/
  5. Interpolate SA Euler angles to all full frames
  6. e2tomogram_stage2_Oct7th2025v2.py  ->  tomograms/
Frame ranges are auto-loaded from info/pick_ranges.json if present.
"""
    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    parser.add_argument("--basename",     type=str,  default="",    help="Short label for output files, e.g. CT07 or 08391_08393. If omitted, names are derived from the movie filenames.", guitype='strbox', row=0, col=0, rowspan=1, colspan=3)
    parser.add_argument("--movie_pos",    type=str,  default="",    help="Positive tilt movie (gain-corrected).", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=1, col=0, rowspan=1, colspan=3)
    parser.add_argument("--movie_neg",    type=str,  default="",    help="Negative tilt movie (gain-corrected).", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=2, col=0, rowspan=1, colspan=3)
    parser.add_argument("--gainfile",     type=str,  default="GainFiducial.hdf", help="Gain reference HDF file.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=3, col=0, rowspan=1, colspan=3)
    parser.add_argument("--tiltRange",    type=int,  default=120,   help="Total tilt range in degrees.", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
    parser.add_argument("--compressbits", type=int,  default=6,     help="Compression bits for output.", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
    parser.add_argument("--startTltNeg",  type=int,  default=0,     help="Start frame — negative tilt.", guitype='intbox', row=5, col=0, rowspan=1, colspan=1)
    parser.add_argument("--endTltNeg",    type=int,  default=9999,  help="End frame — negative tilt.", guitype='intbox', row=5, col=1, rowspan=1, colspan=1)
    parser.add_argument("--startTltPos",  type=int,  default=0,     help="Start frame — positive tilt.", guitype='intbox', row=6, col=0, rowspan=1, colspan=1)
    parser.add_argument("--endTltPos",    type=int,  default=9999,  help="End frame — positive tilt.", guitype='intbox', row=6, col=1, rowspan=1, colspan=1)
    parser.add_argument("--avgseq",       type=int,   default=60,   help="Frames to average per tilt for sub-averaged tilt series.", guitype='intbox', row=7, col=0, rowspan=1, colspan=1)
    parser.add_argument("--threads",      type=int,   default=12,   help="CPU threads for e2tomogram.py --stage2prep.", guitype='intbox', row=7, col=1, rowspan=1, colspan=1)
    parser.add_argument("--clipz",        type=int,   default=768,  help="Z clip size (pixels) for e2tomogram.py --stage2prep.", guitype='intbox', row=7, col=2, rowspan=1, colspan=1)
    parser.add_argument("--tltkeep",      type=float, default=1.0,  help="Fraction of tilt images to keep (1.0 = keep all).", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_lst",       default=False, action="store_true", help="Skip step 1: build combined LST.", guitype='boolbox', row=9, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_sa",        default=False, action="store_true", help="Skip step 2: sub-averaged tilt series.", guitype='boolbox', row=9, col=1, rowspan=1, colspan=1)
    parser.add_argument("--skip_tomogram",  default=False, action="store_true", help="Skip step 3: e2tomogram.py --stage2prep.", guitype='boolbox', row=9, col=2, rowspan=1, colspan=1)
    parser.add_argument("--skip_fullsh4",   default=False, action="store_true", help="Skip step 4: full-frame Sh4 stack.", guitype='boolbox', row=10, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_interp",    default=False, action="store_true", help="Skip step 5: Euler interpolation.", guitype='boolbox', row=10, col=1, rowspan=1, colspan=1)
    parser.add_argument("--skip_stage2",    default=False, action="store_true", help="Skip step 6: stage-2 reconstruction.", guitype='boolbox', row=10, col=2, rowspan=1, colspan=1)

    parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross-platform PPID.")

    (options, args) = parser.parse_args()

    logid = E2init(sys.argv, options.ppid)

    # Load frame ranges from pick_ranges.json if present
    json_path = os.path.join(os.getcwd(), "info", "pick_ranges.json")
    if os.path.isfile(json_path):
        try:
            with open(json_path) as f:
                d = json.load(f)
            if not d.get("unidirectional"):
                if options.startTltNeg == 0:    options.startTltNeg = d.get("neg_start", 0)
                if options.endTltNeg   == 9999: options.endTltNeg   = d.get("neg_end",   9999)
                if options.startTltPos == 0:    options.startTltPos = d.get("pos_start", 0)
                if options.endTltPos   == 9999: options.endTltPos   = d.get("pos_end",   9999)
            print(f"Loaded frame ranges from {json_path}")
        except Exception as e:
            print(f"Warning: could not read {json_path}: {e}")

    if not (options.movie_neg and options.movie_pos):
        print("Error: --movie_neg and --movie_pos are required.")
        E2end(logid)
        return

    # Compute all output paths upfront
    bn           = options.basename or None
    lst_path     = _lst_path(options.movie_neg, options.movie_pos, bn)
    sa_path      = _sa_path(lst_path, options.avgseq, bn)
    fullsh4_path = _fullsh4_path(lst_path, bn)
    tomo_json    = _tomogram_json(sa_path)
    interp_json  = _interp_json(fullsh4_path, bn)
    interp_log   = _interp_log(fullsh4_path, bn)
    s2_label     = _stage2_label(bn, interp_json)
    s2_log       = os.path.join("tomograms", f"{s2_label}_stage2.txt")

    print(f"\nOutput paths:")
    print(f"  LST          : {lst_path}")
    print(f"  SA tiltseries: {sa_path}")
    print(f"  Tomogram JSON: {tomo_json}")
    print(f"  Full Sh4     : {fullsh4_path}")
    print(f"  Interp JSON  : {interp_json}")
    print(f"  Stage-2 tomo : tomograms/{s2_label}.hdf")

    # Step 1
    if options.skip_lst:
        print(f"\n--- Step 1 skipped (--skip_lst)  [{lst_path}]")
    else:
        make_combined_lst(
            options.movie_neg, options.movie_pos,
            options.startTltNeg, options.endTltNeg,
            options.startTltPos, options.endTltPos,
            lst_path,
        )

    # Step 2
    if options.skip_sa:
        print(f"\n--- Step 2 skipped (--skip_sa)  [{sa_path}]")
    else:
        make_tiltseries_sa(lst_path, sa_path, options.avgseq, options.compressbits)

    # Step 3
    if options.skip_tomogram:
        print(f"\n--- Step 3 skipped (--skip_tomogram)  [{sa_path}]")
    else:
        run_tomogram(sa_path, options.tiltRange, options.threads, options.clipz, options.tltkeep)

    # Step 4
    if options.skip_fullsh4:
        print(f"\n--- Step 4 skipped (--skip_fullsh4)  [{fullsh4_path}]")
    else:
        make_fullframe_sh4(lst_path, fullsh4_path, options.compressbits)

    # Step 5
    if options.skip_interp:
        print(f"\n--- Step 5 skipped (--skip_interp)  [{interp_json}]")
    else:
        run_interpolate_eulers(tomo_json, fullsh4_path, interp_json, interp_log)

    # Step 6
    if options.skip_stage2:
        print(f"\n--- Step 6 skipped (--skip_stage2)  [tomograms/{s2_label}.hdf]")
    else:
        run_stage2(interp_json, s2_label, s2_log)

    E2end(logid)

if __name__ == "__main__":
    main()
