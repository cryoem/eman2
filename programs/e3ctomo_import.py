#!/usr/bin/env python3
"""
e3ctomo_import.py  --  Continuous-tilt movie import: copy, thumbnail, time-average,
                        gain-correct, and CCF-align.

Workflow
--------
  i)   Copy input movies into outdir/movies/ using e2proc2d.py
         (sets apix, meanshrink, compressbits, compresslevel=1)
  iii) Create heavily downsampled thumbnails in outdir/thumbnails/
         (subaveraged by ntave, shrunk to 512x512)
  iv)  For each movie: rolling time-average + gain correction  -> *_GainTimeAve18.hdf
       then CCF-align each frame to its time-average           -> *_GainShifted.hdf

Note: step ii (gain calculation) is handled separately by e3movie.py.
      The gain reference is supplied here via --gain.

Usage
-----
  e3ctomo_import.py <neg1> <pos1> [<neg2> <pos2> …] \\
      --gain <gainfile> \\
      [--apix=3.714] [--meanshrink=4] [--compressbits=4] \\
      [--outdir=.] [--ntave=18] [--chunk=50] \\
      [--skip-copy] [--skip-thumbs]

Movies are supplied as interleaved neg/pos pairs but each is processed
independently through the pipeline.
"""

import sys
import os
import argparse
import glob
import datetime
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ── EMAN2/3 double-patch fix ──────────────────────────────────────────────────
import EMAN2_cppwrap as _cppwrap
_c_emdata_init      = _cppwrap.EMData.__init__
_c_get_image_count  = _cppwrap.EMUtil.get_image_count
_c_read_image       = _cppwrap.EMData.read_image
_c_read_images      = _cppwrap.EMData.read_images
_c_write_image      = _cppwrap.EMData.write_image
_c_write_images     = _cppwrap.EMData.write_images

from EMAN3 import EMData, EMUtil, Region

EMData.__initc           = _c_emdata_init
EMUtil.get_image_count_c  = staticmethod(_c_get_image_count)
EMData.read_image_c       = _c_read_image
EMData.read_images_c      = staticmethod(_c_read_images)
EMData.write_image_c      = _c_write_image
EMData.write_images_c     = staticmethod(_c_write_images)
# ─────────────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ct_functions import do_time_average_chunked_rolling_gain, center_clip, hdf_nimg

VALID_MEANSHRINK = {1, 2, 4}


def ts():
    print(datetime.datetime.now().strftime("  [%Y-%m-%d %H:%M:%S]"))


# ── Step i: copy & spatially downsample movies ────────────────────────────────

def step_i_copy_movies(input_movies, outdir, apix, meanshrink, compressbits):
    movies_dir = os.path.join(outdir, "movies")
    os.makedirs(movies_dir, exist_ok=True)

    for src in input_movies:
        basename = os.path.basename(src)
        dst = os.path.join(movies_dir, basename)
        if os.path.isfile(dst):
            os.remove(dst)
            print(f"Removed existing movie: {dst}")
        cmd = (f"e2proc2d.py {src} {dst}:4 "
               f"--apix={apix} --meanshrink={meanshrink} "
               f"--compressbits={compressbits}")
        print(cmd)
        os.system(cmd)


# ── Step iii: thumbnails ──────────────────────────────────────────────────────

def step_iii_thumbnails(input_movies, outdir, meanshrink, aveseqInt=18):
    thumbs_dir = os.path.join(outdir, "thumbnails")
    os.makedirs(thumbs_dir, exist_ok=True)

    for src in input_movies:
        basename = os.path.basename(src)
        nimg = hdf_nimg(src)
        img0 = EMData(src, 0)
        nx = img0.get_xsize()
        shrink_factor = max(1, nx // 512)

        dst = os.path.join(thumbs_dir, basename[:-4] + "_thumb.hdf")
        if os.path.isfile(dst):
            os.remove(dst)
            print(f"Removed existing thumbnail: {dst}")

        print(f"Creating thumbnails: {src} -> {dst}  "
              f"(avgseq={aveseqInt}, shrink={shrink_factor})")

        out_idx = 0
        for start in range(0, nimg, aveseqInt):
            stop = min(start + aveseqInt, nimg)
            chunk = EMData.read_images(src, list(range(start, stop)))
            avg = chunk[0].copy()
            avg.mult(0.0)
            for fr in chunk:
                avg += fr
            avg.mult(1.0 / len(chunk))
            if shrink_factor > 1:
                avg.process_inplace("math.meanshrink", {"n": shrink_factor})
            avg.write_image(dst + ":6", out_idx)
            out_idx += 1

        print(f"  wrote {out_idx} thumbnail frames")


# ── Step iii (post): combine neg+pos thumbnail pairs ─────────────────────────

def step_iii_combine_thumbnails(movie_files, outdir):
    """For each neg/pos pair combine into one file: neg reversed + pos appended."""
    thumbs_dir = os.path.join(outdir, "thumbnails")

    pairs = [(movie_files[i], movie_files[i + 1])
             for i in range(0, len(movie_files) - 1, 2)]

    for neg_movie, pos_movie in pairs:
        neg_base = os.path.basename(neg_movie)[:-4]
        pos_base = os.path.basename(pos_movie)[:-4]
        neg_thumb = os.path.join(thumbs_dir, neg_base + "_thumb.hdf")
        pos_thumb = os.path.join(thumbs_dir, pos_base + "_thumb.hdf")
        combined  = os.path.join(thumbs_dir, neg_base + "_combined_thumb.hdf")

        if not os.path.isfile(neg_thumb):
            print(f"WARNING: neg thumbnail not found, skipping pair: {neg_thumb}")
            continue
        if not os.path.isfile(pos_thumb):
            print(f"WARNING: pos thumbnail not found, skipping pair: {pos_thumb}")
            continue

        if os.path.isfile(combined):
            os.remove(combined)
            print(f"  Removed existing: {os.path.basename(combined)}")

        n_neg = hdf_nimg(neg_thumb)
        n_pos = hdf_nimg(pos_thumb)
        print(f"Combining thumbnails: {neg_base} ({n_neg} frames, reversed)"
              f" + {pos_base} ({n_pos} frames)")
        print(f"  -> {os.path.basename(combined)}")

        out_idx = 0
        for i in range(n_neg - 1, -1, -1):
            EMData(neg_thumb, i).write_image(combined + ":6", out_idx)
            out_idx += 1
        for i in range(n_pos):
            EMData(pos_thumb, i).write_image(combined + ":6", out_idx)
            out_idx += 1

        print(f"  wrote {out_idx} frames total")


# ── Step iv: rolling time-average + gain correction + CCF alignment ───────────

def step_iv_gain_and_align(input_movies, gain_file, apix,
                           number_to_be_averaged=18, chunkframes=50):
    gain_data = EMData(gain_file)

    for jMovie, hdf_file in enumerate(input_movies):
        print(f"\n=== Movie {jMovie}: {hdf_file} ===")
        ts()

        base = os.path.basename(hdf_file)
        gain_time_ave_file = base[:-4] + f"_GainTimeAve{number_to_be_averaged}.hdf"
        gain_shifted_file  = base[:-4] + "_GainShifted.hdf"
        nimg = hdf_nimg(hdf_file)
        print(f"  {nimg} frames")
        print(f"  Output GainTimeAve : {gain_time_ave_file}")
        print(f"  Output GainShifted : {gain_shifted_file}")

        for f in (gain_time_ave_file, gain_shifted_file):
            if os.path.isfile(f):
                os.remove(f)
                print(f"  Removed existing: {f}")

        # 4a: rolling time-average + gain correction
        do_time_average_chunked_rolling_gain(
            hdf_file,
            number_to_be_averaged=number_to_be_averaged,
            chunk=chunkframes,
            out=gain_time_ave_file,
            verbose=True,
            gain_fsp=gain_file,
        )
        ts()

        # 4b: CCF align each frame to its time-average
        movie_string = base[9:14] if len(base) > 14 else base[:5]

        os.makedirs("peaks", exist_ok=True)
        fn_shifts       = f"peaks/Shift_XY_{base[:-4]}.npy"
        fn_shifts_after = f"peaks/Shift_again_XY_{base[:-4]}.npy"

        peak_coordsx_Vec         = np.zeros(nimg)
        peak_coordsy_Vec         = np.zeros(nimg)
        peak_coordsx_Vec_Shifted = np.zeros(nimg)
        peak_coordsy_Vec_Shifted = np.zeros(nimg)

        # Gaussian low-pass filter for CCF (computed once per movie)
        rsz = 32
        filt = EMData(rsz * 2, rsz * 2, 1)
        filt.to_one()
        filt.process_inplace("mask.gaussian", {"outer_radius": 2})
        filt.process_inplace("xform.phaseorigin.tocorner")
        filttf = filt.do_fft()

        img0 = EMData(hdf_file, 0)
        clip_size = min(int(0.75 * img0.get_xsize()), 3072 // 2)
        start_x = (img0.get_xsize() - clip_size) // 2
        start_y = (img0.get_ysize() - clip_size) // 2
        region_xy = Region(start_x, start_y, clip_size, clip_size)

        try:
            os.remove(gain_shifted_file)
        except FileNotFoundError:
            pass

        for start in range(0, nimg, chunkframes):
            stop = min(start + chunkframes, nimg)

            if start % 200 == 0:
                print(f"  [align chunk] {start}:{stop}")

            now_chunk  = EMData.read_images(hdf_file,           range(start, stop))
            tave_chunk = EMData.read_images(gain_time_ave_file,  range(start, stop))

            shifted_chunk = []

            for j_local, (img_now, img_tave) in enumerate(zip(now_chunk, tave_chunk)):
                j_global = start + j_local

                # Gain-correct the raw frame
                img_now.process_inplace(
                    "math.fixgain.counting",
                    {"gain": gain_data, "gainmin": 2, "gainmax": 2},
                )

                # Pre-shift CCF
                img_clip  = img_now.get_clip(region_xy)
                img_clip.process_inplace("normalize.edgemean")
                tave_clip = img_tave.get_clip(region_xy)
                tave_clip.process_inplace("normalize.edgemean")

                ccf  = tave_clip.calc_ccf(img_clip)
                ccfB = ccf.process("xform.phaseorigin.tocorner")
                cccfBc = center_clip(ccfB, (2 * rsz, 2 * rsz))
                ccfsf  = cccfBc.do_fft()
                ccfsreal = (ccfsf * filttf).do_ift()

                peak = ccfsreal.calc_max_location()
                dx = int(peak[0] - rsz)
                dy = int(peak[1] - rsz)

                peak_coordsx_Vec[j_global] = dx
                peak_coordsy_Vec[j_global] = dy

                # Shift
                img_now.translate(dx, dy, 0)
                img_now = img_now + 0.001
                img_now.process_inplace("mask.dampedzeroedgefill",
                                        {"mode": "edgetapermean"})
                img_now = img_now - 0.001
                img_now["apix_x"] = apix
                img_now["apix_y"] = apix
                img_now["apix_z"] = apix
                shifted_chunk.append(img_now)

                # Post-shift CCF check
                img_clip2  = img_now.get_clip(region_xy)
                img_clip2.process_inplace("normalize.edgemean")
                tave_clip2 = img_tave.get_clip(region_xy)
                tave_clip2.process_inplace("normalize.edgemean")

                ccf2 = tave_clip2.calc_ccf(img_clip2)
                ccfB2 = ccf2.process("xform.phaseorigin.tocorner")
                cccfBc2 = center_clip(ccfB2, (2 * rsz, 2 * rsz))
                ccfsf2  = cccfBc2.do_fft()
                ccfsreal2 = (ccfsf2 * filttf).do_ift()

                peak2 = ccfsreal2.calc_max_location()
                peak_coordsx_Vec_Shifted[j_global] = peak2[0] - rsz
                peak_coordsy_Vec_Shifted[j_global] = peak2[1] - rsz

            EMData.write_images(gain_shifted_file + ":4", shifted_chunk, start)
            np.savetxt(fn_shifts,
                       np.column_stack((peak_coordsx_Vec, peak_coordsy_Vec)),
                       fmt="%.6f")
            np.savetxt(fn_shifts_after,
                       np.column_stack((peak_coordsx_Vec_Shifted,
                                        peak_coordsy_Vec_Shifted)),
                       fmt="%.6f")

        # Shift histograms
        for vec_x, vec_y, label in [
            (peak_coordsx_Vec,         peak_coordsy_Vec,         "Shifts"),
            (peak_coordsx_Vec_Shifted, peak_coordsy_Vec_Shifted, "ShiftsAfter"),
        ]:
            plt.figure(figsize=(15, 5))
            plt.subplot(1, 2, 1)
            plt.hist(vec_x, 50, log=True)
            plt.title(f"x {label} {movie_string}")
            plt.subplot(1, 2, 2)
            plt.hist(vec_y, 50, log=True)
            plt.title(f"y {label}")
            plt.savefig(f"peaks/{label}_{movie_string}.jpg")
            plt.close()

        ts()
        print(f"  Written: {gain_shifted_file}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    import EMAN2
    # EMAN2 import overwrites all _c backup slots — re-apply the fix
    EMData.__initc           = _c_emdata_init
    EMUtil.get_image_count_c = staticmethod(_c_get_image_count)
    EMData.read_image_c      = _c_read_image
    EMData.read_images_c     = staticmethod(_c_read_images)
    EMData.write_image_c     = _c_write_image
    EMData.write_images_c    = staticmethod(_c_write_images)

    progname = os.path.basename(sys.argv[0])
    usage = """e3ctomo_import.py movies [options]

Continuous-tilt import: copy, thumbnail, time-average, gain-correct, and CCF-align movies.
Supply movies as interleaved neg/pos pairs: neg1 pos1 neg2 pos2 ...
"""
    parser = EMAN2.EMArgumentParser(usage=usage, version=EMAN2.EMANVERSION)
    parser.add_pos_argument(name="movies", help="Input movie HDF files (interleaved neg/pos pairs).", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=True)", row=0, col=0, rowspan=1, colspan=3)
    parser.add_argument("--gain", type=str, default="GainFiducial.hdf", help="Gain reference HDF file.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=1, col=0, rowspan=1, colspan=3)
    parser.add_argument("--outdir", type=str, default=".", help="Project output directory.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", filecheck=False, row=2, col=0, rowspan=1, colspan=3)
    parser.add_header(name="apixhdr", help="", title=u"Å/pixel", row=3, col=0, rowspan=1, colspan=1)
    parser.add_argument("--apix", type=float, default=3.714, help=u"Å/pixel — Angstroms per pixel.", guitype='floatbox', row=3, col=1, rowspan=1, colspan=1)
    parser.add_argument("--meanshrink", type=int, default=1, help="Spatial shrink factor (1, 2, or 4; default: 1).", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
    parser.add_argument("--compressbits", type=int, default=4, help="Compression bits for movie copies (default: 4).", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
    parser.add_argument("--ntave", type=int, default=18, help="Rolling time-average window size (default: 18).", guitype='intbox', row=4, col=2, rowspan=1, colspan=1)
    parser.add_argument("--chunk", type=int, default=50, help="I/O chunk size for processing (default: 50).")
    parser.add_argument("--skip_copy", default=False, help="Skip step i — movies already present in outdir/movies/.", action="store_true", guitype='boolbox', row=5, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_thumbs", default=False, help="Skip step iii — thumbnail creation.", action="store_true", guitype='boolbox', row=5, col=1, rowspan=1, colspan=1)
    parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross-platform PPID.")
    (options, args) = parser.parse_args()

    logid = EMAN2.E2init(sys.argv, options.ppid)

    input_movies = [os.path.abspath(m) for m in args] if args else []
    gain_file    = os.path.abspath(options.gain)
    outdir       = os.path.abspath(options.outdir)

    if options.meanshrink not in VALID_MEANSHRINK:
        print(f"ERROR: --meanshrink must be one of {sorted(VALID_MEANSHRINK)}, "
              f"got {options.meanshrink}")
        sys.exit(1)

    if not os.path.isfile(gain_file):
        print(f"ERROR: gain file not found: {gain_file}")
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)
    print(f"Working directory: {os.getcwd()}")
    ts()

    # movie_files always mirrors input_movies, just inside outdir/movies/
    movies_dir  = os.path.join(outdir, "movies")
    movie_files = [os.path.join(movies_dir, os.path.basename(m)) for m in input_movies]

    # ── Step i ────────────────────────────────────────────────────────────────
    if not options.skip_copy:
        print("\n--- Step i: copying movies to movies/ ---")
        step_i_copy_movies(input_movies, outdir, options.apix, options.meanshrink, options.compressbits)
    else:
        print("Skipping step i (--skip_copy).")
        missing = [f for f in movie_files if not os.path.isfile(f)]
        if missing:
            for f in missing:
                print(f"ERROR: expected movie not found: {f}")
            sys.exit(1)

    # ── Step iii ──────────────────────────────────────────────────────────────
    if not options.skip_thumbs:
        print("\n--- Step iii: creating thumbnails ---")
        step_iii_thumbnails(movie_files, outdir, options.meanshrink, options.ntave)
        if len(movie_files) >= 2:
            print("\n--- Step iii (combine): merging neg+pos thumbnail pairs ---")
            step_iii_combine_thumbnails(movie_files, outdir)
    else:
        print("Skipping step iii (--skip_thumbs).")

    # ── Step iv ───────────────────────────────────────────────────────────────
    print("\n--- Step iv: rolling time-average + gain correction + CCF alignment ---")
    step_iv_gain_and_align(
        movie_files,
        gain_file=gain_file,
        apix=options.apix,
        number_to_be_averaged=options.ntave,
        chunkframes=options.chunk,
    )

    print("\nAll done.")
    ts()

    EMAN2.E2end(logid)


if __name__ == "__main__":
    main()
