"""
ct_functions.py  --  Helper functions for continuous-tilt cryo-EM processing.
"""

import os
import numpy as np
import h5py
from EMAN3 import EMData, EMUtil, Averagers, Region

def hdf_nimg(path):
    with h5py.File(path, 'r') as f:
        return len(f['MDF']['images'])


def do_time_average_chunked_rolling_gain(
    fsp,
    number_to_be_averaged,
    chunk=50,
    out=None,
    verbose=True,
    gain_fsp=None,
):
    """
    Local time-average of a large stack using overlapped chunked I/O
    and a rolling sum (center excluded). Optional gain correction of outputs.

    Parameters
    ----------
    fsp : str
        Input HDF stack filename.
    number_to_be_averaged : int
        Total window size INCLUDING center (e.g., 18 -> half window = 9).
        Center is excluded from the mean when computing the output.
        Should be even for symmetry.
    chunk : int
        Number of *center* frames to process per iteration (no overlap counted here).
    out : str or None
        Output HDF filename; if None, auto-generate from input.
    verbose : bool
        Print progress every ~100 frames.
    gain_fsp : str or None
        Optional path to gain reference HDF file.
    """
    nimg = hdf_nimg(fsp)
    half = number_to_be_averaged // 2
    if out is None:
        out = fsp[:-4] + f"_GainTimeAve{number_to_be_averaged}.hdf"

    gain_im = EMData(gain_fsp) if gain_fsp is not None else None

    for start in range(0, nimg, chunk):
        stop = min(start + chunk, nimg)

        lo = max(0, start - half)
        hi = min(nimg, stop + half)
        idx = list(range(lo, hi))

        imgs = EMData.read_images(fsp, idx)
        m = len(imgs)

        def window_bounds(j):
            j_off = j - lo
            wlo = max(0, j_off - half)
            whi = min(m, j_off + half + 1)
            return j_off, wlo, whi

        outs = []
        outs_start = start

        j = start
        j_off, wlo, whi = window_bounds(j)

        sum_all = imgs[wlo].copy()
        sum_all.mult(0.0)
        for k in range(wlo, whi):
            sum_all += imgs[k]

        prev_wlo, prev_whi = wlo, whi

        count = (whi - wlo) - 1
        out_img = sum_all.copy()
        out_img -= imgs[j_off]
        if count > 0:
            out_img.mult(1.0 / float(count))

        if gain_im is not None:
            out_img.process_inplace("math.fixgain.counting", {"gain": gain_im, "gainmin": 2, "gainmax": 2})

        outs.append(out_img)
        if verbose and (j % 100 == 0):
            print(f"Processed frame {j}/{nimg}")

        for j in range(start + 1, stop):
            j_off, wlo, whi = window_bounds(j)

            while prev_wlo < wlo:
                sum_all -= imgs[prev_wlo]
                prev_wlo += 1
            while prev_wlo > wlo:
                prev_wlo -= 1
                sum_all += imgs[prev_wlo]
            while prev_whi < whi:
                sum_all += imgs[prev_whi]
                prev_whi += 1
            while prev_whi > whi:
                prev_whi -= 1
                sum_all -= imgs[prev_whi]

            count = (whi - wlo) - 1
            out_img = sum_all.copy()
            out_img -= imgs[j_off]
            if count > 0:
                out_img.mult(1.0 / float(count))

            if gain_im is not None:
                out_img.process_inplace("math.fixgain.counting", {"gain": gain_im, "gainmin": 2, "gainmax": 2})

            outs.append(out_img)
            if verbose and (j % 100 == 0):
                print(f"Processed frame {j}/{nimg}")

        EMData.write_images(out + ":6", outs, outs_start)

    if verbose:
        print("Done.")


def do_time_average_chunked_rolling(
    fsp,
    number_to_be_averaged,
    chunk=50,
    out=None,
    verbose=True,
    compress_level=1
):
    """
    Local time-average of a large HDF stack using overlapped chunked I/O
    and a rolling sum (center excluded).

    Parameters
    ----------
    fsp : str
        Input HDF stack filename.
    number_to_be_averaged : int
        Total window size INCLUDING center (e.g., 18 -> half window = 9).
        Center is excluded from the mean when computing the output.
        Should be even for symmetry.
    chunk : int
        Number of *center* frames to process per iteration.
    out : str or None
        Output HDF filename; if None, auto-generate from input.
    verbose : bool
        Print progress.
    """
    nimg = hdf_nimg(fsp)
    half = number_to_be_averaged // 2
    if out is None:
        out = fsp[:-4] + f"TimeAve{number_to_be_averaged}.hdf"

    for start in range(0, nimg, chunk):
        stop = min(start + chunk, nimg)

        lo = max(0, start - half)
        hi = min(nimg, stop + half)
        idx = list(range(lo, hi))

        imgs = EMData.read_images(fsp, idx)
        m = len(imgs)

        def window_bounds(j):
            j_off = j - lo
            wlo = max(0, j_off - half)
            whi = min(m, j_off + half + 1)
            return j_off, wlo, whi

        outs = []
        outs_start = start

        j = start
        j_off, wlo, whi = window_bounds(j)

        sum_all = imgs[wlo].copy()
        sum_all.mult(0.0)
        for k in range(wlo, whi):
            sum_all += imgs[k]

        prev_wlo, prev_whi = wlo, whi

        count = (whi - wlo) - 1
        out_img = sum_all.copy()
        out_img -= imgs[j_off]
        if count > 0:
            out_img.mult(1.0 / float(count))
        outs.append(out_img)

        if verbose and (j % chunk == 0):
            print(f"Processed frame {j}/{nimg}")

        for j in range(start + 1, stop):
            j_off, wlo, whi = window_bounds(j)

            while prev_wlo < wlo:
                sum_all -= imgs[prev_wlo]
                prev_wlo += 1
            while prev_wlo > wlo:
                prev_wlo -= 1
                sum_all += imgs[prev_wlo]
            while prev_whi < whi:
                sum_all += imgs[prev_whi]
                prev_whi += 1
            while prev_whi > whi:
                prev_whi -= 1
                sum_all -= imgs[prev_whi]

            count = (whi - wlo) - 1
            out_img = sum_all.copy()
            out_img -= imgs[j_off]
            if count > 0:
                out_img.mult(1.0 / float(count))
            outs.append(out_img)

            if verbose and (j % chunk == 0):
                print(f"Processed frame {j}/{nimg}")

        EMData.write_images(out + ":6", outs, outs_start, compress_level=compress_level)

    if verbose:
        print("Done.")


def subtract_interior_mean(shrunk_img):
    arr = shrunk_img.numpy()
    h, w = arr.shape
    cy, cx = h // 2, w // 2
    hh, hw = h // 4, w // 4
    interior = arr[cy - hh: cy + hh, cx - hw: cx + hw]
    interior_mean = interior.mean()
    return shrunk_img - interior_mean


def gain_correct(input_filename, gain_filename):
    gc_output_filename = input_filename[:-4] + "_gain.hdf:6"
    print(gc_output_filename)
    gain_output_EMData = EMData(gain_filename)

    Count = -1
    while 1:
        Count += 1
        try:
            avg_EMdata = EMData(input_filename, Count)
            gain_corrected = gain_output_EMData * avg_EMdata
            gain_corrected.write_image(gc_output_filename, -1)
        except:
            print(Count)
            break


def AverageSequentialFrames_GC_Shrink(fname, fnameRootwoSub, sub_average, gain_filename):
    fnameRoot = fnameRootwoSub + str(sub_average)
    Output_fname_ave = fnameRoot + ".hdf"
    CommandNow = "e2proc2d.py {} {} --avgseq {}".format(fname, Output_fname_ave, sub_average)
    print(CommandNow)
    os.system(CommandNow)

    gain_correct(Output_fname_ave, gain_filename)
    gc_output_filename = Output_fname_ave[:-4] + "_gain.hdf"
    print(gc_output_filename)

    zero_output_fn_gainCorrected_AndShrunk = gc_output_filename[:-4] + "Sh.hdf"
    CommandNow = "e2proc2d.py {} {} --fouriershrink 2".format(gc_output_filename, zero_output_fn_gainCorrected_AndShrunk)
    print(CommandNow)
    os.system(CommandNow)


def DoTimeAverage(FileNameToBeTimeAveraged, NumberToBeAveraged):
    nimg = hdf_nimg(FileNameToBeTimeAveraged)
    NumberOnEachSide = NumberToBeAveraged // 2
    FileNameTimeAveraged18 = FileNameToBeTimeAveraged[:-4] + 'TimeAve' + str(NumberToBeAveraged) + '.hdf'

    try:
        os.unlink(FileNameTimeAveraged18)
    except:
        pass

    print('nimg is ' + str(nimg))

    for jImg in range(nimg):
        if jImg % 10 == 0:
            print('jImg is ' + str(jImg))

        jLo = jImg - NumberOnEachSide
        jMin = max(0, jLo)
        jMax = jMin + NumberToBeAveraged + 1

        if jMax > nimg:
            jMax = nimg
            jMin = jMax - NumberToBeAveraged - 1

        TimeAvg = Averagers.get("mean")
        for jLoop in range(NumberToBeAveraged + 1):
            jIndex = jMin + jLoop
            if jIndex == jImg:
                continue
            Frame4Kby4K = EMData(FileNameToBeTimeAveraged, jIndex)
            TimeAvg.add_image(Frame4Kby4K)

        Data2WriteOut = TimeAvg.finish()
        Data2WriteOut.write_image(FileNameTimeAveraged18 + ":6", jImg)


def center_clip(emdata_obj, size):
    """
    Clips the center of an EMData object to the specified size.

    Parameters
    ----------
    emdata_obj : EMData
    size : int or tuple of (width, height)

    Returns
    -------
    EMData clipped to center region.
    """
    try:
        size = (int(size), int(size))
    except TypeError:
        size = (int(size[0]), int(size[1]))

    if size[0] < 1 or size[1] < 1:
        raise ValueError("center_clip(size) must be called with a positive integer")

    x_offset = (emdata_obj.get_xsize() - size[0]) // 2
    y_offset = (emdata_obj.get_ysize() - size[1]) // 2

    if x_offset < 0 or y_offset < 0:
        raise ValueError("Requested size is larger than the original dimensions.")

    return emdata_obj.get_clip(Region(x_offset, y_offset, size[0], size[1]))


def average_image_stack(hdf_filename):
    first_image = EMData(hdf_filename, 0)
    n_images = hdf_nimg(hdf_filename)

    sum_image = EMData(first_image.get_xsize(), first_image.get_ysize(), 1)

    for i in range(n_images):
        image = EMData(hdf_filename, i)
        sum_image.add(image)

    sum_image.mult(1.0 / n_images)
    return sum_image


def modify_hdf_headers(hdf_file, new_header_info):
    """
    Modifies the header information for each image in an HDF stack.

    Parameters
    ----------
    hdf_file : str
        Path to the HDF file.
    new_header_info : dict
        Header attributes and values to modify.
    """
    nimgs = hdf_nimg(hdf_file)
    for i in range(nimgs):
        img = EMData(hdf_file, i)
        for key, value in new_header_info.items():
            img.set_attr(key, value)
        img.write_image(hdf_file + ":6", i)
