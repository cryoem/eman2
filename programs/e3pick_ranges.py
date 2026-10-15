#!/usr/bin/env python
"""
e3pick_ranges.py  --  Pick good frame ranges from a continuous-tilt thumbnail movie.

The thumbnail movie is a single HDF stack: [neg frames, reversed] + [pos frames].

Usage
-----
  e3pick_ranges.py thumbnailmovie [--full_movie MOVIE] [--avgseq N] [--unidirectional]

Controls (bidirectional)
------------------------
  Frame slider     : scrub through all frames
  Left/Right arrow : step one frame
  1  Set Start (Neg)    2  Set End (Neg)
  3  Set Start (Pos)    4  Set End (Pos)
  Q  Done

Controls (unidirectional)
--------------------------
  S  Set Start    E  Set End    Q  Done
"""

import sys, os, json

# Pre-load import_params.json so parser defaults can reference import_params_data
import_params_data = {}
try:
    _ip = os.path.join(os.getcwd(), "info", "import_params.json")
    if os.path.isfile(_ip):
        with open(_ip) as _f:
            import_params_data = json.load(_f)
except Exception:
    pass

import EMAN2_cppwrap as _cppwrap
_c_emdata_init        = _cppwrap.EMData.__init__
_c_get_image_count    = _cppwrap.EMUtil.get_image_count
_c_read_image         = _cppwrap.EMData.read_image
_c_read_images        = _cppwrap.EMData.read_images
_c_write_image        = _cppwrap.EMData.write_image

from EMAN3 import EMData, EMUtil
from eman2_gui.emapplication import EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.valslider import ValSlider

EMData.__initc           = _c_emdata_init
EMUtil.get_image_count_c = staticmethod(_c_get_image_count)
EMData.read_image_c      = _c_read_image
EMData.read_images_c     = staticmethod(_c_read_images)
EMData.write_image_c     = _c_write_image

# Match whatever Qt binding EMApp (and this environment's eman2_gui) is built on --
# this script is symlinked across multiple conda envs that don't necessarily agree
# on PyQt5 vs PySide6, so a hardcoded import breaks whichever env doesn't match.
if EMApp.__bases__[0].__module__.split('.')[0] == "PySide6":
    from PySide6 import QtWidgets, QtCore
else:
    from PyQt5 import QtWidgets, QtCore
import h5py


def hdf_nimg(path):
    with h5py.File(path, 'r') as f:
        return len(f['MDF']['images'])


# ---------------------------------------------------------------------------
# Path helpers (mirrors e3ctomo_reconstruct.py)
# ---------------------------------------------------------------------------

def _lst_path(neg_gain_cor, pos_gain_cor, basename=None):
    if basename:
        return os.path.join("sets", f"{basename}_combined.lst")
    neg_base = os.path.splitext(os.path.basename(neg_gain_cor))[0]
    pos_base = os.path.splitext(os.path.basename(pos_gain_cor))[0]
    return os.path.join("sets", f"{neg_base}__{pos_base}__combined.lst")

def _sa_path(lst_path, avgseq, meanshrink=4, basename=None):
    tag = f"_Sh{meanshrink}" if meanshrink > 1 else ""
    if basename:
        return os.path.join("tiltseries", f"{basename}__SA{avgseq}{tag}.hdf")
    stem = os.path.splitext(os.path.basename(lst_path))[0]
    return os.path.join("tiltseries", f"{stem}__SA{avgseq}{tag}.hdf")


# ---------------------------------------------------------------------------
# Post-done pipeline: steps 1 and 2
# ---------------------------------------------------------------------------

def _even_chunks(n, avgseq):
    """Split n consecutive frames into chunks whose sizes are all within one of avgseq.

    A fixed stride (range(0, n, avgseq)) leaves the remainder of the half as its
    own chunk. Since the neg/pos seam fix subaverages each half separately, that
    now happens twice -- once per half -- and the leftover can be a single frame.
    Fiducial2 (pos = 2701 = 45*60 + 1) produced a final tilt built from ONE frame,
    measured at 5.9x the noise of its neighbours; Fiducial4 (neg = 3961 = 66*60 + 1)
    would have done the same at the tilt nearest zero degrees.

    Choosing the chunk COUNT by rounding and then distributing the frames as evenly
    as possible keeps every chunk within one frame of avgseq, discards nothing, and
    reduces to plain fixed-stride chunking whenever avgseq divides n exactly.

    Returns a list of (first, last) inclusive index pairs relative to the half.
    """
    if n <= 0:
        return []
    k = max(1, int(round(n / float(avgseq))))
    base, rem = divmod(n, k)          # rem chunks get base+1 frames, the rest base
    sizes = [base + 1] * rem + [base] * (k - rem)
    out, s = [], 0
    for sz in sizes:
        out.append((s, s + sz - 1))
        s += sz
    return out

def run_steps_1_2(options):
    """Build combined LST (step 1) and SA tilt series (step 2)."""
    import json

    neg_gain_cor = os.path.abspath(options.neg_gain_cor)
    pos_gain_cor = os.path.abspath(options.pos_gain_cor)
    bn        = options.basename or None
    avgseq    = options.avgseq_sa
    cbits     = options.compressbits
    shrink    = options.meanshrink

    # Prefer THIS dataset's own archived ranges. info/pick_ranges.json is
    # overwritten by whichever dataset was picked last, so a deferred or batched
    # run (see --nogui) must not read it blindly -- that silently builds the LST
    # for the wrong dataset. Same basename-aware rule as e3ctomo_reconstruct.py.
    json_path = os.path.join("info", "pick_ranges.json")
    if bn:
        _archive = os.path.join("info", f"pick_ranges_{bn}.json")
        if os.path.isfile(_archive):
            json_path = _archive
    with open(json_path) as f:
        d = json.load(f)
    _d_bn = d.get("basename")
    if bn and _d_bn and _d_bn != bn:
        print(f"Error: {json_path} holds ranges for '{_d_bn}', not '{bn}'. "
              f"Re-pick ranges for {bn}, or restore info/pick_ranges_{bn}.json.")
        return
    print(f"  frame ranges from {json_path}")
    neg_start = d.get("neg_start", 0)
    neg_end   = d.get("neg_end",   9999)
    pos_start = d.get("pos_start", 0)
    pos_end   = d.get("pos_end",   9999)

    lst  = _lst_path(neg_gain_cor, pos_gain_cor, bn)
    sa   = _sa_path(lst, avgseq, shrink, bn)

    # Step 1 — combined LST
    n_neg = EMUtil.get_image_count(neg_gain_cor)
    neg_start = max(0, neg_start)
    neg_end   = min(n_neg - 1, neg_end)
    pos_end_max = n_neg + EMUtil.get_image_count(pos_gain_cor) - 1
    pos_start = max(n_neg, pos_start)
    pos_end   = min(pos_end_max, pos_end)
    if getattr(options, 'skip_lst', False):
        print(f"\n--- Step 1 skipped (--skip_lst)  [{lst}]")
    else:
        print("\n=== Step 1 — build combined LST ===")
        os.makedirs("sets", exist_ok=True)
        with open(lst, "w") as f:
            f.write("#LST\n")
            for i in range(neg_start, neg_end + 1):
                f.write(f"{n_neg - 1 - i}\t{neg_gain_cor}\n")
            for i in range(pos_start, pos_end + 1):
                f.write(f"{i - n_neg}\t{pos_gain_cor}\n")
        print(f"  Wrote {lst}")
        print(f"  Neg frames: {neg_end - neg_start + 1}  "
              f"(combined {neg_start}–{neg_end} → orig neg {n_neg-1-neg_start}–{n_neg-1-neg_end})")
        print(f"  Pos frames: {pos_end - pos_start + 1}  "
              f"(pos movie {pos_start-n_neg}–{pos_end-n_neg})")

    # Step 2 — sub-average (+ optional spatial shrink)
    if getattr(options, 'skip_sa', False):
        print(f"\n--- Step 2 skipped (--skip_sa)  [{sa}]")
    else:
        print(f"\n=== Step 2 — sub-average (meanshrink={shrink}) ===")
        os.makedirs("tiltseries", exist_ok=True)
        if os.path.isfile(sa):
            os.remove(sa)
            print(f"  Removed existing: {sa}")
        hdr = EMData(neg_gain_cor, 0, True)
        output_apix = hdr["apix_x"] * shrink
        shrink_arg = f" --meanshrink={shrink}" if shrink > 1 else ""

        # Subaverage each half SEPARATELY.
        #
        # e2proc2d --avgseq chunks the list in fixed groups, so a single pass over
        # the combined LST produces one chunk that straddles the neg/pos boundary
        # whenever the neg half is not a whole multiple of avgseq. That chunk
        # averages frames from two different acquisitions, separated by a stage
        # return and whatever drift occurred between them: a double image, sitting
        # at the tilt nearest zero degrees where it carries the most weight in the
        # alignment. Observed on Fiducial1 and Fiducial3 (both 31 neg + 29 pos),
        # where it left the two halves misregistered by 195 px and 60 px because
        # the only frame linking them was the corrupted one. Fiducial4 (1 neg +
        # 59 pos) was unaffected and registered to 2 px.
        #
        # The averaging is done here rather than by e2proc2d --avgseq so that the
        # chunk boundaries come from _even_chunks(): a fixed stride leaves each
        # half's remainder as its own chunk, which can be a single frame. See
        # _even_chunks() for the two datasets that hit this.
        n_sel_neg = neg_end - neg_start + 1
        n_sel_pos = pos_end - pos_start + 1
        with open(lst) as _f:
            _rows = [ln for ln in _f if not ln.startswith("#")]

        chunks = []          # LST row ranges, recorded for the angle model
        ok = True
        _idx = 0
        for _base, _n in ((0, n_sel_neg), (n_sel_neg, n_sel_pos)):
            if _n <= 0:
                continue
            _ranges = _even_chunks(_n, avgseq)
            print(f"  {'neg' if _base == 0 else 'pos'} half: {_n} frames -> "
                  f"{len(_ranges)} chunks of "
                  f"{min(b-a+1 for a,b in _ranges)}-{max(b-a+1 for a,b in _ranges)}")
            for _a, _b in _ranges:
                _acc = None
                for _row in _rows[_base + _a: _base + _b + 1]:
                    _fld = _row.split()
                    _im = EMData(_fld[1], int(_fld[0]))
                    if _acc is None:
                        _acc = _im
                    else:
                        _acc.add(_im)
                if _acc is None:
                    print(f"  WARNING: empty chunk at rows {_base+_a}-{_base+_b}")
                    ok = False
                    continue
                _acc.mult(1.0 / float(_b - _a + 1))
                if shrink > 1:
                    _acc.process_inplace("math.meanshrink", {"n": shrink})
                _acc["apix_x"] = _acc["apix_y"] = _acc["apix_z"] = output_apix
                _acc.write_image(f"{sa}:{cbits}", _idx)
                chunks.append([_base + _a, _base + _b])
                _idx += 1

        _sz = [c[1] - c[0] + 1 for c in chunks]
        _nneg = sum(1 for c in chunks if c[0] < n_sel_neg)
        print(f"  Subaveraged each half separately: {len(chunks)} tilts "
              f"({_nneg} neg + {len(chunks) - _nneg} pos), "
              f"chunk sizes {min(_sz)}-{max(_sz)} (target {avgseq}), "
              f"no chunk crosses the neg/pos boundary")
        if _idx != len(chunks):
            print(f"  WARNING: wrote {_idx} images but expected {len(chunks)} chunks")

        # Record the chunk boundaries so the angle model uses the same scheme.
        # Absent this file, e3ctomo_reconstruct.py falls back to uniform chunking
        # over the combined list, which is what older data was built with.
        if ok and bn:
            os.makedirs("info", exist_ok=True)
            _cj = os.path.join("info", f"sa_chunks_{bn}.json")
            with open(_cj, "w") as _f:
                json.dump({"basename": bn, "avgseq_sa": avgseq,
                           "n_sel_neg": n_sel_neg, "n_sel_pos": n_sel_pos,
                           "scheme": "per_half", "chunks": chunks}, _f, indent=1)
            print(f"  Wrote {_cj}")

    # Step 2b — selected thumbnail
    thumb = getattr(options, 'thumbnailmovie', None)
    if getattr(options, 'skip_thumb', False):
        print(f"\n--- Step 2b skipped (--skip_thumb)")
    elif bn and thumb:
        ntave = options.ntave or 18
        n_neg_thumb = n_neg // ntave
        t_neg_start = neg_start // ntave
        t_neg_end   = neg_end   // ntave
        t_pos_start = (pos_start - n_neg) // ntave + n_neg_thumb
        t_pos_end   = (pos_end   - n_neg) // ntave + n_neg_thumb
        print(f"\n=== Step 2b — selected thumbnail ===")
        make_selected_thumb(thumb, t_neg_start, t_neg_end, t_pos_start, t_pos_end, bn)


def make_selected_thumb(thumbnailmovie, neg_start, neg_end, pos_start, pos_end, basename):
    """Write thumbnails/<basename>_selected_thumb.hdf with only the selected frames."""
    if not thumbnailmovie or not os.path.isfile(thumbnailmovie):
        print("  Selected thumbnail skipped (combined thumbnail not found).")
        return
    os.makedirs("thumbnails", exist_ok=True)
    out = os.path.join("thumbnails", f"{basename}_selected_thumb.hdf")
    if os.path.isfile(out):
        os.remove(out)
    out_idx = 0
    for i in range(neg_start, neg_end + 1):
        img = EMData(thumbnailmovie, i)
        img.write_image_c(out, out_idx)
        out_idx += 1
    for i in range(pos_start, pos_end + 1):
        img = EMData(thumbnailmovie, i)
        img.write_image_c(out, out_idx)
        out_idx += 1
    print(f"  Selected thumbnail → {out}  ({out_idx} frames)")


def load_stack(path):
    n = hdf_nimg(path)
    print(f"  Loading {n} frames from {os.path.basename(path)} ...")
    frames = [EMData(path, i) for i in range(n)]
    return frames


class PickRangesWindow(QtWidgets.QMainWindow):

    def __init__(self, frames, args):
        super().__init__()
        self.frames = frames
        self.args = args
        self.done_clicked = False
        self.processing = False
        self.unidirectional = args.unidirectional
        self.avgseq = args.ntave or 1
        self.idx = 0
        self._markers_set = set()
        n = len(frames)

        if self.unidirectional:
            self.start = 0
            self.end = n - 1
        else:
            half = n // 2
            self.neg_start = 0
            self.neg_end   = half - 1
            self.pos_start = half
            self.pos_end   = n - 1

        self.setWindowTitle("e3pick_ranges  —  frame range picker")
        self.resize(700, 820)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        vbox = QtWidgets.QVBoxLayout(central)
        vbox.setSpacing(4)

        self.fname_label = QtWidgets.QLabel(os.path.basename(args.thumbnailmovie))
        self.fname_label.setAlignment(QtCore.Qt.AlignCenter)
        self.fname_label.setStyleSheet("font-size: 11pt; font-weight: bold;")
        vbox.addWidget(self.fname_label)

        self.image_widget = EMImage2DWidget()
        self.image_widget.setMinimumSize(512, 512)
        vbox.addWidget(self.image_widget, stretch=1)

        self.frame_label = QtWidgets.QLabel("")
        self.frame_label.setAlignment(QtCore.Qt.AlignCenter)
        self.frame_label.setStyleSheet("font-size: 12pt; font-weight: bold;")
        vbox.addWidget(self.frame_label)

        slider_row = QtWidgets.QHBoxLayout()
        slider_row.setContentsMargins(0, 0, 0, 0)
        self.frame_slider = ValSlider(rng=(0, n - 1), label="Frame", value=0, labelwidth=50)
        self.frame_slider.text.setFixedWidth(40)   # half of ValSlider default (80px)
        self.frame_slider.valueChanged.connect(self._on_slider)
        slider_row.addWidget(self.frame_slider)
        vbox.addLayout(slider_row)

        # "Kept Frames" label width matches ValSlider's left portion:
        #   label(50) + spacing(6) + textbox(40) + spacing(6) = 102px
        slider_row2 = QtWidgets.QHBoxLayout()
        slider_row2.setContentsMargins(0, 0, 0, 0)
        slider_row2.setSpacing(0)
        kept_label = QtWidgets.QLabel("Kept Frames")
        kept_label.setFixedWidth(102)
        kept_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.range_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.range_slider.setRange(0, n - 1)
        self.range_slider.setValue(0)
        self.range_slider.setEnabled(False)
        self.range_slider.valueChanged.connect(self._on_range_slider)
        slider_row2.addWidget(kept_label)
        slider_row2.addWidget(self.range_slider)
        vbox.addLayout(slider_row2)

        self.range_label = QtWidgets.QLabel("")
        self.range_label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(self.range_label)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.setSpacing(6)

        self._btn_active = {}   # btn -> bool

        if self.unidirectional:
            self.btn_start = QtWidgets.QPushButton("Set Start  [S]")
            self.btn_end   = QtWidgets.QPushButton("Set End  [E]")
            self.btn_done  = QtWidgets.QPushButton("Done  [Q]")
            self.btn_start.clicked.connect(self._set_start)
            self.btn_end.clicked.connect(self._set_end)
            self.btn_done.clicked.connect(self._done)
            self._btn_active[self.btn_start] = False
            self._btn_active[self.btn_end]   = False
            for b in (self.btn_start, self.btn_end, self.btn_done):
                btn_row.addWidget(b)
            hint_text = "Left/Right arrow keys: step one frame  |  S = start    E = end    Q = done"
        else:
            self.btn_neg_start = QtWidgets.QPushButton("Set Start (Neg)  [1]")
            self.btn_neg_end   = QtWidgets.QPushButton("Set End (Neg)  [2]")
            self.btn_pos_start = QtWidgets.QPushButton("Set Start (Pos)  [3]")
            self.btn_pos_end   = QtWidgets.QPushButton("Set End (Pos)  [4]")
            self.btn_done      = QtWidgets.QPushButton("Done  [Q]")
            self.btn_neg_start.clicked.connect(self._set_neg_start)
            self.btn_neg_end.clicked.connect(self._set_neg_end)
            self.btn_pos_start.clicked.connect(self._set_pos_start)
            self.btn_pos_end.clicked.connect(self._set_pos_end)
            self.btn_done.clicked.connect(self._done)
            self._btn_active[self.btn_neg_start] = False
            self._btn_active[self.btn_neg_end]   = False
            self._btn_active[self.btn_pos_start] = False
            self._btn_active[self.btn_pos_end]   = False
            for b in (self.btn_neg_start, self.btn_neg_end,
                      self.btn_pos_start, self.btn_pos_end, self.btn_done):
                btn_row.addWidget(b)
            hint_text = "Left/Right arrow keys: step one frame  |  1=neg start  2=neg end  3=pos start  4=pos end  Q=done"

        vbox.addLayout(btn_row)

        hint = QtWidgets.QLabel(hint_text)
        hint.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(hint)

        self._try_load_saved_ranges()
        self._update_display()

    # ------------------------------------------------------------------

    def _try_load_saved_ranges(self):
        """Reload previously saved ranges from info/pick_ranges_<basename>.json or pick_ranges.json."""
        json_path = None
        if self.args.basename:
            archive = os.path.join("info", f"pick_ranges_{self.args.basename}.json")
            if os.path.isfile(archive):
                json_path = archive
        if json_path is None:
            flat = os.path.join("info", "pick_ranges.json")
            if os.path.isfile(flat):
                json_path = flat
        if json_path is None:
            return
        try:
            with open(json_path) as f:
                d = json.load(f)
        except Exception:
            return
        if d.get("unidirectional") != self.unidirectional:
            return
        a = self.avgseq
        n = len(self.frames)
        if self.unidirectional:
            if "start" in d and "end" in d:
                self.start = max(0, min(n - 1, d["start"] // a))
                self.end   = max(0, min(n - 1, d["end"]   // a))
                self._markers_set.update({"start", "end"})
                self._btn_active[self.btn_start] = True
                self._btn_active[self.btn_end]   = True
                self.btn_start.setStyleSheet(self._COLOR_START)
                self.btn_end.setStyleSheet(self._COLOR_END)
        else:
            keys = ("neg_start", "neg_end", "pos_start", "pos_end")
            if all(k in d for k in keys):
                half = n // 2
                self.neg_start = max(0,      min(half - 1, d["neg_start"] // a))
                self.neg_end   = max(0,      min(half - 1, d["neg_end"]   // a))
                self.pos_start = max(half,   min(n - 1,    d["pos_start"] // a))
                self.pos_end   = max(half,   min(n - 1,    d["pos_end"]   // a))
                self._markers_set.update({"neg_start", "neg_end", "pos_start", "pos_end"})
                self._btn_active[self.btn_neg_start] = True
                self._btn_active[self.btn_neg_end]   = True
                self._btn_active[self.btn_pos_start] = True
                self._btn_active[self.btn_pos_end]   = True
                self.btn_neg_start.setStyleSheet(self._COLOR_START)
                self.btn_neg_end.setStyleSheet(self._COLOR_END)
                self.btn_pos_start.setStyleSheet(self._COLOR_START)
                self.btn_pos_end.setStyleSheet(self._COLOR_END)
        print(f"  Reloaded ranges from {json_path}")

    # ------------------------------------------------------------------
    _COLOR_START = "background-color: #90EE90;"   # light green
    _COLOR_END   = "background-color: #FF9999;"   # light red
    _COLOR_OFF   = ""

    def _toggle_btn(self, btn, color):
        """Toggle button active state and update its color."""
        active = not self._btn_active.get(btn, False)
        self._btn_active[btn] = active
        btn.setStyleSheet(color if active else self._COLOR_OFF)

    def _all_markers_set(self):
        if self.unidirectional:
            return {'start', 'end'} <= self._markers_set
        return {'neg_start', 'neg_end', 'pos_start', 'pos_end'} <= self._markers_set

    def _valid_frames(self):
        """Sorted list of accessible frame indices (inside selected ranges, gaps excluded)."""
        if self.unidirectional:
            lo, hi = sorted((self.start, self.end))
            return list(range(lo, hi + 1))
        neg_lo, neg_hi = sorted((self.neg_start, self.neg_end))
        pos_lo, pos_hi = sorted((self.pos_start, self.pos_end))
        return list(range(neg_lo, neg_hi + 1)) + list(range(pos_lo, pos_hi + 1))

    def _update_range_slider(self):
        enabled = self._all_markers_set()
        self.range_slider.setEnabled(enabled)
        if not enabled:
            return
        valid = self._valid_frames()
        if not valid:
            return
        n_valid = len(valid)
        if self.idx in valid:
            pos = valid.index(self.idx)
        else:
            pos = min(range(n_valid), key=lambda i: abs(valid[i] - self.idx))
        self.range_slider.blockSignals(True)
        self.range_slider.setRange(0, n_valid - 1)
        self.range_slider.setValue(pos)
        self.range_slider.blockSignals(False)

    def _update_display(self):
        n = len(self.frames)

        self.frame_label.setText(f"Thumbnail frame  {self.idx} / {n-1}")
        self.image_widget.set_data(self.frames[self.idx])

        self.frame_slider.blockSignals(True)
        self.frame_slider.setValue(self.idx)
        self.frame_slider.blockSignals(False)

        if self.unidirectional:
            self.range_label.setText(
                f"Start: {self.start}      End: {self.end}"
            )
        else:
            ns, ne = self.neg_start, self.neg_end
            ps, pe = self.pos_start, self.pos_end
            self.range_label.setText(
                f"Neg:  start {ns}  end {ne}      "
                f"Pos:  start {ps}  end {pe}"
            )

        self._update_range_slider()

    def _on_range_slider(self, val):
        valid = self._valid_frames()
        if not valid:
            return
        self.idx = valid[max(0, min(val, len(valid) - 1))]
        self._update_display()

    def _on_slider(self, val):
        n = len(self.frames)
        self.idx = max(0, min(int(round(val)), n - 1))
        self._update_display()

    def _set_start(self):
        self.start = self.idx
        self._markers_set.add('start')
        self._toggle_btn(self.btn_start, self._COLOR_START)
        self._update_display()

    def _set_end(self):
        self.end = self.idx
        self._markers_set.add('end')
        self._toggle_btn(self.btn_end, self._COLOR_END)
        self._update_display()

    def _set_neg_start(self):
        self.neg_start = self.idx
        self._markers_set.add('neg_start')
        self._toggle_btn(self.btn_neg_start, self._COLOR_START)
        self._update_display()

    def _set_neg_end(self):
        self.neg_end = self.idx
        self._markers_set.add('neg_end')
        self._toggle_btn(self.btn_neg_end, self._COLOR_END)
        self._update_display()

    def _set_pos_start(self):
        self.pos_start = self.idx
        self._markers_set.add('pos_start')
        self._toggle_btn(self.btn_pos_start, self._COLOR_START)
        self._update_display()

    def _set_pos_end(self):
        self.pos_end = self.idx
        self._markers_set.add('pos_end')
        self._toggle_btn(self.btn_pos_end, self._COLOR_END)
        self._update_display()

    def _done(self):
        self.done_clicked = True
        self._print_results()
        # Disable all interactive controls while processing
        self.processing = True
        for btn in list(self._btn_active.keys()) + [self.btn_done]:
            btn.setEnabled(False)
        self.frame_slider.setEnabled(False)
        self.range_slider.setEnabled(False)
        self.setWindowTitle("e3pick_ranges  —  building LST and SA tilt series...")
        self.frame_label.setText("Building tilt series — see terminal for progress")
        self.frame_label.setStyleSheet("font-size: 12pt; font-weight: bold; color: darkorange;")
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.quit()

    def keyPressEvent(self, event):
        if self.processing:
            return
        key = event.key()
        n = len(self.frames)
        if key == QtCore.Qt.Key_Right:
            self.idx = min(self.idx + 1, n - 1)
            self._update_display()
        elif key == QtCore.Qt.Key_Left:
            self.idx = max(self.idx - 1, 0)
            self._update_display()
        elif key == QtCore.Qt.Key_Q:
            self._done()
        elif self.unidirectional:
            if key == QtCore.Qt.Key_S:
                self._set_start()
            elif key == QtCore.Qt.Key_E:
                self._set_end()
            else:
                super().keyPressEvent(event)
        else:
            if key == QtCore.Qt.Key_1:
                self._set_neg_start()
            elif key == QtCore.Qt.Key_2:
                self._set_neg_end()
            elif key == QtCore.Qt.Key_3:
                self._set_pos_start()
            elif key == QtCore.Qt.Key_4:
                self._set_pos_end()
            else:
                super().keyPressEvent(event)

    def closeEvent(self, event):
        if not self.done_clicked:
            self._print_results()
        event.accept()

    def _write_json(self):
        import json
        info_dir = os.path.join(os.getcwd(), "info")
        os.makedirs(info_dir, exist_ok=True)
        path = os.path.join(info_dir, "pick_ranges.json")
        a = self.avgseq
        if self.unidirectional:
            data = {"unidirectional": True,
                    "basename":   self.args.basename,
                    "avgseq_sa":  self.args.avgseq_sa,
                    "meanshrink": self.args.meanshrink,
                    "start": self.start * a,
                    "end":   self.end   * a}
        else:
            data = {"unidirectional": False,
                    "basename":   self.args.basename,
                    "avgseq_sa":  self.args.avgseq_sa,
                    "meanshrink": self.args.meanshrink,
                    "neg_start": self.neg_start * a,
                    "neg_end":   self.neg_end   * a,
                    "pos_start": self.pos_start * a,
                    "pos_end":   self.pos_end   * a}
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        archive = os.path.join(info_dir, f"pick_ranges_{self.args.basename}.json")
        with open(archive, "w") as f:
            json.dump(data, f, indent=2)
        print(f"Ranges written to {path}")

    def _print_results(self):
        self._write_json()
        a = self.avgseq
        if self.unidirectional:
            s, e = self.start, self.end
            print("\n=== Selected ranges (thumbnail frames) ===")
            print(f"  start: {s}  end: {e}  (avgseq={a})")
            print("\n=== Converted to full-stack frames ===")
            print(f"  start: {s*a}  end: {e*a}")
        else:
            ns, ne = self.neg_start, self.neg_end
            ps, pe = self.pos_start, self.pos_end
            print("\n=== Selected ranges (thumbnail frames) ===")
            print(f"  neg: {ns} → {ne}  (avgseq={a})")
            print(f"  pos: {ps} → {pe}  (avgseq={a})")
            print("\n=== Converted to full-stack frames ===")
            print(f"  neg: {ns*a} → {ne*a}")
            print(f"  pos: {ps*a} → {pe*a}")


def main():
    import EMAN2
    progname = os.path.basename(sys.argv[0])
    usage = """e3pick_ranges.py thumbnailmovie [options]

Pick good frame ranges from a continuous-tilt thumbnail movie.
The thumbnail movie is a single HDF stack: [neg frames, reversed] + [pos frames].
"""
    parser = EMAN2.EMArgumentParser(usage=usage, version=EMAN2.EMANVERSION)
    parser.add_pos_argument(name="thumbnailmovie", help="Thumbnail HDF stack (neg reversed + pos concatenated).", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=0, col=0, rowspan=1, colspan=3)
    parser.add_argument("--neg_gain_cor",  type=str, default="", help="Negative tilt GainShifted movie.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=1, col=0, rowspan=1, colspan=3)
    parser.add_argument("--pos_gain_cor",  type=str, default="", help="Positive tilt GainShifted movie.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=2, col=0, rowspan=1, colspan=3)
    parser.add_argument("--basename",   type=str, default="", help="Short label for output files (e.g. CT07).", guitype='strbox', row=3, col=0, rowspan=1, colspan=3)
    parser.add_argument("--ntave",      type=int, default=18, help="Frames averaged per thumbnail — auto-loaded from Gain Correction (do not change).", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
    parser.add_argument("--avgseq_sa",  type=int, default=60, help="Frames to average per tilt for SA tilt series.", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
    parser.add_argument("--meanshrink", type=int, default=4,  help="Spatial shrink for SA tilt series (4 for 4K cameras → 1K; 2 for 2K cameras → 1K).", guitype='intbox', row=4, col=2, rowspan=1, colspan=1)
    parser.add_argument("--compressbits", type=int, default=6, help="Compression bits for SA tilt series.", guitype='intbox', row=5, col=0, rowspan=1, colspan=1)
    parser.add_argument("--unidirectional", default=False, help="Single-direction acquisition (no neg/pos split).", action="store_true", guitype='boolbox', row=5, col=2, rowspan=1, colspan=1)
    parser.add_argument("--skip_lst",   default=False, action="store_true", help="Skip step 1: build combined LST.", guitype='boolbox', row=6, col=0, rowspan=1, colspan=1)
    parser.add_argument("--skip_sa",    default=False, action="store_true", help="Skip step 2: build SA tilt series.", guitype='boolbox', row=6, col=1, rowspan=1, colspan=1)
    parser.add_argument("--skip_thumb", default=False, action="store_true", help="Skip step 2b: write selected thumbnail.", guitype='boolbox', row=6, col=2, rowspan=1, colspan=1)
    parser.add_argument("--nogui", default=False, action="store_true", help="Do not open the picker. Build the LST and SA tilt series from the ranges already saved in info/pick_ranges_<basename>.json (written whenever a picker window is closed, Done or not). Requires --basename, --neg_gain_cor and --pos_gain_cor. Lets the interactive picking and the slow LST/SA builds be run separately, e.g. batched overnight.")
    parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross-platform PPID.")
    (options, args) = parser.parse_args()

    logid = EMAN2.E2init(sys.argv, options.ppid)

    # add_pos_argument does not create the attribute when the positional is
    # omitted, which it legitimately is under --nogui (no picker, no thumbnail
    # needed) -- so default it rather than assuming it exists.
    thumbnailmovie = args[0] if args else getattr(options, "thumbnailmovie", "")
    options.thumbnailmovie = thumbnailmovie

    # Override ntave from import_params.json if present
    import json as _json
    _ip = os.path.join("info", "import_params.json")
    if os.path.isfile(_ip):
        try:
            options.ntave = _json.load(open(_ip)).get("ntave", options.ntave)
        except Exception:
            pass
    options.ntave = options.ntave or 18

    # Headless: the ranges were already saved by an earlier picker session
    # (closing the window writes info/pick_ranges_<basename>.json even when
    # Done was never clicked), so just run the heavy steps. Deliberately placed
    # before load_stack -- reading the whole thumbnail movie is pointless here.
    if options.nogui:
        if not (options.neg_gain_cor and options.pos_gain_cor):
            print("Error: --nogui requires --neg_gain_cor and --pos_gain_cor.")
            EMAN2.E2end(logid)
            return
        if not options.basename:
            print("Error: --nogui requires --basename, so the right saved ranges are used.")
            EMAN2.E2end(logid)
            return
        print(f"--nogui: skipping the picker, building LST + SA for {options.basename}")
        run_steps_1_2(options)
        EMAN2.E2end(logid)
        return

    print(f"Loading {thumbnailmovie} ...")
    frames = load_stack(thumbnailmovie)

    app = EMApp()
    win = PickRangesWindow(frames, options)
    win.show()
    app.exec_()

    if win.done_clicked and options.neg_gain_cor and options.pos_gain_cor:
        run_steps_1_2(options)
    elif win.done_clicked and not (options.neg_gain_cor and options.pos_gain_cor):
        print("Note: --neg_gain_cor/--pos_gain_cor not supplied; skipping LST and SA steps.")

    EMAN2.E2end(logid)


if __name__ == "__main__":
    main()
