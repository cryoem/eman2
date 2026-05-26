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

import sys, os

import EMAN2_cppwrap as _cppwrap
_c_emdata_init        = _cppwrap.EMData.__init__
_c_get_image_count    = _cppwrap.EMUtil.get_image_count
_c_read_image         = _cppwrap.EMData.read_image
_c_read_images        = _cppwrap.EMData.read_images

from EMAN3 import EMData, EMUtil
from eman2_gui.emapplication import EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.valslider import ValSlider

EMData.__initc           = _c_emdata_init
EMUtil.get_image_count_c = staticmethod(_c_get_image_count)
EMData.read_image_c      = _c_read_image
EMData.read_images_c     = staticmethod(_c_read_images)

from PyQt5 import QtWidgets, QtCore
import h5py


def hdf_nimg(path):
    with h5py.File(path, 'r') as f:
        return len(f['MDF']['images'])


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
        self.unidirectional = args.unidirectional
        self.avgseq = args.avgseq or 1
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

        self._update_display()

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
        self._print_results()
        QtWidgets.QApplication.quit()

    def keyPressEvent(self, event):
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
                    "start": self.start * a,
                    "end":   self.end   * a}
        else:
            data = {"unidirectional": False,
                    "neg_start": self.neg_start * a,
                    "neg_end":   self.neg_end   * a,
                    "pos_start": self.pos_start * a,
                    "pos_end":   self.pos_end   * a}
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        print(f"Ranges written to {path}")

    def _print_results(self):
        self._write_json()
        a = self.avgseq
        movie = self.args.full_movie or \
                os.path.basename(self.args.thumbnailmovie).replace('_thumb', '')

        if self.unidirectional:
            s, e = self.start, self.end
            print("\n=== Selected ranges (thumbnail frames) ===")
            print(f"  start: {s}  end: {e}  (avgseq={a})")
            print("\n=== Converted to full-stack frames ===")
            print(f"  start: {s*a}  end: {e*a}")
            print("\n=== Command ===")
            print(
                f"e3ContinuousTilt.py \\\n"
                f"    --movie {movie} \\\n"
                f"    --gainfile=GainFiducial.hdf \\\n"
                f"    --tiltRange=120 --compressbits=6 \\\n"
                f"    --startTlt={s*a} --endTlt={e*a}"
            )
        else:
            ns, ne = self.neg_start, self.neg_end
            ps, pe = self.pos_start, self.pos_end
            print("\n=== Selected ranges (thumbnail frames) ===")
            print(f"  neg: {ns} → {ne}  (avgseq={a})")
            print(f"  pos: {ps} → {pe}  (avgseq={a})")
            print("\n=== Converted to full-stack frames ===")
            print(f"  neg: {ns*a} → {ne*a}")
            print(f"  pos: {ps*a} → {pe*a}")
            print("\n=== Command ===")
            print(
                f"e3ContinuousTilt.py \\\n"
                f"    --movie {movie} \\\n"
                f"    --gainfile=GainFiducial.hdf \\\n"
                f"    --tiltRange=120 --compressbits=6 \\\n"
                f"    --startTltNeg={ns*a} --endTltNeg={ne*a} \\\n"
                f"    --startTltPos={ps*a} --endTltPos={pe*a}"
            )


def main():
    import EMAN2
    progname = os.path.basename(sys.argv[0])
    usage = """e3pick_ranges.py thumbnailmovie [options]

Pick good frame ranges from a continuous-tilt thumbnail movie.
The thumbnail movie is a single HDF stack: [neg frames, reversed] + [pos frames].
"""
    parser = EMAN2.EMArgumentParser(usage=usage, version=EMAN2.EMANVERSION)
    parser.add_pos_argument(name="thumbnailmovie", help="Thumbnail HDF stack (neg reversed + pos concatenated).", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", row=0, col=0, rowspan=1, colspan=3)
    parser.add_argument("--full_movie", type=str, default="", help="Full-resolution movie for avgseq inference (optional).", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", filecheck=False, row=1, col=0, rowspan=1, colspan=3)
    parser.add_argument("--avgseq", type=int, default=0, help="Thumbnail averaging factor (0 = infer from full_movie, else 1).", guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
    parser.add_argument("--unidirectional", default=False, help="Single-direction acquisition (no neg/pos split).", action="store_true", guitype='boolbox', row=3, col=0, rowspan=1, colspan=1)
    parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross-platform PPID.")
    (options, args) = parser.parse_args()

    logid = EMAN2.E2init(sys.argv, options.ppid)

    thumbnailmovie = args[0] if args else options.thumbnailmovie
    full_movie = options.full_movie if options.full_movie else None

    avgseq = options.avgseq if options.avgseq else None
    if avgseq is None and full_movie:
        n_thumb = hdf_nimg(thumbnailmovie)
        n_full  = hdf_nimg(full_movie)
        avgseq = round(n_full / n_thumb)
        print(f"  avgseq inferred: {n_full} full / {n_thumb} thumb = {avgseq}")

    options.thumbnailmovie = thumbnailmovie
    options.full_movie = full_movie
    options.avgseq = avgseq or 1

    print(f"Loading {thumbnailmovie} ...")
    frames = load_stack(thumbnailmovie)

    app = EMApp()
    win = PickRangesWindow(frames, options)
    win.show()
    app.exec_()

    EMAN2.E2end(logid)


if __name__ == "__main__":
    main()
