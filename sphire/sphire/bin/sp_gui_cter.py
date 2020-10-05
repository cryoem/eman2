#!/usr/bin/env python
from __future__ import print_function, division
from past.utils import old_div

#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2015-2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2015-2019 Max Planck Institute of Molecular Physiology
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

import EMAN2
import EMAN2_cppwrap
import EMAN2_meta
import PyQt5
import PyQt5.QtCore
import PyQt5.QtWidgets
import collections
import eman2_gui.emapplication
import eman2_gui.emimage2d
import eman2_gui.emplot2d
import eman2_gui.emshape
import numpy
import optparse
import os
import scipy.interpolate
from ..libpy import sp_global_def
from ..libpy import sp_statistics
from ..libpy import sp_utilities
import sys
import traceback
import eman2_gui.valslider

# from numpy import array,arange

# try:
# 	from PyQt5 import QtCore, QtGui, QtOpenGL
# 	from PyQt5.QtCore import Qt
# 	from PyQt5.QtCore import QTimer
# 	from PyQt5.QtWidgets import QListWidget, QLabel, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QSizePolicy, QComboBox, QMessageBox, QFileDialog, QListWidgetItem, qApp
# except:
# except ImportError:
# 	try:
# 		from PyQt4 import QtCore, QtGui, QtOpenGL
# 		from PyQt4.QtCore import Qt
# 		from PyQt4.QtCore import QTimer
# 		from PyQt4.QtGui import QListWidget, QLabel, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QSizePolicy, QComboBox, QMessageBox, QFileDialog, QListWidgetItem, qApp
# 	except:
# 		print("Warning: PyQt4 or 5 must be installed")
# 		sys.exit(1)


"""
Scipy now calls numpy 1.15, which generates numerous warnings of the form 
"RuntimeWarning: numpy.dtype size changed, may indicate binary incompatibility. Expected 96, got 88".
Filterwarnings suppreses this message.
"""
numpy.warnings.filterwarnings("ignore", message="numpy.dtype")


def run():
    progname = os.path.basename(sys.argv[0])
    usage = (
        progname
        + """  cter_ctf_file 
	This GUI application is designed for the evaluation of micrographs using the parameters outputed by CTER.
	"""
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)

    parser.add_option(
        "--ctffind",
        action="store_true",
        default=False,
        help="use CTFFIND outputs (PWROT_DIR/*_avrot.txt, POWER2D_DIR/*.mrc) ",
    )
    parser.add_option(
        "--pwrot_dir",
        default="pwrot",
        help="directory for 1D profiles (default: pwrot)",
    )
    parser.add_option(
        "--power2d_dir",
        default="power2d",
        help="directory for 2D power spectra (default: power2d)",
    )
    parser.add_option(
        "--micthumb_dir",
        default="micthumb",
        help="directory for micrograph thumbnails (default: micthumb)",
    )

    (options, args) = parser.parse_args(sys.argv[1:])
    ####print('options', options)
    ####sys.exit()

    if len(args) > 2:
        print("see usage " + usage)
        sys.exit()

    app = eman2_gui.emapplication.EMApp()

    cter_ctf_file = None
    if len(args) == 1 and args[0] != "":
        cter_ctf_file = args[0]
    # else: # Do nothing

    # Make sure main window is shown, raised, and activated upon startup.
    ####gui = SXGuiCter(use_ctffind=options.ctffind, pwrot_dir=options.pwrot_dir, power2d_dir=options.power2d_dir, micthumb_dir=options.micthumb_dir)
    gui = SXGuiCter(options)
    gui.show()
    gui.raise_()
    gui.activateWindow()

    # read CTER partres file if necessary
    if cter_ctf_file != None:
        gui.readCterPartresFile(os.path.relpath(cter_ctf_file))

    # NOTE: 2016/03/08 Toshio Moriya
    # Unfortunately, EMApp.execute() will print logid (e.g. "None") upon exit...
    app.execute()


class SXListWidget(PyQt5.QtWidgets.QListWidget):
    """Exactly like a normal list widget but intercepts a few keyboard events"""

    def keyPressEvent(self, event):
        if event.key() in (PyQt5.QtCore.Qt.Key_Up, PyQt5.QtCore.Qt.Key_Down):
            PyQt5.QtWidgets.QListWidget.keyPressEvent(self, event)
            return

        self.emit(PyQt5.QtCore.SIGNAL("keypress"), event)


class SXPlot2DWidget(eman2_gui.emplot2d.EMPlot2DWidget):

    mouseup = PyQt5.QtCore.pyqtSignal(PyQt5.QtGui.QMouseEvent)

    def full_refresh(self):
        """
		This function is called from resizeGL and from the inspector when somebody toggles the display of a line
		"""
        self.needupd = 1
        self.del_shapes(
            ("xcross", "ycross", "lcross", "Circle")
        )  # for EMPlot2DInspector
        self.del_shapes(
            ("error_astig", "error_def", "error_ctf")
        )  # for SXGuiCter.wplotrotavgcoarse & SXGuiCter.wplotrotavgfine

    def mouseReleaseEvent(self, event):
        eman2_gui.emplot2d.EMPlot2DWidget.mouseReleaseEvent(self, event)
        if event.button() == PyQt5.QtCore.Qt.LeftButton:
            self.mouseup.emit(event)  # self.emit(QtCore.SIGNAL("mouseup"),event)


class SXLogoButton(PyQt5.QtWidgets.QLabel):
    def __init__(self, imagename, logo_width, parent=None, keep_aspect=False):
        super(SXLogoButton, self).__init__(parent)

        # Width of logo image
        logo_file_path = "{0}{1}".format(EMAN2.get_image_directory(), imagename)

        # Non-square option
        if keep_aspect:
            logo_map = PyQt5.QtGui.QPixmap(logo_file_path)
            logo_height = round(
                old_div(logo_map.height() * logo_width, logo_map.width())
            )
        else:
            logo_height = logo_width

        # Style of widget
        self.setFixedSize(logo_width, logo_height)
        self.customButtonStyle = """
			SXLogoButton {{background-color: rgba(0, 0, 0, 0); border: 0px solid black; border-radius: 0px; image: url("{0}");}}
			""".format(
            logo_file_path
        )

        # Set style and add click event
        self.setStyleSheet(self.customButtonStyle)

    def add_sxmenu_item_btn_widget(self, sxmenu_item_btn_subarea_widget):
        sxmenu_item_btn_subarea_widget.addWidget(self)


class SXCterParam:
    """
	Class containing the following attributes:
		idx_cter : index
		label : text label
		widget: Qt widget
	"""

    def __init__(self, idx_cter, label=None, widget=None):
        self.idx_cter = idx_cter
        self.label = label
        self.widget = widget

    def list_info(self):
        info = []
        info.append(self.idx_cter)
        info.append(self.label)
        info.append(self.widget)

        return info


class SXCterSorted:
    """
	Class containing the following attributes:
		idx_cter : index from SXCterParam
		idx_sort : sorted index (in current class)
		param_name : parameter name
	"""

    def __init__(self, idx_cter, idx_sort, param_name):
        self.idx_cter = idx_cter
        self.idx_sort = idx_sort
        self.param_name = param_name

    def list_info(self):
        info = []
        info.append(self.idx_cter)
        info.append(self.idx_sort)
        info.append(self.param_name)

        return info


class SXHistMap:
    """
	Class containing the following attributes:
		param_name: parameter name
		idx_cter : index from SXCterParam
		idx_sort : index from SXCterSorted
		val_min : minimum value over all entries
		val_max : maximum value over all entries
		unapply_threshold_lower : unapplied lower threshold
		unapply_threshold_upper : unapplied upper threshold
		unapply_widget_lower : Qt widget for unapplied lower threshold
		unapply_widget_upper : Qt widget for unapplied upper hreshold
		apply_threshold_lower : applied lower threshold
		apply_threshold_upper : applied upper threshold
		apply_widget_lower : Qt widget for applied lower threshold
		apply_widget_upper : Qt widget for applied upper threshold
		label : text label
	"""

    def __init__(
        self,
        param_name,
        idx_cter,
        idx_sort,
        val_min,
        val_max,
        unapply_threshold_lower,
        unapply_threshold_upper,
        unapply_widget_lower,
        unapply_widget_upper,
        apply_threshold_lower,
        apply_threshold_upper,
        apply_widget_lower,
        apply_widget_upper,
        label,
    ):
        self.param_name = param_name
        self.idx_cter = idx_cter
        self.idx_sort = idx_sort
        self.val_min = val_min
        self.val_max = val_max
        self.unapply_threshold_lower = unapply_threshold_lower
        self.unapply_threshold_upper = unapply_threshold_upper
        self.unapply_widget_lower = unapply_widget_lower
        self.unapply_widget_upper = unapply_widget_upper
        self.apply_threshold_lower = apply_threshold_lower
        self.apply_threshold_upper = apply_threshold_upper
        self.apply_widget_lower = apply_widget_lower
        self.apply_widget_upper = apply_widget_upper
        self.label = label

    def list_info(self):
        info = []
        info.append(self.param_name)
        info.append(self.idx_cter)
        info.append(self.idx_sort)
        info.append(self.val_min)
        info.append(self.val_max)
        info.append(self.unapply_threshold_lower)
        info.append(self.unapply_threshold_upper)
        info.append(self.unapply_widget_lower)
        info.append(self.unapply_widget_upper)
        info.append(self.apply_threshold_lower)
        info.append(self.apply_threshold_upper)
        info.append(self.apply_widget_lower)
        info.append(self.apply_widget_upper)
        info.append(self.label)

        return info


class SXCheckboxMap:
    """
	Class containing the following attributes:
		item_name : item name
		item_label : text label
		idx_rotinf : index
		item_widget : Qt Widget
	"""

    def __init__(self, item_name, item_label, idx_rotinf, item_widget):
        self.item_name = item_name
        self.item_label = item_label
        self.idx_rotinf = idx_rotinf
        self.item_widget = item_widget

    def list_info():
        info = []
        info.append(self.item_name)
        info.append(self.item_label)
        info.append(self.idx_rotinf)
        info.append(self.item_widget)

        return info


class SXThresholdMap:
    """
	Class containing the following attributes: label, color, index
	"""

    def __init__(self, label, color, index):
        self.label = label
        self.color = color
        self.index = index

    def list_info():
        info = []
        info.append(self.label)
        info.append(self.color)
        info.append(self.index)

        return info


class SXGuiCter(PyQt5.QtWidgets.QWidget):
    ####def __init__(self, use_ctffind=None, pwrot_dir='pwrot', power2d_dir='power2d', micthumb_dir='micthumb'):
    def __init__(self, options):
        """Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,quality)
		'parms' is [box size,ctf,box coord,set of excluded boxnums,quality,oversampling]
		"""

        super(SXGuiCter, self).__init__(None)

        ###		# MRK_TEST: Flags to test experimental functions
        ###		self.is_enable_max_power = False # MRK_TEST: Can this be an option in future?

        # NOTE: 2016/03/15 Toshio Moriya
        # Due to the representation error of float number,
        # default thresholds using max/min value of each parameter can be different after
        # saving/loading the threshold values from the file.
        # This problem is solved by round() function for thresholding parameter values
        #
        # About Precision:
        # Double precision numbers have 53 bits (16 digits ~~ 15.95 digits = 53*log10(2)) of precision and
        # regular floats have 24 bits (8 digits ~~ 7.22 digits = 24*log10(2)) of precision.
        # The floating point in python uses double precision to store the values, and
        # usually taken to be 15 for practical purposes
        #
        self.round_ndigits = 15

        self.set_ctffind(options)
        self.setWindowIcon(
            PyQt5.QtGui.QIcon(EMAN2.get_image_directory() + "sparxicon.png")
        )
        self.installEventFilter(self)  # Necessary for self.eventFilter()

        self.n_idx_cter = self.build_cter_params()  # returns length of list
        self.n_idx_cter_extra = 2
        self.build_cter_params_sort()
        self.n_idx_hist = self.build_hist_params()  # returns length of list
        self.build_rot1d_indices()
        self.build_threshold_control_map()
        self.build_graph_map()
        self.enumerate_threshold_status()
        self.build_threshold_status_labels()

        self.initialize_paths()
        self.initialize_sorting()
        self.initialize_display()
        self.initialize_popup_windows()

        # Emit signals
        self.whistparam.mouseup.connect(
            self.histparammouseup
        )  # self.whistparam.connect(self.whistparam,QtCore.SIGNAL("mouseup"),self.histparammouseup)
        self.wscatterparam.mouseup.connect(
            self.plotparammouseup
        )  # self.wscatterparam.connect(self.wscatterparam,QtCore.SIGNAL("mouseup"),self.plotparammouseup)

        # Manage windows
        self.draw_main_window()
        self.signal_handler()
        self.set_size_popup_windows()

        # Try to recover sizes & positions of windows of the previous GUI session
        EMAN2.E2loadappwin("sp_gui_cter", "main", self)
        EMAN2.E2loadappwin("sp_gui_cter", "fft", self.wfft.qt_parent)
        EMAN2.E2loadappwin("sp_gui_cter", "imgmicthumb", self.wimgmicthumb.qt_parent)
        EMAN2.E2loadappwin(
            "sp_gui_cter", "plotcoarse", self.wplotrotavgcoarse.qt_parent
        )
        EMAN2.E2loadappwin("sp_gui_cter", "plotfine", self.wplotrotavgfine.qt_parent)
        EMAN2.E2loadappwin("sp_gui_cter", "histparam", self.whistparam.qt_parent)
        EMAN2.E2loadappwin("sp_gui_cter", "plotparam", self.wscatterparam.qt_parent)

        # This section is responsible for background updates
        self.checkedpwrot = False
        self.busy = False
        self.needredisp = False

        self.timer = PyQt5.QtCore.QTimer()
        self.timer.timeout.connect(self.timeOut)
        self.timer.start(100)

    def set_ctffind(self, options):
        self.use_ctffind = options.ctffind
        self.pwrot_dir = options.pwrot_dir
        self.power2d_dir = options.power2d_dir
        self.micthumb_dir = options.micthumb_dir

        if not self.use_ctffind:
            self.pwrot_suffix = "_rotinf.txt"
            self.power2d_suffix = "_pws.hdf"
            self.micthumb_suffix = "_thumb.hdf"

        # CTFFIND extensions hardwired by CTFFIND
        else:
            self.pwrot_suffix = "_avrot.txt"  # '_rotinftxt'
            self.power2d_suffix = ".mrc"  # '_pws.hdf'

            # Thumbnails not generated
            self.micthumb_suffix = ".mrc"

    def build_cter_params(self):
        self.cter_params = collections.OrderedDict()

        self.add_cter_param("id", label="CTER ID")
        self.add_cter_param("select", label="Select")
        self.add_cter_param("def", label="Defocus [um]")
        self.add_cter_param("cs", label="Cs [mm]")
        self.add_cter_param("vol", label="Voltage [kV]")
        self.add_cter_param("apix", label="Pixel Size [A]")
        self.add_cter_param("bfactor", label="B-factor [A^2]")
        self.add_cter_param("total_ac", label="Total Amp. Contrast [%]")
        self.add_cter_param("astig_amp", label="Astig. Amp.[um]")
        self.add_cter_param("astig_ang", label="Astig. Ang.[deg]")
        self.add_cter_param("sd_def", label="Defocus SD [um]")
        self.add_cter_param("sd_total_ac", label="Total Amp. Contrast SD [%]")
        self.add_cter_param("sd_astig_amp", label="Astig. Amp. SD [um]")
        self.add_cter_param("sd_astig_ang", label="Astig. Ang. SD [deg]")
        self.add_cter_param("cv_def", label="Defocus CV [%]")
        self.add_cter_param("cv_astig_amp", label="Astig. Amp. CV [%]")
        self.add_cter_param("error_def", label="Defocus Freq. Limit [1/A]")
        self.add_cter_param("error_astig", label="Astig. Freq. Limit [1/A]")
        self.add_cter_param("error_ctf", label="CTF Freq. Limit [1/A]")
        self.add_cter_param("max_freq", label="Max Freq. [A]")
        self.add_cter_param("reserved", label="Reserved")
        self.add_cter_param("const_ac", label="Const. Amp. Contrast [%]")
        self.add_cter_param("phase_shift", label="Phase Shift [deg]")
        self.add_cter_param("mic_name", label="Micrograph")
        ###if self.is_enable_max_power == True:
        ###self.add_cter_param('max_power ', label='PW. Rot. File')

        return len(self.cter_params)

    def add_cter_param(self, param, label=None, widget=None):
        idx_cter = len(self.cter_params)
        self.cter_params[param] = SXCterParam(idx_cter, label, widget)

    def build_cter_params_sort(self):
        # It's allowed to allow sorting by a parameter below, but not computing a histogram, e.g., ID and micrograph.

        self.cter_params_sort = collections.OrderedDict()

        self.add_cter_param_sort("id")  # 0
        self.add_cter_param_sort("mic_name")  # 1
        self.add_cter_param_sort("def")  # 2
        self.add_cter_param_sort("total_ac")  # 3
        self.add_cter_param_sort("astig_amp")  # 4
        self.add_cter_param_sort("astig_ang")  # 5
        self.add_cter_param_sort("sd_def")
        # self.add_cter_param_sort('cv_def')			# 6
        # self.add_cter_param_sort('sd_total_ac')		# 7
        # self.add_cter_param_sort('sd_astig_amp')
        # self.add_cter_param_sort('cv_astig_amp')	# 8
        # self.add_cter_param_sort('sd_astig_ang')	# 9
        self.add_cter_param_sort("error_def")  # 10
        self.add_cter_param_sort("error_astig")  # 11
        # self.add_cter_param_sort('error_ctf')		# 12
        # self.add_cter_param_sort('max_freq')		# 13
        # self.add_cter_param_sort('reserved')		# 14
        self.add_cter_param_sort("phase_shift")  # 15
        ###if self.is_enable_max_power == True:
        ###self.add_cter_param_sort('max_power')

        return len(self.cter_params_sort)

    def add_cter_param_sort(self, param):
        idx_sort = len(self.cter_params_sort)
        self.cter_params_sort[param] = SXCterSorted(
            self.cter_params[param].idx_cter, idx_sort, param
        )

    def build_hist_params(self):
        """
		These are the parameters which will be displayed in the histogram pulldown in the GUI.
		If not listed in cter_params_sort above, a warning will be printed to the screen.
		"""

        self.hist_params = collections.OrderedDict()

        self.add_hist_param("def", 0, 5, 0, 5, None, None, 0, 5, None, None)
        self.add_hist_param("total_ac", 0, 100, 0, 100, None, None, 0, 100, None, None)
        self.add_hist_param("astig_amp", 0, 1, 0, 1, None, None, 0, 1, None, None)
        self.add_hist_param("astig_ang", 0, 180, 0, 180, None, None, 0, 180, None, None)
        self.add_hist_param("sd_def", 0, 5, 0, 5, None, None, 0, 5, None, None)
        # self.add_hist_param('cv_def',       0,   5, 0,   5, None, None, 0,   5, None, None)
        # self.add_hist_param('sd_total_ac',  0, 100, 0, 100, None, None, 0, 100, None, None)
        # self.add_hist_param('sd_astig_amp', 0,   1, 0,   1, None, None, 0,   1, None, None)
        # self.add_hist_param('cv_astig_amp', 0,   1, 0,   1, None, None, 0,   1, None, None)
        # self.add_hist_param('sd_astig_ang', 0, 180, 0, 180, None, None, 0, 180, None, None)
        self.add_hist_param("error_def", 0, 10, 0, 10, None, None, 0, 10, None, None)
        self.add_hist_param("error_astig", 0, 10, 0, 10, None, None, 0, 10, None, None)
        # self.add_hist_param('error_ctf',    0,  10, 0,  10, None, None, 0,  10, None, None)
        # self.add_hist_param('max_freq',     0,  10, 0,  10, None, None, 0,  10, None, None)
        # self.add_hist_param('reserved',     0,  10, 0,  10, None, None, 0,  10, None, None)
        self.add_hist_param(
            "phase_shift", 0, 180, 0, 180, None, None, 0, 180, None, None
        )
        ###		if self.is_enable_max_power == True:
        ###self.add_hist_param('max_power', 0, 99999, 0, 99999, None, None, 0, 99999, None, None)

        return len(self.hist_params)

    def add_hist_param(
        self,
        param,
        val_min,
        val_max,
        unapply_threshold_lower,
        unapply_threshold_upper,
        unapply_widget_lower,
        unapply_widget_upper,
        apply_threshold_lower,
        apply_threshold_upper,
        apply_widget_lower,
        apply_widget_upper,
    ):
        n_idx_hist = len(self.hist_params)
        if not param in self.cter_params_sort:
            print(
                "MRK_DEBUG: parameter %s not listed in cter_params_sort. Skipping..."
                % param
            )
        else:
            self.hist_params[param] = SXHistMap(
                param,
                self.cter_params[param].idx_cter,
                self.cter_params_sort[param].idx_cter,
                val_min,
                val_max,
                unapply_threshold_lower,
                unapply_threshold_upper,
                unapply_widget_lower,
                unapply_widget_upper,
                apply_threshold_lower,
                apply_threshold_upper,
                apply_widget_lower,
                apply_widget_upper,
                self.cter_params[param].label,
            )

    def build_threshold_control_map(self):
        self.threshold_control_map = collections.OrderedDict()
        self.add_threshold_control_map("lower", "Lower (blue)", "blue")
        self.add_threshold_control_map("upper", "Upper (red)", "red")
        self.add_threshold_control_map("edit_only", "Edit Only", "black")

        return len(self.threshold_control_map)

    def add_threshold_control_map(self, entry, label, color):
        index = len(self.threshold_control_map)
        self.threshold_control_map[entry] = SXThresholdMap(label, color, index)

    def build_rot1d_indices(self):
        self.rot1d_indices = collections.OrderedDict()

        if not self.use_ctffind:
            self.add_rot1d_indices("cter_id")
            self.add_rot1d_indices("freq")
            self.add_rot1d_indices("exp_no_astig")
            self.add_rot1d_indices("fit_no_astig")
            self.add_rot1d_indices("exp_with_astig")
            self.add_rot1d_indices("fit_with_astig")
            self.add_rot1d_indices("exp_background")
            self.add_rot1d_indices("fit_envelope")
        else:
            self.add_rot1d_indices("freq")
            self.add_rot1d_indices("exp_no_astig")
            self.add_rot1d_indices("exp_with_astig")
            self.add_rot1d_indices("fit_with_astig")
            self.add_rot1d_indices("ccf_exp_fit")
            self.add_rot1d_indices("sigma_squared")

        return len(self.rot1d_indices)

    def add_rot1d_indices(self, item):
        idx_rotinf = len(self.rot1d_indices)
        self.rot1d_indices[item] = idx_rotinf

    def build_graph_map(self):
        self.graph_map = collections.OrderedDict()

        if not self.use_ctffind:
            self.add_graph_map("exp_no_astig", "Exp. No Astig (Black)")
            self.add_graph_map("fit_no_astig", "Fit. No Astig (Blue)")
            self.add_graph_map("exp_with_astig", "Exp. with Astig (Red)")
            self.add_graph_map("fit_with_astig", "Fit. with Astig (Green)")
            self.add_graph_map("exp_background", "Exp. Enhanced (Olive)")
            self.add_graph_map("fit_envelope", "Fit. Enhanced (Cyan)")
        else:
            self.add_graph_map("exp_with_astig", "Exp. with Astig (Black)")
            self.add_graph_map("fit_with_astig", "Fit. with Astig (Blue)")
            self.add_graph_map("ccf_exp_fit", "Exp. CC. Fitted (Red)")
            self.add_graph_map("sigma_squared", "Sigma*2 Noise (Green)")

        return len(self.graph_map)

    def add_graph_map(self, item_name, item_label, item_widget=None):
        idx_rotinf = self.rot1d_indices[item_name]
        self.graph_map[item_name] = SXCheckboxMap(
            item_name, item_label, idx_rotinf, item_widget
        )

    def enumerate_threshold_status(self):
        # Define enumerators for threshold apply status
        i_enum = -1
        i_enum += 1
        self.idx_thresholdset_unapplied = i_enum
        i_enum += 1
        self.idx_thresholdset_applied = i_enum
        i_enum += 1
        self.n_idx_thresholdset = i_enum

    def build_threshold_status_labels(self):
        self.threshold_status_labels = collections.OrderedDict()
        self.threshold_status_labels["unapplied"] = "Unapplied"
        self.threshold_status_labels["applied"] = "Applied"

        return len(self.threshold_status_labels)

    def initialize_paths(self):
        self.cter_partres_file_path = None
        self.cter_entry_list = None
        self.cter_mic_file_path = None
        self.cter_micthumb_file_path = None
        self.cter_pwrot_file_path = None
        self.cter_fft_file_path = None

    def initialize_sorting(self):
        self.curentry = None
        self.cursortidx = 0
        self.cursortorder = False
        self.cursortselect = False
        self.curhistidx = 0
        self.curthresholdcontrol = 0
        self.curentryperbin = 10
        self.cursyncsort = False
        self.curthresholdset = 0

    def initialize_display(self):
        self.curplotrotavgdisplay = False  # True  # (now off by default)
        self.curplotrotzoomdisplay = True
        self.curimgmicthumbdisplay = False  # True  # (now off by default)
        # self.curhistdisable = False
        self.curhistogramdisplay = True
        self.curscatterdisplay = True
        self.curplotfixscale = (
            1.1
        )  # 5  (applied envelope and subtracted background -- can still override from GUI)
        self.curfftdisplay = False

    def initialize_popup_windows(self):
        # NOTE: 2016/03/09 Toshio Moriya
        # To set window flags of EMGLWidget (SXPlot2DWidget and EMImage2DWidget) window,
        # we have to go through its qt_parent attribute to call setWindowTitle()...
        #
        self.wfft = eman2_gui.emimage2d.EMImage2DWidget()
        self.wfft.setWindowTitle("sp_gui_cter - 2D FFT")
        # self.wfft.mmode="app"  # NOT USED?
        self.wfft.qt_parent.setWindowFlags(
            (self.wfft.qt_parent.windowFlags() | PyQt5.QtCore.Qt.CustomizeWindowHint)
            & ~PyQt5.QtCore.Qt.WindowMinimizeButtonHint
        )  # Disabled minimize icon button in window title bar
        self.is_wfft_minimized = False

        self.wimgmicthumb = eman2_gui.emimage2d.EMImage2DWidget()
        self.wimgmicthumb.setWindowTitle("sp_gui_cter - Micrograph Thumbnail")
        # self.wimgmicthumb.mmode="app"  # NOT USED?
        self.wimgmicthumb.qt_parent.setWindowFlags(
            (
                self.wimgmicthumb.qt_parent.windowFlags()
                | PyQt5.QtCore.Qt.CustomizeWindowHint
            )
            & ~PyQt5.QtCore.Qt.WindowMinimizeButtonHint
        )  # Disabled minimize icon button in window title bar
        self.is_wimgmicthumb_minimized = False

        self.wplotrotavgcoarse = SXPlot2DWidget()
        self.wplotrotavgcoarse.setWindowTitle("sp_gui_cter - Plot")
        self.wplotrotavgcoarse.qt_parent.setWindowFlags(
            (
                self.wplotrotavgcoarse.qt_parent.windowFlags()
                | PyQt5.QtCore.Qt.CustomizeWindowHint
            )
            & ~PyQt5.QtCore.Qt.WindowMinimizeButtonHint
        )  # Disabled minimize icon button in window title bar
        self.is_wplotrotavgcoarse_minimized = False

        self.wplotrotavgfine = SXPlot2DWidget()
        self.wplotrotavgfine.setWindowTitle("sp_gui_cter - Plot Zoom")
        self.wplotrotavgfine.qt_parent.setWindowFlags(
            (
                self.wplotrotavgfine.qt_parent.windowFlags()
                | PyQt5.QtCore.Qt.CustomizeWindowHint
            )
            & ~PyQt5.QtCore.Qt.WindowMinimizeButtonHint
        )  # Disabled minimize icon button in window title bar
        self.is_wplotrotavgfine_minimized = False

        self.whistparam = SXPlot2DWidget()
        self.whistparam.setWindowTitle("sp_gui_cter - Histogram")
        self.whistparam.qt_parent.setWindowFlags(
            (
                self.whistparam.qt_parent.windowFlags()
                | PyQt5.QtCore.Qt.CustomizeWindowHint
            )
            & ~PyQt5.QtCore.Qt.WindowMinimizeButtonHint
        )  # Disabled minimize icon button in window title bar
        self.is_whistparam_minimized = False

        self.wscatterparam = SXPlot2DWidget()
        self.wscatterparam.setWindowTitle("sp_gui_cter - Sort Plot")
        self.wscatterparam.qt_parent.setWindowFlags(
            (
                self.wscatterparam.qt_parent.windowFlags()
                | PyQt5.QtCore.Qt.CustomizeWindowHint
            )
            & ~PyQt5.QtCore.Qt.WindowMinimizeButtonHint
        )  # Disabled minimize icon button in window title bar
        self.is_wscatterparam_minimized = False

    def draw_main_window(self):
        """
         _______________________________________________________
		| main                                                 |
		| ______________ ________________ ____________________ |
		| | leftcolumn | | secondcolumn | | biglayout:       | |
		| |            | |              | | ________________ | |
		| |            | |              | | | settings     | | |
		| |            | |              | | |--------------| | |
		| |            | |              | | ________________ | |
		| |            | |              | | | chart        | | |
		| |            | |              | | | ____________ | | |
		| |            | |              | | | | upper    | | | |
		| |            | |              | | | |----------| | | |
		| |            | |              | | | ____________ | | |
		| |            | |              | | | | thresh   | | | |
		| |            | |              | | | |----------| | | |
		| |            | |              | | | ____________ | | |
		| |            | |              | | | | last2row | | | |
		| |            | |              | | | |----------| | | |
		| |            | |              | | |--------------| | |
		| |------------| |--------------| |------------------| |
		|------------------------------------------------------|
		"""

        # Place layout inside QWidget
        templayout = PyQt5.QtWidgets.QHBoxLayout(self)
        templayout.setContentsMargins(0, 0, 0, 0)
        mainwidget = PyQt5.QtWidgets.QWidget(self)
        mainwidget.setObjectName("MainWidgetObject")

        # Color scheme
        background_image_file_path = "{0}sxgui_background.png".format(
            EMAN2.get_image_directory()
        )
        mainwidget.setStyleSheet(
            "QWidget#MainWidgetObject {{background-image: url('{0}')}}".format(
                background_image_file_path
            )
        )
        mainlayout = PyQt5.QtWidgets.QHBoxLayout(mainwidget)
        mainlayout.setContentsMargins(12, 12, 12, 12)
        templayout.addWidget(mainwidget)

        # --------------------------------------------------------------------------------
        # Columns 1-3
        # --------------------------------------------------------------------------------

        editwidth = 100
        sublabelwidth = 140

        leftwidget = PyQt5.QtWidgets.QWidget(self)
        leftwidget.setObjectName("LeftWidgetObject")
        leftwidget.setStyleSheet(
            "QWidget#LeftWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}"
        )

        leftcolumn = PyQt5.QtWidgets.QVBoxLayout(leftwidget)
        leftcolumn.setContentsMargins(0, 10, 0, 10)
        mainlayout.addWidget(leftwidget)

        self.pbopencter = PyQt5.QtWidgets.QPushButton("Open CTER partres file")
        leftcolumn.addWidget(self.pbopencter)

        self.add_centered_label("", leftcolumn)  # spacer
        self.add_centered_label("<b>Display Windows:</b>", leftcolumn)

        self.cbrotavgdisplay = eman2_gui.valslider.CheckBox(
            None, None, self.curplotrotavgdisplay
        )
        self.add_label_with_checkbox(
            "Rot. Avg. Plot", self.cbrotavgdisplay, leftcolumn, labelwidth=sublabelwidth
        )

        self.cbrotzoomdisplay = eman2_gui.valslider.CheckBox(
            None, None, self.curplotrotzoomdisplay
        )
        self.add_label_with_checkbox(
            "Rot. Avg. Plot Zoom",
            self.cbrotzoomdisplay,
            leftcolumn,
            labelwidth=sublabelwidth,
        )

        self.cbhistogramdisplay = eman2_gui.valslider.CheckBox(
            None, None, self.curhistogramdisplay
        )
        self.add_label_with_checkbox(
            "Histogram", self.cbhistogramdisplay, leftcolumn, labelwidth=sublabelwidth
        )

        self.cbscatterdisplay = eman2_gui.valslider.CheckBox(
            None, None, self.curscatterdisplay
        )
        self.add_label_with_checkbox(
            "Sort Plot", self.cbscatterdisplay, leftcolumn, labelwidth=sublabelwidth
        )

        self.cbmicthumbdisplay = eman2_gui.valslider.CheckBox(
            None, None, self.curimgmicthumbdisplay
        )
        self.add_label_with_checkbox(
            "Micrograph Thumbnail",
            self.cbmicthumbdisplay,
            leftcolumn,
            labelwidth=sublabelwidth,
        )

        self.cbfftdisplay = eman2_gui.valslider.CheckBox(None, None, self.curfftdisplay)
        self.add_label_with_checkbox(
            "2D Power Spectrum", self.cbfftdisplay, leftcolumn, labelwidth=sublabelwidth
        )

        self.add_centered_label("", leftcolumn)  # spacer
        self.add_centered_label("<b>Display Curves:</b>", leftcolumn)

        # Add checkboxes for each plot
        for item in self.graph_map:
            self.graph_map[item].item_widget = eman2_gui.valslider.CheckBox(
                None, None, True
            )
            self.add_label_with_checkbox(
                self.graph_map[item].item_label,
                self.graph_map[item].item_widget,
                leftcolumn,
                labelwidth=sublabelwidth,
            )

        self.vbplotfixscale = eman2_gui.valslider.ValBox(
            self, (0, 99999), None, self.curplotfixscale
        )  # default <- self.curplotfixscale
        self.add_label_with_value(
            "Plot Fix Scale",
            self.vbplotfixscale,
            leftcolumn,
            labelwidth=sublabelwidth,
            enabled=True,
            intonly=False,
            style_sheet="color: rgb(0,0,0);",
        )

        self.add_centered_label("", leftcolumn)  # spacer

        self.pbrefreshgraphs = PyQt5.QtWidgets.QPushButton("Refresh Graphs")
        self.pbrefreshgraphs.setEnabled(False)
        leftcolumn.addWidget(self.pbrefreshgraphs)

        leftcolumn.addStretch(1)

        logo = SXLogoButton(
            "logo_transphire_thick.png", 256, parent=self, keep_aspect=True
        )
        logo.add_sxmenu_item_btn_widget(leftcolumn)

        # --------------------------------------------------------------------------------
        # 2nd column
        # --------------------------------------------------------------------------------

        secondcolumn = PyQt5.QtWidgets.QVBoxLayout()
        secondcolumn.setContentsMargins(0, 0, 0, 0)
        mainlayout.addLayout(secondcolumn)

        # plot list and plot mode combobox
        row_span_entry_list = 27  # length of file list (QListWidget)
        self.lbentry = SXListWidget(self)
        self.lbentry.setSizePolicy(
            PyQt5.QtWidgets.QSizePolicy.Preferred, PyQt5.QtWidgets.QSizePolicy.Expanding
        )
        self.lbentry.setMinimumWidth(220)
        secondcolumn.addWidget(self.lbentry)

        labelwidth = 180
        editwidth = 100
        sublabelwidth = editwidth
        sidemargin = 10
        topmargin = 3
        bottommargin = 3
        borderwidth = 1

        # --------------------------------------------------------------------------------
        # The main big layout (VBox): I'll split up into 2 HBoxLayouts: settings & chart
        # --------------------------------------------------------------------------------

        biglayout = PyQt5.QtWidgets.QVBoxLayout()
        biglayout.setContentsMargins(0, 0, 0, 0)
        mainlayout.addLayout(biglayout)

        # ---------------------------------
        # Settings layout
        # ---------------------------------

        selectEMwidget = PyQt5.QtWidgets.QWidget(self)
        selectEMwidget.setObjectName("SelectWidgetObject")
        selectEMwidget.setStyleSheet(
            "QWidget#SelectWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}"
        )

        selectEMlayout = PyQt5.QtWidgets.QHBoxLayout(selectEMwidget)
        biglayout.addWidget(selectEMwidget)

        # There will be two VBox layouts: select and em
        selectlayout = PyQt5.QtWidgets.QVBoxLayout()
        selectlayout.setContentsMargins(0, 0, 0, 0)
        selectEMlayout.addLayout(selectlayout)

        labelwidth = 150

        self.add_centered_label("<b>Selection Summary:</b>", selectlayout)

        self.vbnentry = eman2_gui.valslider.ValBox(self, (0, 10000), None, 0)
        self.add_label_with_value(
            "Num. of entries",
            self.vbnentry,
            selectlayout,
            labelwidth=labelwidth,
            style_sheet="color: rgb(0,0,0);",
        )

        self.vbuncheckcounts = eman2_gui.valslider.ValBox(self, (0, 1000000), None, 0)
        self.add_label_with_value(
            "Unchecked",
            self.vbuncheckcounts,
            selectlayout,
            labelwidth=labelwidth,
            style_sheet="color: rgb(0,0,0);",
        )

        self.vbuncheckratio = eman2_gui.valslider.ValBox(self, (0, 1.0), None, 0)
        self.add_label_with_value(
            "Ratio",
            self.vbuncheckratio,
            selectlayout,
            labelwidth=labelwidth,
            style_sheet="color: rgb(0,0,0);",
            intonly=False,
        )

        # Electron microscopy settings
        labelwidth = 110

        emlayout = PyQt5.QtWidgets.QVBoxLayout()
        emlayout.setContentsMargins(
            0, 0, 0, 0
        )  # I don't know why these are needed to align the layouts in the HBox
        selectEMlayout.addLayout(emlayout)

        self.add_centered_label("<b>Electron Microscopy:</b>", emlayout)

        # Voltage
        self.cter_params["vol"].widget = self.add_label_with_param(
            self.cter_params["vol"].label,
            emlayout,
            0,
            500,
            intonly=False,
            style_sheet="color: rgb(127,127,127);",
            labelwidth=labelwidth,
        )

        # Spherical aberration
        self.cter_params["cs"].widget = self.add_label_with_param(
            self.cter_params["cs"].label,
            emlayout,
            0,
            5,
            intonly=False,
            style_sheet="color: rgb(127,127,127);",
            labelwidth=labelwidth,
        )

        # Pixel size
        self.cter_params["apix"].widget = self.add_label_with_param(
            self.cter_params["apix"].label,
            emlayout,
            0,
            500,
            intonly=False,
            style_sheet="color: rgb(127,127,127);",
            labelwidth=labelwidth,
        )

        # Add logo
        cterlogolayout = PyQt5.QtWidgets.QVBoxLayout()
        cterlogolayout.setContentsMargins(0, 0, 0, 0)
        selectEMlayout.addLayout(cterlogolayout)
        ##logo = SXLogoButton("sxgui_pictograph_cryolo.png", 80, parent=self)
        ##logo.add_sxmenu_item_btn_widget(cterlogolayout)
        logo = SXLogoButton("sxgui_pictograph_cter.png", 80, parent=self)
        logo.add_sxmenu_item_btn_widget(cterlogolayout)

        # --------------------------------------------------------------------------
        # Chart layout will be split up into 3 VBoxLayouts: upper, thresh, last2row
        # --------------------------------------------------------------------------
        labelwidth = 160

        cterwidget = PyQt5.QtWidgets.QWidget(self)
        cterwidget.setObjectName("CterWidgetObject")
        cterwidget.setStyleSheet(
            "QWidget#CterWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}"
        )
        biglayout.addWidget(cterwidget)

        chartlayout = PyQt5.QtWidgets.QVBoxLayout(cterwidget)
        chartlayout.setContentsMargins(0, 10, 10, 10)

        upperlayout = PyQt5.QtWidgets.QHBoxLayout()
        upperlayout.setContentsMargins(0, 0, 0, 0)
        chartlayout.addLayout(upperlayout)

        first2rowslayout = PyQt5.QtWidgets.QVBoxLayout()
        first2rowslayout.setContentsMargins(0, 0, 0, 0)
        upperlayout.addLayout(first2rowslayout)

        self.add_centered_label(
            "<b>Current Entry Info:</b>", first2rowslayout, labelwidth=350
        )  # hardwired labelwidth

        self.ssortedid = eman2_gui.valslider.ValBox(self, (0, 10000), None, 0)
        self.add_label_with_value(
            "Sorted ID",
            self.ssortedid,
            first2rowslayout,
            style_sheet="color: rgb(0,0,0);",
            labelwidth=labelwidth,
            editwidth=editwidth,
        )

        self.cter_params["id"].widget = self.add_label_with_param(
            self.cter_params["id"].label,
            first2rowslayout,
            0,
            10000,
            intonly=True,
            style_sheet="color: rgb(127,127,127);",
            labelwidth=labelwidth,
            editwidth=80,
            maxwidth=100,
            enabled=False,
        )

        # Layout for thresholds and column labels
        threshlayout = PyQt5.QtWidgets.QHBoxLayout()
        threshlayout.maximumSize()
        chartlayout.addLayout(threshlayout)

        # Draw borderless box to preserve spacing
        cterFrame = PyQt5.QtWidgets.QWidget(self)
        cterFrame.setContentsMargins(0, 0, 0, 0)
        cterFrame.setStyleSheet("border: %spx solid transparent;" % borderwidth)
        threshlayout.addWidget(cterFrame)

        # Layout for CTER columns: select through phase shift
        cterlayout = PyQt5.QtWidgets.QVBoxLayout(cterFrame)
        cterlayout.setContentsMargins(0, 0, 0, bottommargin)

        # Selection flag
        self.cter_params["select"].widget = self.add_label_with_param(
            self.cter_params["select"].label,
            cterlayout,
            0,
            1,
            intonly=True,
            labelwidth=labelwidth,
            editwidth=editwidth,
            style_sheet="color: rgb(127,127,127);",
            maxwidth=100,
            enabled=False,
        )

        # Draw boxes around limits
        unappliedLimitFrame = PyQt5.QtWidgets.QWidget(
            self
        )  # QGridLayout doesn't have borders, so I'm enclosing it inside a QFrame, which can
        unappliedLimitFrame.setContentsMargins(0, 0, 0, 0)
        unappliedLimitFrame.setStyleSheet(
            "border: %spx solid rgb(0,0,0);" % borderwidth
        )
        threshlayout.addWidget(unappliedLimitFrame)

        # Layout for unapplied threshholds
        unappliedlayout = PyQt5.QtWidgets.QVBoxLayout(unappliedLimitFrame)
        unappliedlayout.setObjectName("unappliedlayout")
        unappliedlayout.setContentsMargins(
            sidemargin, topmargin, sidemargin, bottommargin
        )

        self.add_centered_label(
            "<b>Unapplied Thresholds:</b>", unappliedlayout, style_sheet="border: 0px;"
        )

        # Applied limits
        appliedLimitFrame = PyQt5.QtWidgets.QWidget(
            self
        )  # QGridLayout doesn't have borders, so I'm enclosing it inside a QFrame, which can
        appliedLimitFrame.setStyleSheet("border: %spx solid rgb(0,0,0);" % borderwidth)
        threshlayout.addWidget(appliedLimitFrame)

        # Layout for applied threshholds
        appliedlayout = PyQt5.QtWidgets.QVBoxLayout(appliedLimitFrame)
        appliedlayout.setObjectName("appliedlayout")
        appliedlayout.setContentsMargins(
            sidemargin, topmargin, sidemargin, bottommargin
        )

        self.add_centered_label(
            "<b>Applied Thresholds:</b>", appliedlayout, style_sheet="border: 0px;"
        )

        # Draw limits, initially min & max
        for idx, (param, curr_hist_params) in enumerate(self.hist_params.items()):
            self.cter_params[param].widget = self.add_label_with_param(
                self.cter_params[param].label,
                cterlayout,
                curr_hist_params.val_min,
                curr_hist_params.val_max,
                intonly=False,
                labelwidth=labelwidth,
                editwidth=editwidth,
                maxwidth=100,
            )
            self.add_threshold_widgets(
                param, unappliedlayout, appliedlayout, labelwidth=80, editwidth=80
            )

        # Last two CTER rows
        last2rowFrame = PyQt5.QtWidgets.QWidget(self)
        last2rowFrame.setContentsMargins(0, 0, 0, 0)
        last2rowFrame.setStyleSheet("border: %spx solid transparent;" % borderwidth)
        chartlayout.addWidget(last2rowFrame)
        last2rowlayout = PyQt5.QtWidgets.QVBoxLayout(last2rowFrame)
        last2rowlayout.setContentsMargins(0, 0, 0, 0)

        # Amplitude contrast and B-factor (These used to be grayed out.)
        self.cter_params["const_ac"].widget = self.add_label_with_param(
            self.cter_params["const_ac"].label,
            last2rowlayout,
            0,
            1600,
            labelwidth=labelwidth,
            editwidth=editwidth,
            enabled=False,
            style_sheet="color: rgb(127,127,127);",
            intonly=False,
            maxwidth=100,
        )
        self.cter_params["bfactor"].widget = self.add_label_with_param(
            self.cter_params["bfactor"].label,
            last2rowlayout,
            0,
            1600,
            labelwidth=labelwidth,
            editwidth=editwidth,
            enabled=False,
            style_sheet="color: rgb(127,127,127);",
            maxwidth=100,
        )

        # ---------------------------------
        # Settings layout
        # ---------------------------------

        settingswidget = PyQt5.QtWidgets.QWidget(self)
        settingswidget.setObjectName("SettingsWidgetObject")
        settingswidget.setStyleSheet(
            "QWidget#SettingsWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}"
        )

        settingslayout = PyQt5.QtWidgets.QHBoxLayout(settingswidget)
        biglayout.addWidget(settingswidget)

        # There will be three VBox layouts: sort, histogram/plot, and save/load
        sortlayout = PyQt5.QtWidgets.QVBoxLayout()
        sortlayout.setContentsMargins(0, 0, 0, 0)
        settingslayout.addLayout(sortlayout)

        self.add_centered_label("<b>Sort CTER Partres Entries:</b>", sortlayout)

        # ---------------------------------
        # Pulldown menu options for sorting
        # ---------------------------------

        self.ssort = PyQt5.QtWidgets.QComboBox(self)
        for param in self.cter_params_sort.keys():
            self.ssort.addItem(self.cter_params[param].label)
        self.ssort.setCurrentIndex(self.cursortidx)
        sortlayout.addWidget(self.ssort)

        # Checkbox to reverse order
        self.cbsortorder = eman2_gui.valslider.CheckBox(None, None, self.cursortorder)
        self.add_label_with_checkbox(
            "Descending", self.cbsortorder, sortlayout, labelwidth=sublabelwidth
        )

        # Checkbox to list deselected files first
        self.cbsortselect = eman2_gui.valslider.CheckBox(None, None, self.cursortselect)
        self.add_label_with_checkbox(
            "Sort Select", self.cbsortselect, sortlayout, labelwidth=sublabelwidth
        )

        sortlayout.addStretch(1)

        self.pbreapplysort = PyQt5.QtWidgets.QPushButton("Reapply Sort")
        self.pbreapplysort.setEnabled(False)
        sortlayout.addWidget(self.pbreapplysort)

        # ---------------------------------
        # Histogram & Plot Settings
        # ---------------------------------

        histplotlayout = PyQt5.QtWidgets.QVBoxLayout()
        histplotlayout.setContentsMargins(
            0, 0, 0, 0
        )  # I don't know why these are needed to align the layouts in the HBox
        settingslayout.addLayout(histplotlayout)

        self.add_centered_label("<b>Histogram & Plot Settings:</b>", histplotlayout)

        # Pulldown menu options for histogram
        self.shist = PyQt5.QtWidgets.QComboBox(self)
        self.shist.setMaximumWidth(250)
        for param in self.hist_params.keys():
            self.shist.addItem(self.cter_params[param].label)
        self.shist.setCurrentIndex(self.curhistidx)
        histplotlayout.addWidget(self.shist)

        # Pulldown menu to move lower/upper threshold
        self.add_label_with_pulldown(
            "Move Threshold", histplotlayout, labelwidth=100, menuwidth=140
        )

        # Checkbox to sort according to histogrammed/plotted values
        self.cbsyncsort = eman2_gui.valslider.CheckBox(None, None, self.cursyncsort)
        self.add_label_with_checkbox(
            "Sync. w/Sort", self.cbsyncsort, histplotlayout, labelwidth=95
        )

        self.vsentryperbin = self.add_label_with_param(
            "counts/bin", histplotlayout, 0, 10000, labelwidth=95, intonly=True
        )

        histplotlayout.addStretch(1)

        self.pbapplyallthreshold = PyQt5.QtWidgets.QPushButton("Apply All Thresholds")
        self.pbapplyallthreshold.setMaximumWidth(250)
        self.pbapplyallthreshold.setEnabled(False)
        histplotlayout.addWidget(self.pbapplyallthreshold)

        # ---------------------------------
        # Save/Load thresholds
        # ---------------------------------

        saveloadlayout = PyQt5.QtWidgets.QVBoxLayout()
        saveloadlayout.setContentsMargins(0, 0, 0, 0)
        saveloadlayout.setSpacing(8.5)
        settingslayout.addLayout(saveloadlayout)

        self.add_centered_label("<b>Save/Load Thresholds:</b>", saveloadlayout)

        # Pulldown menu to save unapplied/applied thresholds
        self.sthresholdset = PyQt5.QtWidgets.QComboBox(self)
        for status in self.threshold_status_labels:
            self.sthresholdset.addItem(self.threshold_status_labels[status])
        self.sthresholdset.setCurrentIndex(self.curthresholdset)
        saveloadlayout.addWidget(self.sthresholdset)

        # Save/Load threshold buttons
        self.pbsavethresholdset = PyQt5.QtWidgets.QPushButton("Save")
        self.pbloadthresholdset = PyQt5.QtWidgets.QPushButton("Load")
        self.add_two_buttons(
            self.pbsavethresholdset, self.pbloadthresholdset, saveloadlayout
        )

        self.add_centered_label("<b>Save Selection:</b>", saveloadlayout)

        # Prefix for output files
        self.vfileprefix = eman2_gui.valslider.StringBox(self, None, "Trial00")
        self.add_label_with_textbox("File Prefix", self.vfileprefix, saveloadlayout)

        saveloadlayout.addStretch(1)

        self.pbsaveselection = PyQt5.QtWidgets.QPushButton("Save Selection")
        self.pbsaveselection.setEnabled(False)
        saveloadlayout.addWidget(self.pbsaveselection)

        # Pad at the bottom
        biglayout.addStretch(1)

        self.setWindowTitle("sp_gui_cter - Control Panel")

    def signal_handler(self):
        self.pbopencter.clicked[bool].connect(self.openCterPartres)

        self.cbrotavgdisplay.valueChanged.connect(self.newRotAvgDisplay)
        self.cbrotzoomdisplay.valueChanged.connect(self.newRotZoomDisplay)
        self.cbmicthumbdisplay.valueChanged.connect(self.newMicThumbDisplay)
        self.cbhistogramdisplay.valueChanged.connect(self.newHistogramDisplay)
        self.cbscatterdisplay.valueChanged.connect(self.newScatterDisplay)
        self.cbfftdisplay.valueChanged.connect(self.newFFTDisplay)

        for item in self.graph_map:
            self.graph_map[item].item_widget.valueChanged.connect(
                self.updatePlotVisibility
            )

        self.vbplotfixscale.valueChanged.connect(self.newPlotFixScale)
        self.pbrefreshgraphs.clicked[bool].connect(self.refreshGraphs)

        self.lbentry.currentRowChanged[int].connect(self.newEntry)
        self.lbentry.itemChanged.connect(self.updateEntrySelect)

        self.ssort.currentIndexChanged[int].connect(self.newSort)
        self.cbsortorder.valueChanged.connect(self.newSortOrder)
        self.cbsortselect.valueChanged.connect(self.newSortSelect)
        self.pbreapplysort.clicked[bool].connect(self.reapplySort)

        for (param, curr_hist_params) in self.hist_params.items():
            curr_hist_params.unapply_widget_lower.valueChanged.connect(
                self.newThresholdLower
            )
            curr_hist_params.unapply_widget_upper.valueChanged.connect(
                self.newThresholdUpper
            )

        self.shist.currentIndexChanged[int].connect(self.newHistogramRow)
        self.sthresholdcontrol.currentIndexChanged[int].connect(
            self.newThresholdControl
        )
        self.cbsyncsort.valueChanged.connect(self.newSyncSort)
        self.vsentryperbin.valueChanged.connect(self.newEntryPerBin)
        self.pbapplyallthreshold.clicked[bool].connect(self.applyAllThresholds)

        self.sthresholdset.currentIndexChanged[int].connect(self.newThresholdSet)
        self.pbsavethresholdset.clicked[bool].connect(self.saveThresholdSet)
        self.pbloadthresholdset.clicked[bool].connect(self.loadThresholdSet)

        self.pbsaveselection.clicked[bool].connect(self.saveSelection)

    def set_size_popup_windows(self):
        # Set default sizes & positions of windows in case this is the first time to run in this project directory
        # (I figured these values out by printing the width and height in resize event)
        main_win_width = 1166
        main_win_height = 726
        child_win_width = 1166
        child_win_height = 726
        img_size = 512

        win_left = 0
        win_top = 0
        win_left_shift = 30
        win_top_shift = 30

        self.resize(main_win_width, main_win_height)
        self.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift
        self.wscatterparam.qt_parent.resize(child_win_width, child_win_height)
        self.wscatterparam.qt_parent.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift
        self.whistparam.qt_parent.resize(child_win_width, child_win_height)
        self.whistparam.qt_parent.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift
        self.wplotrotavgfine.qt_parent.resize(child_win_width, child_win_height)
        self.wplotrotavgfine.qt_parent.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift
        self.wplotrotavgcoarse.qt_parent.resize(child_win_width, child_win_height)
        self.wplotrotavgcoarse.qt_parent.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift

        self.wimgmicthumb.set_data(
            sp_utilities.model_blank(img_size, img_size, bckg=1.0)
        )  # resize does not work if no image is set, img_size appears to affect centering

        self.wimgmicthumb.qt_parent.resize(child_win_width, child_win_height)
        self.wimgmicthumb.scroll_to(-1 * img_size, -1 * img_size)
        self.wimgmicthumb.qt_parent.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift

        self.wfft.set_data(
            sp_utilities.model_blank(img_size, img_size, bckg=1.0)
        )  # resize does not work if no image is set
        self.wfft.qt_parent.resize(child_win_width, child_win_height)
        self.wfft.qt_parent.move(win_left, win_top)
        win_left += win_left_shift
        win_top += win_top_shift
        self.wfft.scroll_to(-1 * img_size, -1 * img_size)

    def add_centered_label(self, labeltext, target, labelwidth=None, style_sheet=None):
        temp_label = PyQt5.QtWidgets.QLabel(labeltext, self)
        temp_label.setAlignment(
            PyQt5.QtCore.Qt.AlignHCenter | PyQt5.QtCore.Qt.AlignVCenter
        )
        if labelwidth:
            temp_label.setMaximumWidth(labelwidth)
        if style_sheet:
            temp_label.setStyleSheet(style_sheet)  # "border: 0px;")  #
        else:
            temp_label.setStyleSheet("border: 0px;")  # style_sheet)  #
        target.addWidget(temp_label)

    def add_label_with_value(
        self,
        labeltext,
        valueentry,
        target,
        labelwidth=90,
        enabled=False,
        intonly=True,
        style_sheet="color: rgb(127,127,127);",
        editwidth=80,
        maxwidth=100,
    ):

        # Label
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)
        temp_label = PyQt5.QtWidgets.QLabel(labeltext, self)
        temp_label.setAlignment(
            PyQt5.QtCore.Qt.AlignRight | PyQt5.QtCore.Qt.AlignVCenter
        )
        temp_label.setFixedSize(PyQt5.QtCore.QSize(labelwidth, 20))
        temp_hbox.addWidget(temp_label)

        # Value
        valueentry.setEnabled(enabled)
        valueentry.intonly = intonly
        valueentry.text.setStyleSheet(style_sheet)
        valueentry.text.setMinimumSize(PyQt5.QtCore.QSize(editwidth, 0))
        valueentry.text.setMaximumWidth(maxwidth)
        temp_hbox.addWidget(valueentry)
        temp_hbox.addStretch(1)
        target.addLayout(temp_hbox, 0)

    def add_label_with_param(
        self,
        labeltext,
        target,
        val_min,
        val_max,
        intonly=False,
        style_sheet="color: rgb(0,0,0);",
        labelwidth=80,
        editwidth=80,
        maxwidth=100,
        enabled=False,
    ):

        # Initialize ValBox
        val_default = val_min
        valueentry = eman2_gui.valslider.ValBox(
            self, (val_min, val_max), None, val_default
        )

        self.add_label_with_value(
            labeltext,
            valueentry,
            target,
            labelwidth=labelwidth,
            enabled=enabled,
            intonly=intonly,
            style_sheet=style_sheet,
            editwidth=editwidth,
            maxwidth=maxwidth,
        )

        return valueentry  # returns a ValBox

    def add_label_with_checkbox(self, labeltext, checkbutton, target, labelwidth=140):
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)
        temp_label = PyQt5.QtWidgets.QLabel(labeltext, self)
        temp_label.setAlignment(
            PyQt5.QtCore.Qt.AlignRight | PyQt5.QtCore.Qt.AlignVCenter
        )
        temp_label.setFixedSize(PyQt5.QtCore.QSize(labelwidth, 20))
        temp_hbox.addWidget(temp_label)

        temp_hbox.addWidget(checkbutton)
        target.addLayout(temp_hbox, 0)

    def add_threshold_widgets(
        self, param, target1, target2, labelwidth=80, editwidth=80
    ):
        curr_hist_params = self.hist_params[param]

        val_min = curr_hist_params.val_min
        val_max = curr_hist_params.val_max

        # Label
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)

        # Add widget for unapplied thresholds
        curr_hist_params.unapply_widget_lower = eman2_gui.valslider.ValBox(
            self, (val_min, val_max), None, val_min, labelwidth
        )
        self.add_threshold_widget(
            curr_hist_params.unapply_widget_lower,
            temp_hbox,
            stylesheet="color: rgb(0,0,255);",
        )

        curr_hist_params.unapply_widget_upper = eman2_gui.valslider.ValBox(
            self, (val_min, val_max), None, val_max, labelwidth
        )
        self.add_threshold_widget(
            curr_hist_params.unapply_widget_upper,
            temp_hbox,
            stylesheet="color: rgb(255,0,0);",
        )
        target1.addLayout(temp_hbox, 0)

        # Label
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)

        # Add widget for applied thresholds
        curr_hist_params.apply_widget_lower = eman2_gui.valslider.ValBox(
            self, (val_min, val_max), None, val_min, labelwidth
        )
        self.add_threshold_widget(
            curr_hist_params.apply_widget_lower,
            temp_hbox,
            stylesheet="color: rgb(0,0,255);",
        )

        curr_hist_params.apply_widget_upper = eman2_gui.valslider.ValBox(
            self, (val_min, val_max), None, val_max, labelwidth
        )
        self.add_threshold_widget(
            curr_hist_params.apply_widget_upper,
            temp_hbox,
            stylesheet="color: rgb(255,0,0);",
        )
        target2.addLayout(temp_hbox, 0)

    def add_threshold_widget(
        self, widget, target, stylesheet="color: rgb(0,0,0);", editwidth=80
    ):
        # Add widget for unapplied thresholds
        widget.setEnabled(False)
        widget.text.setStyleSheet(stylesheet)
        widget.text.setMinimumSize(PyQt5.QtCore.QSize(editwidth, 0))
        target.addWidget(widget)

    def add_label_with_pulldown(self, labeltext, target, labelwidth=90, menuwidth=100):
        # Label
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)
        temp_label = PyQt5.QtWidgets.QLabel(labeltext, self)
        temp_label.setAlignment(
            PyQt5.QtCore.Qt.AlignRight | PyQt5.QtCore.Qt.AlignVCenter
        )
        temp_label.setFixedSize(PyQt5.QtCore.QSize(labelwidth, 20))
        temp_hbox.addWidget(temp_label)

        # Pulldown/ComboBox (may want to generalize this someday)
        self.sthresholdcontrol = PyQt5.QtWidgets.QComboBox(self)
        self.sthresholdcontrol.setMaximumWidth(menuwidth)
        for index, (control, map_entry) in enumerate(
            self.threshold_control_map.items()
        ):
            self.sthresholdcontrol.addItem(map_entry.label)
            self.sthresholdcontrol.setItemData(
                map_entry.index,
                PyQt5.QtGui.QColor(map_entry.color),
                PyQt5.QtCore.Qt.TextColorRole,
            )
        self.sthresholdcontrol.setCurrentIndex(self.curthresholdcontrol)
        temp_hbox.addWidget(self.sthresholdcontrol)

        target.addLayout(temp_hbox, 0)

    def add_two_buttons(self, button1, button2, target):
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)
        button1.setEnabled(False)
        temp_hbox.addWidget(button1)
        button2.setEnabled(False)
        temp_hbox.addWidget(button2)
        target.addLayout(temp_hbox, 0)

    def add_label_with_textbox(
        self, labeltext, stringbox, target, labelwidth=140
    ):  # , valueentry, target, labelwidth=90):
        # Label
        temp_hbox = PyQt5.QtWidgets.QHBoxLayout()
        temp_hbox.setContentsMargins(0, 0, 0, 0)
        temp_label = PyQt5.QtWidgets.QLabel(labeltext, self)
        temp_label.setAlignment(
            PyQt5.QtCore.Qt.AlignRight | PyQt5.QtCore.Qt.AlignVCenter
        )
        # temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
        temp_hbox.addWidget(temp_label)
        temp_hbox.addWidget(stringbox)
        target.addLayout(temp_hbox, 0)

    def readCterPartresFile(self, file_path):
        """Read all entries from a CTER partres file into the list box"""

        if not os.path.exists(file_path):
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                "Cannot find CTER partres file (%s). Please check the file path."
                % (file_path),
            )
            return

        if os.path.basename(file_path).find("partres") == -1:
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                'Invalid file name for CTER partres file (%s). The file name must contain "partres".'
                % (file_path),
            )
            return

        if file_path[-1 * len(".txt") :] != ".txt":
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                'Invalid file extension for CTER partres file (%s). The file extension must be ".txt".'
                % (file_path),
            )
            return

        new_entry_list = sp_utilities.read_text_row(file_path)
        if len(new_entry_list) == 0:
            PyQt5.QtWidgets.QMessageBox.warning(
                self,
                "Warning",
                "Specified CTER partres file (%s) does not contain any entry. Please check the file."
                % (file_path),
            )
            return
        assert len(new_entry_list) > 0, "MRK_DEBUG"

        if len(new_entry_list[0]) != self.n_idx_cter - self.n_idx_cter_extra:
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                "The number of columns (%d) has to be %d in %s."
                % (
                    len(new_entry_list[0]),
                    self.n_idx_cter - self.n_idx_cter_extra,
                    file_path,
                ),
            )
            return

        self.busy = True

        for cter_id in range(len(new_entry_list)):
            #
            # NOTE: 2017/11/22 Toshio Moriya
            # This CTER file format is the latest
            #
            if len(new_entry_list[cter_id]) == self.n_idx_cter - self.n_idx_cter_extra:
                # Add extra items first to make sure indices match
                extended_entry = []
                extended_entry = extended_entry + [cter_id]
                extended_entry = extended_entry + [1]

                extended_entry = (
                    extended_entry + new_entry_list[cter_id]
                )  # original entry

                ###				if self.is_enable_max_power == True:
                ###					extended_entry = extended_entry + [0.0] # MRK_TEST: self.idx_cter_max_power, <extra> maximum power in experimental rotational average (with astigmatism)

                # Store the extended entry to entry list
                new_entry_list[cter_id] = extended_entry

                ###				# MRK_TEST: Set max value of pwrot related to this micrograph
                ###				if self.is_enable_max_power == True: # MRK_TEST:
                ###					new_cter_mic_file_path = new_entry_list[cter_id][self.idx_cter_mic_name]
                ###					mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))
                ###					new_cter_pwrot_file_path = os.path.join(os.path.dirname(file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
                ###					new_rotinf_table = read_text_file(new_cter_pwrot_file_path, ncol=-1) # MRK_TEST:
                ###					new_entry_list[cter_id][self.idx_cter_max_power] = max(new_rotinf_table[self.idx_rotinf_exp_with_astig]) # MRK_TEST:

                # Always set selection state to 1 (selected)
                new_entry_list[cter_id][self.cter_params["select"].idx_cter] = 1

            else:
                assert False, (
                    "MRK_DEBUG: Found Invalid number of columns (%d) in %s"
                    % (len(new_entry_list[0]), file_path)
                )

        # now set the new status

        self.cter_partres_file_path = file_path
        self.cter_entry_list = new_entry_list

        # Loop through histogram parameters
        for idx_hist, (param, curr_hist_params) in enumerate(self.hist_params.items()):
            col_num = curr_hist_params.idx_cter
            val_min = min(self.cter_entry_list, key=lambda x: x[col_num])[col_num]
            val_min = round(val_min, self.round_ndigits)
            val_max = max(self.cter_entry_list, key=lambda x: x[col_num])[col_num]
            val_max = round(val_max, self.round_ndigits)
            curr_hist_params.val_min = val_min
            curr_hist_params.val_max = val_max
            curr_hist_params.unapply_threshold_lower = val_min
            curr_hist_params.unapply_threshold_upper = val_max
            curr_hist_params.unapply_widget_lower.setValue(val_min)
            curr_hist_params.unapply_widget_upper.setValue(val_max)
            curr_hist_params.apply_threshold_lower = val_min
            curr_hist_params.apply_threshold_upper = val_max
            curr_hist_params.apply_widget_lower.setValue(val_min)
            curr_hist_params.apply_widget_upper.setValue(val_max)

            idx_sort = self.cter_params_sort[param].idx_sort
            if val_min == val_max:
                self.shist.model().item(idx_hist).setEnabled(False)
                self.ssort.model().item(idx_sort).setEnabled(False)
                curr_hist_params.unapply_widget_lower.text.setStyleSheet(
                    "color: rgb(127,127,127);"
                )
                curr_hist_params.unapply_widget_upper.text.setStyleSheet(
                    "color: rgb(127,127,127);"
                )
                curr_hist_params.apply_widget_lower.text.setStyleSheet(
                    "color: rgb(127,127,127);"
                )
                curr_hist_params.apply_widget_upper.text.setStyleSheet(
                    "color: rgb(127,127,127);"
                )
            else:
                assert val_min < val_max, "MRK_DEBUG"
                self.shist.model().item(idx_hist).setEnabled(True)
                self.ssort.model().item(idx_sort).setEnabled(True)
                curr_hist_params.unapply_widget_lower.text.setStyleSheet(
                    "color: rgb(0,0,255);"
                )
                curr_hist_params.unapply_widget_upper.text.setStyleSheet(
                    "color: rgb(255,0,0);"
                )
                curr_hist_params.apply_widget_lower.text.setStyleSheet(
                    "color: rgb(0,0,255);"
                )
                curr_hist_params.apply_widget_upper.text.setStyleSheet(
                    "color: rgb(255,0,0);"
                )

        # Set disable status of histogram
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key, [1] is entry
        if selected_hist_param.val_min == selected_hist_param.val_max:
            idx_cter = selected_hist_param.idx_cter
            # self.curhistdisable=True
            self.curhistogramdisplay = False
            if self.whistparam.isVisible():
                self.whistparam.hide()
            if self.wscatterparam.isVisible():
                self.wscatterparam.hide()
            # Error message of this condition should be displayed at the end of this function for smooth visual presentation

            param_label = self.cter_params[selected_hist_param.param_name].label
            PyQt5.QtWidgets.QMessageBox.information(
                self,
                "Information",
                "All entries have the same selected parameter values (%s). \n\nParameter Histogram & Plot will not be shown."
                % (param_label),
            )

        self.updateEntryList()

        # Set the number of entries
        self.vbnentry.setValue(len(self.cter_entry_list))

        # Set the range of histogram bin
        self.vsentryperbin.setValue(self.curentryperbin)

        self.updateUncheckCounts()

        # Enable buttons
        self.pbrefreshgraphs.setEnabled(True)
        self.pbreapplysort.setEnabled(True)
        self.pbapplyallthreshold.setEnabled(True)
        self.pbsavethresholdset.setEnabled(True)
        self.pbloadthresholdset.setEnabled(True)
        self.pbsaveselection.setEnabled(True)

        cter_pwrot_dir = os.path.join(
            os.path.dirname(self.cter_partres_file_path), self.pwrot_dir
        )
        if os.path.exists(cter_pwrot_dir):
            # if not self.cbrotavgdisplay.getEnabled(): # MRK_NOTE: 2017/11/22 Toshio Moriya: This method does not work as I expected
            self.cbrotavgdisplay.setEnabled(True)
            for item in self.graph_map:
                self.graph_map[item].item_widget.setEnabled(True)

            self.vbplotfixscale.setEnabled(True)
        else:
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                'Cannot find "%s" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nPower spectrum rotational average plots display option is disabled for this session.'
                % (self.pwrot_dir, self.cter_partres_file_path),
            )

            # if self.cbrotavgdisplay.getEnabled(): # MRK_NOTE: 2017/11/22 Toshio Moriya: This method does not work as I expected
            self.cbrotavgdisplay.setEnabled(False)
            for item in self.graph_map:
                self.graph_map[item].item_widget.setEnabled(False)
            self.vbplotfixscale.setEnabled(False)
            # Error message of this condition should be displayed at the end of this function for smooth visual presentation
            # QMessageBox.warning(None,"Warning","Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nPower spectrum rotational average plots display option is disabled for this session." % (cter_pwrot_dir, self.cter_partres_file_path))

        cter_micthumb_dir = os.path.join(
            os.path.dirname(self.cter_partres_file_path), self.micthumb_dir
        )
        # print "MRK_DEBUG: cter_micthumb_dir = \"%s\" in readCterPartresFile() "% (cter_micthumb_dir)
        if os.path.exists(cter_micthumb_dir):
            # if not self.cbmicthumbdisplay.getEnabled(): # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
            self.cbmicthumbdisplay.setEnabled(True)
        else:
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                'Cannot find "%s" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nMicrograph thumbnail display option is disabled for this session.'
                % (cter_micthumb_dir, self.cter_partres_file_path),
            )

            # if self.cbmicthumbdisplay.getEnabled() != self.curimgmicthumbdisplay: # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
            self.cbmicthumbdisplay.setEnabled(False)
            # Error message of this condition should be displayed at the end of this function for smooth visual presentation
            # QMessageBox.warning(None,"Warning","Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nMicrograph thumbnail display option is disabled for this session." % (cter_micthumb_dir, self.cter_partres_file_path))

        # NOTE: 2016/01/03 Toshio Moriya
        # Force update related plots to hide too much scaling delay...
        self.updateImgMicThumb(False)
        self.updateHist()
        self.updatePlotParam()
        self.updateFFT()

        # Coarse plot
        if self.curplotrotavgdisplay:
            if not self.wplotrotavgcoarse.isVisible():
                self.wplotrotavgcoarse.show()
        else:
            assert not self.curplotrotavgdisplay, "MRK_DEBUG"
            if self.wplotrotavgcoarse.isVisible():
                self.wplotrotavgcoarse.hide()

        # Zoomed plot
        if self.curplotrotzoomdisplay:
            if not self.wplotrotavgfine.isVisible():
                self.wplotrotavgfine.show()
        else:
            assert not self.curplotrotavgdisplay, "MRK_DEBUG"
            if self.wplotrotavgfine.isVisible():
                self.wplotrotavgfine.hide()

        if self.curimgmicthumbdisplay:
            if not self.wimgmicthumb.isVisible():
                self.wimgmicthumb.show()
        else:
            assert not self.curimgmicthumbdisplay, "MRK_DEBUG"
            if self.wimgmicthumb.isVisible():
                self.wimgmicthumb.hide()

        if self.curfftdisplay:
            if not self.wfft.isVisible():
                self.wfft.show()
        else:
            assert not self.curfftdisplay, "MRK_DEBUG"
            if self.wfft.isVisible():
                self.wfft.hide()

        self.busy = False
        self.needredisp = True

    def openCterPartres(self, val=None):
        """Open CTER partres file"""
        file_path = PyQt5.QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Open CTER partres File",
            options=PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog,
        )

        # QFileDialog gives a tuple in Qt5, w/unicode elements
        if isinstance(file_path, tuple):
            file_path = file_path[0]
        else:
            file_path = str(file_path)
        if file_path == "":
            return

        self.readCterPartresFile(os.path.relpath(file_path))

    def updateHist(self, error_display=True):
        if self.whistparam == None:
            return  # it's closed/not visible
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected
        # if self.curhistdisable == True: return # do nothing while it is hidden
        if not self.curhistogramdisplay:
            return

        val_list = []

        # Create Histogram for selected parameter
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key, [1] is entry
        idx_cter = selected_hist_param.idx_cter
        for cter_entry in self.cter_entry_list:
            val_list.append(cter_entry[idx_cter])

        n_bin = 1
        if len(self.cter_entry_list) < self.curentryperbin:
            self.curentryperbin = len(self.cter_entry_list)
            self.vsentryperbin.setValue(self.curentryperbin)
            # self.vsentryperbin.setRange(1,len(self.cter_entry_list))
        elif self.curentryperbin < 1:
            self.curentryperbin = 1
            self.vsentryperbin.setValue(self.curentryperbin)
            # self.vsentryperbin.setRange(1,len(self.cter_entry_list))
        n_bin = old_div(
            len(self.cter_entry_list), self.curentryperbin
        )  # needs to be an integer
        assert len(val_list) >= n_bin, "MRK_DEBUG"
        assert n_bin > 0, "MRK_DEBUG"

        # Pad with zero for better visual impression...
        hist_x_list, hist_y_list = sp_statistics.hist_list(val_list, n_bin)
        hist_x_list += [max(val_list)]
        hist_y_list += [0]
        self.whistparam.set_data(
            (hist_x_list, hist_y_list),
            "hist_param",
            quiet=False,
            color=0,
            linetype=0,
            symtype=0,
        )

        # MRK_NOTE: 2015/12/17 Toshio Moriya
        # This may NOT be good place to update the following information...
        param_label = self.cter_params[selected_hist_param.param_name].label

        # Get value for current micrograph
        param_val = self.cter_entry_list[self.curentry][idx_cter]

        # Get ordinate range in order to overlay current micrograph's value on histogram
        val_min = round(min(hist_y_list), self.round_ndigits)
        val_max = round(max(hist_y_list), self.round_ndigits)

        unapply_threshold_lower_val = round(
            selected_hist_param.unapply_threshold_lower, self.round_ndigits
        )
        apply_threshold_lower_val = round(
            selected_hist_param.apply_threshold_lower, self.round_ndigits
        )
        unapply_threshold_upper_val = round(
            selected_hist_param.unapply_threshold_upper, self.round_ndigits
        )
        apply_threshold_upper_val = round(
            selected_hist_param.apply_threshold_upper, self.round_ndigits
        )

        self.whistparam.set_data(
            ([param_val, param_val], [val_min, val_max]),
            "selected_val",
            quiet=False,
            color=3,
            linetype=1,
        )  # default linetype was hard to see
        self.whistparam.set_data(
            (
                [unapply_threshold_lower_val, unapply_threshold_lower_val],
                [val_min, val_max],
            ),
            "unapply_threshold_lower_val",
            quiet=False,
            color=1,
            linetype=1,
        )
        self.whistparam.set_data(
            (
                [apply_threshold_lower_val, apply_threshold_lower_val],
                [val_min, val_max],
            ),
            "apply_threshold_lower_val",
            quiet=False,
            color=1,
        )
        self.whistparam.set_data(
            (
                [unapply_threshold_upper_val, unapply_threshold_upper_val],
                [val_min, val_max],
            ),
            "unapply_threshold_upper_val",
            quiet=False,
            color=2,
            linetype=1,
        )
        self.whistparam.set_data(
            (
                [apply_threshold_upper_val, apply_threshold_upper_val],
                [val_min, val_max],
            ),
            "apply_threshold_upper_val",
            quiet=False,
            color=2,
        )

        self.whistparam.setAxisParms(param_label, "Image Counts")
        # x_margin = (hist_x_list[-1] - hist_x_list[0]) * 0.05
        # NOTE: 2016/01/02 Toshio Moriya
        # Disable manual rescale for now and use autoscale
        # self.whistparam.rescale(min(val_list),max(val_list),0,max(hist_y_list) * 1.05)
        self.whistparam.autoscale(True)

    def updatePlotParam(self, error_display=True):
        if self.wscatterparam == None:
            return  # it's closed/not visible
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected
        if not self.curscatterdisplay and not self.curhistogramdisplay:
            return

        x_list = []
        y_list = []

        # Create graph for selected parameter
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key, [1] is entry
        idx_cter = selected_hist_param.idx_cter
        for cter_id in range(len(self.cter_entry_list)):
            x_list.append(cter_id)
            y_list.append(self.cter_entry_list[cter_id][idx_cter])
        # self.wscatterparam.set_data((x_list,y_list),"plot_param",quiet=False,color=0)
        self.wscatterparam.set_data(
            (x_list, y_list), "plot_param", quiet=False, color=0, linetype=0, symtype=0
        )

        # Create graph for single parameter value of selected entry
        # MRK_NOTE: 2015/12/17 Toshio Moriya
        # This may NOT be good place to update the following information...
        param_label = self.cter_params[selected_hist_param.param_name].label
        param_val = round(
            self.cter_entry_list[self.curentry][idx_cter], self.round_ndigits
        )

        unapply_threshold_lower_val = round(
            selected_hist_param.unapply_threshold_lower, self.round_ndigits
        )
        apply_threshold_lower_val = round(
            selected_hist_param.apply_threshold_lower, self.round_ndigits
        )
        unapply_threshold_upper_val = round(
            selected_hist_param.unapply_threshold_upper, self.round_ndigits
        )
        apply_threshold_upper_val = round(
            selected_hist_param.apply_threshold_upper, self.round_ndigits
        )

        y_list = [param_val] * len(x_list)
        self.wscatterparam.set_data(
            (x_list, y_list), "selected_val", quiet=False, color=3, linetype=1
        )  # default linetype was hard to see
        y_list = [unapply_threshold_lower_val] * len(x_list)
        self.wscatterparam.set_data(
            (x_list, y_list),
            "unapply_threshold_lower_val",
            quiet=False,
            color=1,
            linetype=1,
        )
        y_list = [apply_threshold_lower_val] * len(x_list)
        self.wscatterparam.set_data(
            (x_list, y_list), "apply_threshold_lower_val", quiet=False, color=1
        )
        y_list = [unapply_threshold_upper_val] * len(x_list)
        self.wscatterparam.set_data(
            (x_list, y_list),
            "unapply_threshold_upper_val",
            quiet=False,
            color=2,
            linetype=1,
        )
        y_list = [apply_threshold_upper_val] * len(x_list)
        self.wscatterparam.set_data(
            (x_list, y_list), "apply_threshold_upper_val", quiet=False, color=2
        )

        self.wscatterparam.setAxisParms("Sorted Image ID", param_label)
        # NOTE: 2016/01/02 Toshio Moriya
        # Use autoscale for now
        self.wscatterparam.autoscale(True)

    def updatePlotCurves(self, error_display=True):
        if not self.curplotrotavgdisplay and not self.curplotrotzoomdisplay:
            return  # Both power-spectrum rotational-average plot displays are unchecked
        if self.wplotrotavgcoarse == None and self.wplotrotavgfine == None:
            return  # closed/not visible
        if self.cter_pwrot_file_path == None:
            return  # no cter entry is selected

        # print "MRK_DEBUG: self.cter_pwrot_file_path =\"%s\" in updatePlotCurves() "% (self.cter_pwrot_file_path)
        if not os.path.exists(self.cter_pwrot_file_path):
            if self.wplotrotavgcoarse.isVisible():
                self.wplotrotavgcoarse.hide()
            if self.wplotrotavgfine.isVisible():
                self.wplotrotavgfine.hide()
            if self.wfft.isVisible():
                self.wfft.hide()
            if error_display and os.path.exists(
                os.path.dirname(self.cter_pwrot_file_path)
            ):
                PyQt5.QtWidgets.QMessageBox.warning(
                    None,
                    "Warning",
                    "Cannot find file cter_pwrot_file_path (%s). Please check the contents of pwrot directory. \n\nPlots will not be shown."
                    % (self.cter_pwrot_file_path),
                )
            return
        assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"

        if not self.use_ctffind:
            # Now update the plots
            self.rotinf_table = sp_utilities.read_text_file(
                self.cter_pwrot_file_path, ncol=-1
            )
            ncols = len(
                self.rotinf_table
            )  # old rotinf file has 6 columns, 8 with baseline-fitting

            # Subtract background and apply envelope
            if ncols == 6:
                if not self.checkedpwrot:
                    self.checkedpwrot = True
                newCurveList = self.fitSpline()
                self.rotinf_table = self.rotinf_table + newCurveList

            # Spatial frequency
            spFreqList = self.rotinf_table[self.rot1d_indices["freq"]]

            # Loop through power-spectrum profiles
            for index, (item, graph) in enumerate(self.graph_map.items()):
                columnValues = self.rotinf_table[graph.idx_rotinf]
                self.wplotrotavgcoarse.set_data(
                    (spFreqList, columnValues),
                    graph.item_name,
                    quiet=False,
                    color=index,
                    linetype=0,
                )
                self.wplotrotavgfine.set_data(
                    (spFreqList, columnValues),
                    graph.item_name,
                    quiet=False,
                    color=index,
                    linetype=0,
                )
        else:
            # Now update the plots
            file_object = open(self.cter_pwrot_file_path, "r")
            self.rotinf_table = file_object.readlines()[5:]  # skip lines 1-5
            ncols = 6

            # Spatial frequency
            spFreqList = [
                float(i) for i in self.rotinf_table[0].split()
            ]  # str by default

            # Loop through power-spectrum profiles
            for index, (item, graph) in enumerate(self.graph_map.items()):
                columnValues = [
                    float(i) for i in self.rotinf_table[graph.idx_rotinf].split()
                ]  # str by default
                self.wplotrotavgcoarse.set_data(
                    (spFreqList, columnValues),
                    graph.item_name,
                    quiet=False,
                    color=index,
                    linetype=0,
                )
                self.wplotrotavgfine.set_data(
                    (spFreqList, columnValues),
                    graph.item_name,
                    quiet=False,
                    color=index,
                    linetype=0,
                )

        # NOTE: 2016/01/02 Toshio Moriya
        # Disable manual rescale for now and use autoscale
        # self.wplotrotavgcoarse.rescale(spFreqList[0],spFreqList[-1],0.0,1.0)
        # self.wplotrotavgcoarse.autoscale(True)
        self.wplotrotavgfine.rescale(
            spFreqList[0], spFreqList[-1], 0.0, self.curplotfixscale
        )

        nyquist_freq = spFreqList[-1]

        # Draw limits on plots
        self.draw_limits("error_astig", "Astig. Limit", 0, 0, 0.5, 18, nyquist_freq)
        self.draw_limits("error_def", "Defocus Limit", 0.5, 0, 0, 36, nyquist_freq)
        self.draw_limits("error_ctf", "CTF Limit", 0, 0.5, 0, 54, nyquist_freq)

        self.wplotrotavgcoarse.setAxisParms(
            "frequency (1/" + "$\AA$" + ")", "power spectrum"
        )
        self.wplotrotavgfine.setAxisParms(
            "frequency (1/" + "$\AA$" + ")", "power spectrum"
        )

        # Update plot
        self.updatePlotVisibility()

    def fitSpline(self, bkgdSmooth=0.15, envSmooth=0):
        """Subtract background from experimental curve and apply envelope to theoretical curve
		by fitting relative extrema to a spline"""

        # fit_with_astig
        fitWAstig = self.rotinf_table[self.graph_map["fit_with_astig"].idx_rotinf]

        # exp_with_astig
        expWAstig = self.rotinf_table[self.graph_map["exp_with_astig"].idx_rotinf]

        # Spatial frequency
        spFreqList = self.rotinf_table[self.rot1d_indices["freq"]]

        # initialize extrema
        expMaxX = []
        expMaxY = []
        fitMinX = []
        fitMinY = []
        expMinX = []
        expMinY = []

        # Search experimental and theoretical curves for extrema
        for key in range(2, len(spFreqList) - 2):
            ## Check if fitted curve is at a maximum
            # if not fitWAstig[key] > fitWAstig[key-2] and fitWAstig[key] > fitWAstig[key+2]:
            ## remember value
            # firstMaxKey = key  # not used yet
            ##print('First max.', key, spFreqList[key], fitWAstig[key], expWAstig[key])

            # Check if experimental curve is at a minimum  (skipping low-resolution values)
            if (
                expWAstig[key] < expWAstig[key - 2]
                and expWAstig[key] < expWAstig[key + 2]
            ):
                expMinX.append(spFreqList[key])
                expMinY.append(expWAstig[key])
                # print('Minimum', key, spFreqList[key], fitWAstig[key], expWAstig[key])

            # Check if experimental curve is at a maximum
            if (
                expWAstig[key] > expWAstig[key - 2]
                and expWAstig[key] > expWAstig[key + 2]
            ):
                expMaxX.append(spFreqList[key])
                expMaxY.append(expWAstig[key])
                # print('Maximum',spFreqList[key],expWAstig[key],expWAstig[key-1],expWAstig[key+1])

            # Check if fitted curve is at a minimum
            if (
                fitWAstig[key] < fitWAstig[key - 1]
                and fitWAstig[key] < fitWAstig[key + 1]
            ):
                fitMinX.append(spFreqList[key])
                fitMinY.append(fitWAstig[key])
                # print('Rel. minimum', key, spFreqList[key], fitWAstig[key], expWAstig[key])

        # Change to arrays
        expMaxXNP = numpy.array(expMaxX)
        expMaxYNP = numpy.array(expMaxY)
        fitMinXNP = numpy.array(fitMinX)
        fitMinYNP = numpy.array(fitMinY)
        expMinXNP = numpy.array(expMinX)
        expMinYNP = numpy.array(expMinY)

        xNP = numpy.array(spFreqList)
        # print("xNP",xNP)

        # Sanity check
        if len(fitMinXNP) == 0:
            print("WARNING!! 1D profile has no relative minima. There's probably something weird going on.")
            fitMinXNP = numpy.array(spFreqList)
            fitMinYNP = numpy.array(fitWAstig)
        
        # Splinefit extrema
        fitMinTck = scipy.interpolate.splrep(
            fitMinXNP, fitMinYNP, s=bkgdSmooth
        )  # s==smoothing factor
        ####expMaxTck = scipy.interpolate.splrep(expMaxXNP, expMaxYNP, s=envSmooth)

        # Evaluate splines
        fitMinSpline = scipy.interpolate.splev(xNP, fitMinTck)
        # print("fitMinSpline", type(fitMinSpline), len(fitMinSpline), fitMinSpline)
        ####expMaxSpline = scipy.interpolate.splev(xNP, expMaxTck)  # not used
        expMinSpline = []
        expMin = min(expWAstig[3:])  # value at origin will be zero

        ## Exponential fit (NOT USED)
        ##envelopeExp = optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  expMaxXNP,  expMaxYNP)
        ##print('envelopeExp', float(envelopeExp[0][0]), type(float(envelopeExp[0][0])), float(envelopeExp[0][1]), type(float(envelopeExp[0][1])))
        # envelopeExp = optimize.curve_fit(lambda t,a,b,c: a*np.exp(b*t)+c,  expMaxXNP,  expMaxYNP)[0]  # 2nd element is numpy.array
        ##print('envelopeExp', envelopeExp[0], type(envelopeExp[0]), envelopeExp[1], type(envelopeExp[1]), envelopeExp[2], type(envelopeExp[2]))

        # will subtract background and apply envelope
        expSubtract = []
        theorEnv = []

        # Subtract minimum and multiply by envelope
        for key in range(0, len(spFreqList)):
            expMinSpline.append(expMin)
            diff = expWAstig[key] - expMinSpline[key]
            # print(spFreq[key],float(diff))
            expSubtract.append(diff)

            interp = (fitWAstig[key] - fitMinSpline[key]) * diff

            ##exponential = envelopeExp[0]*np.exp(envelopeExp[1]*spFreqList[key]) + envelopeExp[2]
            # exponential = fitWAstig[key]*envelopeExp[0]*np.exp(envelopeExp[1]*spFreqList[key])
            ##print('exponential', key, spFreqList[key], diff, exponential, type(exponential))

            theorEnv.append(interp)

        ## Make sure splinefit-corrected experimental curve stays above zero (NOT USED, didn't like)
        # minExpSpline = min(expMinSpline[1:])  # origin might have a funny value
        ##print("minExpSpline",minExpSpline)
        ##print(type(expSubtract))
        # maxExpSpline = max(expSubtract[3:])  # origin might have a funny value
        ##print('maxExpSpline',maxExpSpline)
        ##print('firstMaxKey', firstMaxKey, expSubtract[firstMaxKey])
        # maxTheorSubtract = max(theorEnv)
        ##print('maxTheorSubtract',maxTheorSubtract)

        combinedList = []
        combinedList.append(expSubtract)
        combinedList.append(theorEnv)

        return combinedList

    def draw_limits(
        self, error_name, error_label, shape_r, shape_g, shape_b, y_offset, nyquist_freq
    ):
        # (I don't know why I need the following, but on some machines, scrlim is undefined. --Tapu, 2019/04/03)
        if not hasattr(self.wplotrotavgfine, "scrlim"):
            self.wplotrotavgfine.render()
        if not hasattr(self.wplotrotavgcoarse, "scrlim"):
            self.wplotrotavgcoarse.render()
        if not hasattr(self.whistparam, "scrlim"):
            self.whistparam.render()
        if not hasattr(self.wscatterparam, "scrlim"):
            self.wscatterparam.render()

        fineLimits = self.wplotrotavgfine.scrlim
        coarseLimits = self.wplotrotavgcoarse.scrlim

        error_freq = self.cter_entry_list[self.curentry][
            self.cter_params[error_name].idx_cter
        ]
        if error_freq > 0.0 and error_freq <= nyquist_freq:
            error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
            self.wplotrotavgcoarse.add_shape(
                error_name,
                eman2_gui.emshape.EMShape(
                    (
                        "scrline",
                        shape_r,
                        shape_g,
                        shape_b,
                        error_scr_x,
                        coarseLimits[1],
                        error_scr_x,
                        coarseLimits[1] + coarseLimits[3],
                        1,
                    )
                ),
            )
            self.wplotrotavgcoarse.add_shape(
                "%s_freq" % (error_name),
                eman2_gui.emshape.EMShape(
                    (
                        "scrlabel",
                        0,
                        0,
                        0,
                        error_scr_x - 260,
                        coarseLimits[1] + coarseLimits[3] - y_offset,
                        "%s %1.5g (%1.5g)"
                        % (error_label, error_freq, old_div(1.0, error_freq)),
                        120.0,
                        -1,
                    )
                ),
            )
            error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
            self.wplotrotavgfine.add_shape(
                error_name,
                eman2_gui.emshape.EMShape(
                    (
                        "scrline",
                        shape_r,
                        shape_g,
                        shape_b,
                        error_scr_x,
                        fineLimits[1],
                        error_scr_x,
                        fineLimits[1] + fineLimits[3],
                        1,
                    )
                ),
            )
            self.wplotrotavgfine.add_shape(
                "%s_freq" % (error_name),
                eman2_gui.emshape.EMShape(
                    (
                        "scrlabel",
                        0,
                        0,
                        0,
                        error_scr_x - 260,
                        fineLimits[1] + fineLimits[3] - y_offset,
                        "%s %1.5g (%1.5g)"
                        % (error_label, error_freq, old_div(1.0, error_freq)),
                        120.0,
                        -1,
                    )
                ),
            )

    def newRotAvgDisplay(self, val=None):
        """Change rotational average plot display status."""
        # assert self.cbrotavgdisplay.getEnabled() == True, "MRK_DEBUG: 2017/11/22 Toshio Moriya: This method does not work as I expected"
        if self.curplotrotavgdisplay != val:
            # now set the new display status
            self.curplotrotavgdisplay = val
        else:
            assert self.curplotrotavgdisplay == val, "MRK_DEBUG"
            # The status did not change, there is nothing to do nothing
            return

        if self.cter_pwrot_file_path == None:
            return  # no cter entry is selected

        if not os.path.exists(self.cter_pwrot_file_path):
            assert not self.wplotrotavgcoarse.isVisible(), "MRK_DEBUG"
            return
        assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"

        if self.curplotrotavgdisplay and not self.wplotrotavgcoarse.isVisible():
            self.wplotrotavgcoarse.show()
            self.needredisp = True
        elif not self.curplotrotavgdisplay and self.wplotrotavgcoarse.isVisible():
            self.wplotrotavgcoarse.hide()

    def newRotZoomDisplay(self, val=None):
        """Change rotational average plot display status."""

        if self.curplotrotzoomdisplay != val:
            # now set the new display status
            self.curplotrotzoomdisplay = val
        else:
            assert self.curplotrotzoomdisplay == val, "MRK_DEBUG"
            # The status did not change, there is nothing to do
            return

        if self.cter_pwrot_file_path == None:
            return  # no cter entry is selected

        if not os.path.exists(self.cter_pwrot_file_path):
            assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
            return
        assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"

        if self.curplotrotzoomdisplay and not self.wplotrotavgfine.isVisible():
            self.wplotrotavgfine.show()
            self.needredisp = True
        elif not self.curplotrotzoomdisplay and self.wplotrotavgfine.isVisible():
            self.wplotrotavgfine.hide()

    def newHistogramDisplay(self, val=None):
        """Change histogram display status."""

        if self.curhistogramdisplay != val:
            # now set the new display status
            self.curhistogramdisplay = val
        else:
            assert self.curhistogramdisplay == val, "MRK_DEBUG"
            # The status did not change, there is nothing to do
            return

        if self.cter_pwrot_file_path == None:
            return  # no cter entry is selected

        if not os.path.exists(self.cter_pwrot_file_path):
            assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
            return
        assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"

        if self.curhistogramdisplay and not self.whistparam.isVisible():
            self.whistparam.show()
            self.needredisp = True
        elif not self.curhistogramdisplay and self.whistparam.isVisible():
            self.whistparam.hide()

    def newScatterDisplay(self, val=None):
        """Change sort plot display status."""

        if self.curscatterdisplay != val:
            # now set the new display status
            self.curscatterdisplay = val
        else:
            assert self.curscatterdisplay == val, "MRK_DEBUG"
            # The status did not change, there is nothing to do
            return

        if self.cter_pwrot_file_path == None:
            return  # no cter entry is selected

        if not os.path.exists(self.cter_pwrot_file_path):
            assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
            return
        assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"

        if self.curscatterdisplay and not self.wscatterparam.isVisible():
            self.wscatterparam.show()
            self.needredisp = True
        elif not self.curscatterdisplay and self.wscatterparam.isVisible():
            self.wscatterparam.hide()

    def updateEntryList(self):
        """Updated entry list box after sorting of CTER partres entries based on current setting."""

        # sort CTER partres entry list
        assert self.cter_entry_list != None, "MRK_DEBUG"
        param_name = list(self.cter_params_sort.items())[self.cursortidx][
            1
        ].param_name  # cter_params_sort is an OrderedDict, [0] is key, [1] is entry
        if self.cursortidx != self.cter_params_sort["mic_name"].idx_sort:
            self.cter_entry_list = sorted(
                self.cter_entry_list,
                key=lambda x: x[
                    list(self.cter_params_sort.items())[self.cursortidx][1].idx_cter
                ],
                reverse=self.cursortorder,
            )
        else:
            assert self.cursortidx == self.cter_params_sort["mic_name"].idx_sort
            self.cter_entry_list = sorted(
                self.cter_entry_list,
                key=lambda x: os.path.basename(
                    x[list(self.cter_params_sort.items())[self.cursortidx][1].idx_cter]
                ),
                reverse=self.cursortorder,
            )

        if self.cursortselect:
            # Additionaly, sort cter entry list by select state
            self.cter_entry_list = sorted(
                self.cter_entry_list,
                key=lambda x: x[self.cter_params["select"].idx_cter],
            )
        # else: # Do nothing

        # Refresh entry list box
        self.lbentry.clear()
        newItemflags = (
            PyQt5.QtCore.Qt.ItemFlags(PyQt5.QtCore.Qt.ItemIsSelectable)
            | PyQt5.QtCore.Qt.ItemFlags(PyQt5.QtCore.Qt.ItemIsEnabled)
            | PyQt5.QtCore.Qt.ItemFlags(PyQt5.QtCore.Qt.ItemIsUserCheckable)
        )
        for cter_entry in self.cter_entry_list:
            newItem = PyQt5.QtWidgets.QListWidgetItem(
                os.path.basename(cter_entry[self.cter_params["mic_name"].idx_cter])
            )
            newItem.setFlags(newItemflags)
            if cter_entry[self.cter_params["select"].idx_cter] == 1:
                newItem.setCheckState(PyQt5.QtCore.Qt.Checked)
            else:
                assert cter_entry[self.cter_params["select"].idx_cter] == 0, "MRK_DEBUG"
                newItem.setCheckState(PyQt5.QtCore.Qt.Unchecked)
            self.lbentry.addItem(newItem)
            # self.lbentry.addItem(os.path.basename(cter_entry[self.cter_params['mic_name'].idx_cter]))

        self.newEntry(0)
        self.lbentry.setCurrentRow(0)

    def updateImgMicThumb(self, error_display=True):
        if not self.curimgmicthumbdisplay:
            return  # Micrograph thumbnail display is unchecked
        if self.wimgmicthumb == None:
            return  # it's closed/not visible
        if self.cter_micthumb_file_path == None:
            return  # no cter entry is selected

        if not os.path.exists(self.cter_micthumb_file_path):
            if self.wimgmicthumb.isVisible():
                self.wimgmicthumb.hide()
            if error_display and os.path.exists(
                os.path.dirname(self.cter_micthumb_file_path)
            ):
                PyQt5.QtWidgets.QMessageBox.warning(
                    None,
                    "Warning",
                    "Cannot find micrograph thumbnail (%s). Please check your micrograph thumbnail directory. \n\nMicrograph thumbnail will not be shown."
                    % (self.cter_micthumb_file_path),
                )
            return
        assert os.path.exists(self.cter_micthumb_file_path), "MRK_DEBUG"

        # Now update the image
        micthumb_img = EMAN2_cppwrap.EMData(
            self.cter_micthumb_file_path
        )  # read the image from disk
        self.wimgmicthumb.set_data(micthumb_img)
        self.wimgmicthumb.setWindowTitle(
            "sp_gui_cter - Micrograph Thumbnail- %s, %s"
            % (
                os.path.basename(
                    self.cter_entry_list[self.curentry][
                        self.cter_params["mic_name"].idx_cter
                    ]
                ),
                os.path.basename(self.cter_micthumb_file_path),
            )
        )

    def newMicThumbDisplay(self, val=None):
        """Change micrograph thumbnail display status."""
        # assert self.cbmicthumbdisplay.getEnabled() == True, "MRK_DEBUG: 2016/03/22 Toshio Moriya: This method does not work as I expected"
        if self.curimgmicthumbdisplay != val:
            # now set the new display status
            self.curimgmicthumbdisplay = val
        else:
            assert self.curimgmicthumbdisplay == val, "MRK_DEBUG"
            # The status did not change, there is nothing to do nothing
            return

        if self.cter_micthumb_file_path == None:
            return  # no cter entry is selected

        # print "MRK_DEBUG: self.cter_micthumb_file_path =\"%s\" in newMicThumbDisplay() "% (self.cter_micthumb_file_path)
        if not os.path.exists(self.cter_micthumb_file_path):
            assert not self.wimgmicthumb.isVisible(), "MRK_DEBUG"
            return
        assert os.path.exists(self.cter_micthumb_file_path), "MRK_DEBUG"

        if self.curimgmicthumbdisplay and not self.wimgmicthumb.isVisible():
            self.wimgmicthumb.show()
            self.needredisp = True
        elif not self.curimgmicthumbdisplay and self.wimgmicthumb.isVisible():
            self.wimgmicthumb.hide()

    def updateFFT(self, error_display=True):
        if not self.curfftdisplay:
            return  # FFT display is unchecked
        if self.wfft == None:
            return  # it's closed/not visible
        if self.cter_fft_file_path == None:
            # Try directory of partres file
            pwsdir = os.path.join(
                os.path.dirname(self.cter_partres_file_path), self.power2d_dir
            )
            return

        if not os.path.exists(self.cter_fft_file_path):
            if self.wfft.isVisible():
                self.wfft.hide()
            if error_display and os.path.exists(
                os.path.dirname(self.cter_fft_file_path)
            ):
                PyQt5.QtWidgets.QMessageBox.warning(
                    None,
                    "Warning",
                    "Cannot find power spectrum (%s). Please check your power spectrum directory. \n\nPower spectrum will not be shown."
                    % (self.cter_micthumb_file_path),
                )
            return
        assert os.path.exists(self.cter_fft_file_path), "MRK_DEBUG"

        # Now update the image
        fft_img = EMAN2_cppwrap.EMData(
            self.cter_fft_file_path
        )  # read the image from disk
        self.wfft.set_data(fft_img)
        self.wfft.setWindowTitle(
            "sp_gui_cter - 2D FFT %s, %s"
            % (
                os.path.basename(
                    self.cter_entry_list[self.curentry][
                        self.cter_params["mic_name"].idx_cter
                    ]
                ),
                os.path.basename(self.cter_fft_file_path),
            )
        )

    def newFFTDisplay(self, val=None):
        """Change FFT display status."""
        if self.curfftdisplay != val:
            # now set the new display status
            self.curfftdisplay = val
        else:
            assert self.curfftdisplay == val, "MRK_DEBUG"
            # The status did not change, there is nothing to do
            return

        if self.cter_fft_file_path == None:
            return  # no cter entry is selected

        if not os.path.exists(self.cter_fft_file_path):
            assert not self.wfft.isVisible(), "MRK_DEBUG"
            pwsdir = os.path.join(
                os.path.dirname(self.cter_partres_file_path), self.power2d_dir
            )
            PyQt5.QtWidgets.QMessageBox.warning(
                None,
                "Warning",
                'Cannot find "%s" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nPower-spectrum display option is disabled for this session.'
                % (pwsdir, self.cter_partres_file_path),
            )
            self.curfftdisplay = False
            self.cbfftdisplay.setValue(False)
            self.cbfftdisplay.setEnabled(False)
            return
        assert os.path.exists(self.cter_fft_file_path), "MRK_DEBUG"

        if self.curfftdisplay and not self.wfft.isVisible():
            self.wfft.show()
            self.needredisp = True
        elif not self.curfftdisplay and self.wfft.isVisible():
            self.wfft.hide()

    def newEntry(self, currow):
        """called when a new data set is selected from the CTER partres entry list box."""
        assert self.cter_partres_file_path != None, "MRK_DEBUG"
        assert self.cter_entry_list != None, "MRK_DEBUG"

        # always update the current row of cter entry list
        # to get associated micrograph path and pwrot file path
        self.curentry = (
            currow
        )  # row can be the same even after resorting of the cter entry list

        # Get associated micrograph path of current entry
        new_cter_mic_file_path = self.cter_entry_list[self.curentry][
            self.cter_params["mic_name"].idx_cter
        ]

        # Generate associated micthumb & pwrot file path of current entry
        mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))[
            0
        ]
        new_cter_micthumb_file_path = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            self.micthumb_dir,
            mic_basename_root + self.micthumb_suffix,
        )
        new_cter_pwrot_file_path = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            self.pwrot_dir,
            mic_basename_root + self.pwrot_suffix,
        )
        new_cter_fft_file_path = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            self.power2d_dir,
            mic_basename_root + self.power2d_suffix,
        )

        # Changing row does not always change the pwrot file path after resorting of the cter entry list
        # If same, skip the following processes
        if self.cter_pwrot_file_path == new_cter_pwrot_file_path:
            assert (
                self.cter_micthumb_file_path == new_cter_micthumb_file_path
            ), "MRK_DEBUG"
            assert self.cter_mic_file_path == new_cter_mic_file_path, "MRK_DEBUG"
            assert self.cter_fft_file_path == new_cter_fft_file_path, "MRK_DEBUG"
            return

        # now set the new item
        assert self.cter_pwrot_file_path != new_cter_pwrot_file_path, "MRK_DEBUG"
        self.cter_pwrot_file_path = new_cter_pwrot_file_path
        assert self.cter_micthumb_file_path != new_cter_micthumb_file_path, "MRK_DEBUG"
        self.cter_micthumb_file_path = new_cter_micthumb_file_path
        assert self.cter_mic_file_path != new_cter_mic_file_path, "MRK_DEBUG"
        self.cter_mic_file_path = new_cter_mic_file_path
        assert self.cter_fft_file_path != new_cter_fft_file_path, "MRK_DEBUG"
        self.cter_fft_file_path = new_cter_fft_file_path

        self.ssortedid.setValue(self.curentry, True)

        for param in self.cter_params.keys():
            current_param = self.cter_params[param]
            if current_param.widget is not None:
                current_param.widget.setValue(
                    self.cter_entry_list[self.curentry][current_param.idx_cter], True
                )

        # Use blue (lower) & red (higher) fonts to indicate the value is not between applied threshold ranage
        for param in self.hist_params.keys():
            lower_threshold = round(
                self.hist_params[param].apply_threshold_lower, self.round_ndigits
            )
            upper_threshold = round(
                self.hist_params[param].apply_threshold_upper, self.round_ndigits
            )
            idx_cter = self.hist_params[param].idx_cter
            param_val = round(
                self.cter_entry_list[self.curentry][idx_cter], self.round_ndigits
            )
            if lower_threshold < upper_threshold:
                if lower_threshold <= param_val and param_val <= upper_threshold:
                    self.cter_params[
                        self.hist_params[param].param_name
                    ].widget.text.setStyleSheet("color: rgb(0,0,0);")
                elif param_val < lower_threshold:
                    self.cter_params[
                        self.hist_params[param].param_name
                    ].widget.text.setStyleSheet("color: rgb(0,0,255);")
                else:
                    assert upper_threshold < param_val, "MRK_DEBUG"
                    self.cter_params[
                        self.hist_params[param].param_name
                    ].widget.text.setStyleSheet("color: rgb(255,0,0);")
            else:
                assert lower_threshold == upper_threshold, "MRK_DEBUG"
                assert (
                    lower_threshold == param_val and param_val == upper_threshold
                ), "MRK_DEBUG"
                self.cter_params[
                    self.hist_params[param].param_name
                ].widget.text.setStyleSheet("color: rgb(127,127,127);")

        self.wplotrotavgcoarse.setWindowTitle(
            "sp_gui_cter - Plot - %s, %s"
            % (
                os.path.basename(
                    self.cter_entry_list[self.curentry][
                        self.cter_params["mic_name"].idx_cter
                    ]
                ),
                os.path.basename(self.cter_pwrot_file_path),
            )
        )
        self.wplotrotavgfine.setWindowTitle(
            "sp_gui_cter - Plot Zoom- %s, %s"
            % (
                os.path.basename(
                    self.cter_entry_list[self.curentry][
                        self.cter_params["mic_name"].idx_cter
                    ]
                ),
                os.path.basename(self.cter_pwrot_file_path),
            )
        )

        self.needredisp = True

    def updateEntrySelect(self, entry):
        """called when check status of an cter entry in list box is changed."""
        assert self.cter_partres_file_path != None, "MRK_DEBUG"
        assert self.cter_entry_list != None, "MRK_DEBUG"

        newSelect = 1
        if entry.checkState() == PyQt5.QtCore.Qt.Unchecked:
            newSelect = 0
        entry_row = self.lbentry.row(entry)
        self.cter_entry_list[entry_row][self.cter_params["select"].idx_cter] = newSelect

        if self.curentry == entry_row:
            self.cter_params["select"].widget.setValue(
                self.cter_entry_list[self.curentry][
                    self.cter_params["select"].idx_cter
                ],
                True,
            )

        self.updateUncheckCounts()

    def updateUncheckCounts(self):
        """called whenever checked status of cter entries change."""
        assert self.cter_partres_file_path != None, "MRK_DEBUG"
        assert self.cter_entry_list != None, "MRK_DEBUG"

        assert len(self.cter_entry_list) > 0, "MRK_DEBUG"
        n_entry = len(self.cter_entry_list)
        uncheck_counts = n_entry
        for cter_entry in self.cter_entry_list:
            uncheck_counts -= cter_entry[self.cter_params["select"].idx_cter]
        assert uncheck_counts >= 0 and uncheck_counts <= n_entry, "MRK_DEBUG"

        self.vbuncheckcounts.setValue(uncheck_counts, True)
        self.vbuncheckratio.setValue(old_div(float(uncheck_counts), n_entry), True)

    def reapplySort(self, item=None):
        """Called when reapply button is clicked."""
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        self.updateEntryList()

    def newSort(self, cursortidx):
        """Sort CTER partres entries by selected parameter values."""
        if self.cursortidx == cursortidx:
            return

        # now set the new item
        self.cursortidx = cursortidx

        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        self.updateEntryList()

    def newSortOrder(self, sortorder):
        """Change sorting order of CTER partres entries."""
        if self.cursortorder == sortorder:
            return

        # now set the new status
        self.cursortorder = sortorder

        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        self.updateEntryList()

    def newSortSelect(self, sortselect):
        """Change sort select status of CTER partres entries."""
        if self.cursortselect == sortselect:
            return

        # now set the new status
        self.cursortselect = sortselect

        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        self.updateEntryList()

    def newThresholdLower(self):
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key, [1] is entry
        threshold_lower = round(
            selected_hist_param.unapply_widget_lower.getValue(), self.round_ndigits
        )
        if threshold_lower < round(selected_hist_param.val_min, self.round_ndigits):
            threshold_lower = round(selected_hist_param.val_min, self.round_ndigits)
            selected_hist_param.unapply_widget_lower.setValue(threshold_lower)
        elif threshold_lower > round(selected_hist_param.val_max, self.round_ndigits):
            threshold_lower = round(selected_hist_param.val_max, self.round_ndigits)
            selected_hist_param.unapply_widget_lower.setValue(threshold_lower)
        # else: # Do nothing

        # now set the new threshold
        selected_hist_param.unapply_threshold_lower = threshold_lower

        self.needredisp = True

    def newThresholdUpper(self):
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key, [1] is entry
        threshold_upper = round(
            selected_hist_param.unapply_widget_upper.getValue(), self.round_ndigits
        )
        if threshold_upper < round(selected_hist_param.val_min, self.round_ndigits):
            threshold_upper = round(selected_hist_param.val_min, self.round_ndigits)
            selected_hist_param.unapply_widget_upper.setValue(threshold_upper)
        elif threshold_upper > round(selected_hist_param.val_max, self.round_ndigits):
            threshold_upper = round(selected_hist_param.val_max, self.round_ndigits)
            selected_hist_param.unapply_widget_upper.setValue(threshold_upper)
        # else: # Do nothing

        # now set the new threshold
        selected_hist_param.unapply_threshold_upper = threshold_upper

        self.needredisp = True

    def newHistogramRow(self, currow):
        "called when a new row is selected from the Histogram list box"

        if self.curhistidx == currow:
            return

        idx_cter = list(self.cter_params_sort.items())[currow][
            1
        ].idx_cter  # cter_params_sort is an OrderedDict, [0] is key, [1] is entry
        param_name = list(self.cter_params_sort.items())[currow][
            1
        ].param_name  # cter_params_sort is an OrderedDict, [0] is key, [0] is key, [1] is entry
        selected_hist_param = list(self.hist_params.items())[currow][
            1
        ]  # hist_params is an OrderedDict, [0] is key, [0] is key, [1] is entry

        # Disable old item
        if self.curthresholdcontrol == self.threshold_control_map["lower"].index:
            selected_hist_param.unapply_widget_lower.setEnabled(False)
        elif self.curthresholdcontrol == self.threshold_control_map["upper"].index:
            selected_hist_param.unapply_widget_upper.setEnabled(False)
        else:
            assert (
                self.curthresholdcontrol
                == self.threshold_control_map["edit_only"].index
            ), "MRK_DEBUG"
            selected_hist_param.unapply_widget_lower.setEnabled(False)
            selected_hist_param.unapply_widget_upper.setEnabled(False)

        # now set the new item and enable it
        self.curhistidx = currow

        param_name = selected_hist_param.param_name
        param_label = self.cter_params[param_name].label

        # Check if the all selected parameter values are same
        if selected_hist_param.val_min == selected_hist_param.val_max:
            # self.curhistdisable=True
            self.curhistogramdisplay = False
            if self.whistparam.isVisible():
                self.whistparam.hide()
            if self.wscatterparam.isVisible():
                self.wscatterparam.hide()
            PyQt5.QtWidgets.QMessageBox.information(
                self,
                "Information",
                "All entries have the same selected parameter values (%s). \n\nParameter Histogram & Plot will not be shown."
                % (param_label),
            )
        else:
            if self.curthresholdcontrol == self.threshold_control_map["lower"].index:
                selected_hist_param.unapply_widget_lower.setEnabled(True)
            elif self.curthresholdcontrol == self.threshold_control_map["upper"].index:
                selected_hist_param.unapply_widget_upper.setEnabled(True)
            else:
                assert (
                    self.curthresholdcontrol
                    == self.threshold_control_map["edit_only"].index
                ), "MRK_DEBUG"
                selected_hist_param.unapply_widget_lower.setEnabled(True)
                selected_hist_param.unapply_widget_upper.setEnabled(True)

            idx_cter = selected_hist_param.idx_cter

            self.whistparam.setWindowTitle("sp_gui_cter - %s Histogram" % (param_label))

            if self.cursyncsort == True:
                idx_sort = self.cter_params_sort[param_name].idx_sort

                if idx_sort != self.cursortidx:
                    self.newSort(idx_sort)
                    self.ssort.setCurrentIndex(idx_sort)
                # else: assert idx_sort == self.cursortidx, "MRK_DEBUG" # Do nothing
            # else: assert self.cursyncsort == False, "MRK_DEBUG" # Do nothing

            if self.cter_partres_file_path == None:
                return  # no cter ctf file is selected
            if self.cter_entry_list == None:
                return  # no cter ctf file is selected

            # self.curhistdisable=False

            # NOTE: 2016/01/03 Toshio Moriya
            # Force update related plots for scaling delay...
            self.updateHist()
            self.updatePlotParam()
            if self.cursyncsort == True:
                self.updateEntryList()

            self.needredisp = True

    def newThresholdControl(self, currow):
        "called when a new row is selected from the Threshold Control list box"

        if self.curthresholdcontrol == currow:
            return

        # Disable old item
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key
        if self.curthresholdcontrol == self.threshold_control_map["lower"].index:
            selected_hist_param.unapply_widget_lower.setEnabled(False)
        elif self.curthresholdcontrol == self.threshold_control_map["upper"].index:
            selected_hist_param.unapply_widget_upper.setEnabled(False)
        else:
            assert (
                self.curthresholdcontrol
                == self.threshold_control_map["edit_only"].index
            ), "MRK_DEBUG"
            selected_hist_param.unapply_widget_lower.setEnabled(False)
            selected_hist_param.unapply_widget_upper.setEnabled(False)

        # now set the new item and enalble it
        self.curthresholdcontrol = currow

        if self.curthresholdcontrol == self.threshold_control_map["lower"].index:
            selected_hist_param.unapply_widget_lower.setEnabled(True)
        elif self.curthresholdcontrol == self.threshold_control_map["upper"].index:
            selected_hist_param.unapply_widget_upper.setEnabled(True)
        else:
            assert (
                self.curthresholdcontrol
                == self.threshold_control_map["edit_only"].index
            ), "MRK_DEBUG"
            selected_hist_param.unapply_widget_lower.setEnabled(True)
            selected_hist_param.unapply_widget_upper.setEnabled(True)

    def newSyncSort(self, syncsort):
        """Change sync sort enable state."""
        if self.cursyncsort == syncsort:
            return

        # now set the new status
        self.cursyncsort = syncsort
        self.ssort.setEnabled(not self.cursyncsort)

        if self.cursyncsort == True:
            # Get sort index
            hist_param = list(self.hist_params.items())[self.curhistidx][1]
            hist_param_name = hist_param.param_name
            idx_sort = self.cter_params_sort[hist_param_name].idx_sort

            if idx_sort != self.cursortidx:
                self.newSort(idx_sort)
                self.ssort.setCurrentIndex(idx_sort)
            # else: assert idx_sort == self.cursortidx, "MRK_DEBUG" # Do nothing
        # else: assert self.cursyncsort == False, "MRK_DEBUG" # Do nothing

    def newEntryPerBin(self, curentryperbin):
        if self.curentryperbin == curentryperbin:
            return

        # now set the new entry per bin
        self.curentryperbin = curentryperbin

        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        # NOTE: 2016/01/03 Toshio Moriya
        # Force update related plots for scaling delay...
        self.updateHist()
        self.updatePlotParam()

        self.needredisp = True

    def applyAllThresholds(self, item=None):
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        reply = PyQt5.QtWidgets.QMessageBox.question(
            self,
            "Warning",
            "Applying all threshold setting will wipe the previous selection states including manual setting. Do you really want to continue?",
            PyQt5.QtWidgets.QMessageBox.Yes | PyQt5.QtWidgets.QMessageBox.No,
            PyQt5.QtWidgets.QMessageBox.No,
        )
        if reply == PyQt5.QtWidgets.QMessageBox.No:
            return

        # Set set the select status of all cter entries based on the threshold values
        for cter_entry in self.cter_entry_list:
            new_select_state = 1
            for idx_hist, (param, curr_hist_params) in enumerate(
                self.hist_params.items()
            ):
                threshold_lower = round(
                    curr_hist_params.unapply_threshold_lower, self.round_ndigits
                )
                threshold_upper = round(
                    curr_hist_params.unapply_threshold_upper, self.round_ndigits
                )
                curr_hist_params.apply_threshold_lower = threshold_lower
                curr_hist_params.apply_threshold_upper = threshold_upper
                curr_hist_params.apply_widget_lower.setValue(threshold_lower)
                curr_hist_params.apply_widget_upper.setValue(threshold_upper)
                idx_cter = curr_hist_params.idx_cter

                param_val = round(cter_entry[idx_cter], self.round_ndigits)
                if param_val < threshold_lower or threshold_upper < param_val:
                    # print "MRK_DEBUG: Param #%d diselected entry #%04d with (param_val, threshold_lower, threshold_upper) = (%1.15g, %1.15g, %1.15g)" % (idx_hist, idx_cter, param_val, threshold_lower, threshold_upper)
                    new_select_state = 0
                # else: # Do nothing
            cter_entry[self.cter_params["select"].idx_cter] = new_select_state

        self.updateEntryList()
        self.updateUncheckCounts()

    def newThresholdSet(self, currow):
        "called when a new row is selected from the Threshold Set list box"

        if self.curthresholdset == currow:
            return
        # now set the new item and enalble it
        self.curthresholdset = currow

    def writeThresholdSet(self, file_path_out, idx_thresholdset):
        assert self.cter_partres_file_path != None, "MRK_DEBUG"
        assert self.cter_entry_list != None, "MRK_DEBUG"

        file_out = open(file_path_out, "w")

        # Write lines to check consistency upon loading
        file_out.write("# @@@@@ gui_cter thresholds - ")
        file_out.write(
            EMAN2_meta.EMANVERSION + " (GITHUB: " + EMAN2_meta.DATESTAMP + ")"
        )
        file_out.write(" @@@@@ \n")
        file_out.write(
            "# Associated CTER Partres File == %s\n" % (self.cter_partres_file_path)
        )
        file_out.write(
            "# Saved Threshold Set == %s\n"
            % (list(self.threshold_status_labels.items())[idx_thresholdset][1])
        )
        file_out.write(
            "# [Parameter Id] [Parameter Name] [Lower Threshold] [Upper Threshold]\n"
        )

        for idx_hist, (param, map_entry) in enumerate(self.hist_params.items()):
            idx_cter = map_entry.idx_cter
            param_label = map_entry.label

            if idx_thresholdset == self.idx_thresholdset_applied:
                threshold_lower = map_entry.apply_threshold_lower
                threshold_upper = map_entry.apply_threshold_upper
            else:
                threshold_lower = map_entry.unapply_threshold_lower
                threshold_upper = map_entry.unapply_threshold_upper

            # NOTE: 2016/01/26 Toshio Moriya
            # Use the precision for double to minimise precision loss by save & load operations
            file_out.write(
                "%2d %s == %1.15g %1.15g \n"
                % (
                    idx_hist,
                    param_label,
                    round(threshold_lower, self.round_ndigits),
                    round(threshold_upper, self.round_ndigits),
                )
            )

        file_out.close()

    def readThresholdSet(self, file_path_in, idx_thresholdset):
        assert self.cter_partres_file_path != None, "MRK_DEBUG"
        assert self.cter_entry_list != None, "MRK_DEBUG"

        file_in = open(file_path_in, "r")

        # Check if this parameter file is threshold
        line_in = file_in.readline()
        if line_in.find("@@@@@ gui_cter thresholds") != -1:
            # loop through the rest of lines
            for line_in in file_in:
                if line_in[0] == "#":
                    continue

                tokens_in = line_in.split("==")
                assert len(tokens_in) == 2, "MRK_DEBUG"
                tokens_label = tokens_in[0].split()
                assert len(tokens_label) >= 2, "MRK_DEBUG"
                idx_hist = int(tokens_label[0])
                map_entry = list(self.hist_params.items())[idx_hist][
                    1
                ]  # [0] is key, [1] is element

                tokens_val = tokens_in[1].split()
                assert len(tokens_val) == 2, "MRK_DEBUG"
                threshold_lower = round(float(tokens_val[0]), self.round_ndigits)
                threshold_upper = round(float(tokens_val[1]), self.round_ndigits)
                map_entry.unapply_threshold_lower = threshold_lower
                map_entry.unapply_threshold_upper = threshold_upper
                map_entry.unapply_widget_lower.setValue(threshold_lower)
                map_entry.unapply_widget_upper.setValue(threshold_upper)

                self.newThresholdLower()
                self.newThresholdUpper()

            if idx_thresholdset == self.idx_thresholdset_applied:
                self.applyAllThresholds()
        else:
            PyQt5.QtWidgets.QMessageBox.warning(
                self, "Warning", "The specified file is not threshold file."
            )

        file_in.close()

    def saveThresholdSet(self, item=None):
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        title_string = (
            "Save %s Thresholds"
            % list(self.threshold_status_labels.items())[self.curthresholdset][1]
        )
        file_path_out = PyQt5.QtWidgets.QFileDialog.getSaveFileName(
            self, title_string, options=PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        )

        # QFileDialog gives a tuple in Qt5, w/unicode elements
        if isinstance(file_path_out, tuple):
            file_path_out = file_path_out[0]
        else:
            file_path_out = str(file_path_out)

        if file_path_out == "":
            return

        self.writeThresholdSet(os.path.relpath(file_path_out), self.curthresholdset)

    def loadThresholdSet(self, item=None):
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        reply = PyQt5.QtWidgets.QMessageBox.question(
            self,
            "Warning",
            "Loading thresholds will wipe the previous threshold setting. Do you really want to continue?",
            PyQt5.QtWidgets.QMessageBox.Yes | PyQt5.QtWidgets.QMessageBox.No,
            PyQt5.QtWidgets.QMessageBox.No,
        )
        if reply == PyQt5.QtWidgets.QMessageBox.No:
            return

        title_string = (
            "Load %s Thresholds"
            % list(self.threshold_status_labels.items())[self.curthresholdset][1]
        )
        file_path_in = PyQt5.QtWidgets.QFileDialog.getOpenFileName(
            self, title_string, options=PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        )

        # QFileDialog gives a tuple in Qt5, w/unicode elements
        if isinstance(file_path_in, tuple):
            file_path_in = file_path_in[0]
        else:
            file_path_in = str(file_path_in)

        if file_path_in == "":
            return

        self.readThresholdSet(os.path.relpath(file_path_in), self.curthresholdset)

    def saveSelection(self, item=None):
        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        assert (
            os.path.basename(self.cter_partres_file_path).find("partres") != -1
        ), "MRK_DEBUG"
        assert self.cter_partres_file_path[-1 * len(".txt") :] == ".txt", "MRK_DEBUG"

        file_prefix = self.vfileprefix.getValue()
        file_path_out_select = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            "%s_partres_select.txt" % (file_prefix),
        )
        file_path_out_discard = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            "%s_partres_discard.txt" % (file_prefix),
        )
        file_path_out_mic_select = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            "%s_micrographs_select.txt" % (file_prefix),
        )
        file_path_out_mic_discard = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            "%s_micrographs_discard.txt" % (file_prefix),
        )
        file_path_out_thresholds = os.path.join(
            os.path.dirname(self.cter_partres_file_path),
            "%s_thresholds.txt" % (file_prefix),
        )

        existing_file_path = None
        if os.path.exists(file_path_out_select):
            existing_file_path = file_path_out_select
        elif os.path.exists(file_path_out_discard):
            existing_file_path = file_path_out_discard
        elif os.path.exists(file_path_out_mic_select):
            existing_file_path = file_path_out_mic_select
        elif os.path.exists(file_path_out_mic_discard):
            existing_file_path = file_path_out_mic_discard
        elif os.path.exists(file_path_out_thresholds):
            existing_file_path = file_path_out_thresholds
        # else: # Do nothing

        if existing_file_path != None:
            reply = PyQt5.QtWidgets.QMessageBox.question(
                self,
                "Warning",
                "The file (%s) already exists. Do you want to overwrite the file?"
                % (existing_file_path),
                PyQt5.QtWidgets.QMessageBox.Yes | PyQt5.QtWidgets.QMessageBox.No,
                PyQt5.QtWidgets.QMessageBox.No,
            )
            if reply == PyQt5.QtWidgets.QMessageBox.No:
                return

        # Save selection in CTER Format
        file_out_select = open(file_path_out_select, "w")
        file_out_discard = open(file_path_out_discard, "w")

        save_cter_entry_list = sorted(
            self.cter_entry_list, key=lambda x: x[self.cter_params["id"].idx_cter]
        )
        idx_cter_ignore_list = [
            self.cter_params["id"].idx_cter,
            self.cter_params["select"].idx_cter,
        ]

        for cter_entry in save_cter_entry_list:
            file_out = file_out_select
            if cter_entry[self.cter_params["select"].idx_cter] == 0:
                file_out = file_out_discard
            # else: assert cter_entry[self.idx_cter_select] == 1, "MRK_DEBUG" # do nothing

            for idx_cter in range(self.n_idx_cter):
                if idx_cter in idx_cter_ignore_list:
                    # Do nothing
                    continue
                elif idx_cter in [self.cter_params["mic_name"].idx_cter]:
                    file_out.write("  %s" % cter_entry[idx_cter])
                else:
                    file_out.write("  %12.5g" % cter_entry[idx_cter])
            file_out.write("\n")

        file_out_select.close()
        file_out_discard.close()

        # Save selection in Micrograph List Format
        file_out_mic_select = open(file_path_out_mic_select, "w")
        file_out_mic_discard = open(file_path_out_mic_discard, "w")

        for cter_entry in save_cter_entry_list:
            file_out = file_out_mic_select
            if cter_entry[self.cter_params["select"].idx_cter] == 0:
                file_out = file_out_mic_discard
            file_out.write(
                "  %s\n"
                % os.path.basename(cter_entry[self.cter_params["mic_name"].idx_cter])
            )

        file_out_mic_select.close()
        file_out_mic_discard.close()

        # Save the associated applied threshold
        self.writeThresholdSet(file_path_out_thresholds, self.idx_thresholdset_applied)

        PyQt5.QtWidgets.QMessageBox.information(
            self,
            "Information",
            "The following files are saved in %s:\n\nCTER Partres File - Selected: %s\n\nCTER Partres FIle - Discarded: %s\n\nMicrograph - Selected: %s\n\nMicrograph - Discarded: %s\n\nApplied Threshold Set: %s"
            % (
                os.path.dirname(self.cter_partres_file_path),
                os.path.basename(file_path_out_select),
                os.path.basename(file_path_out_discard),
                os.path.basename(file_path_out_mic_select),
                os.path.basename(file_path_out_mic_discard),
                os.path.basename(file_path_out_thresholds),
            ),
        )

    def timeOut(self):
        if self.busy:
            return

        # Redisplay before spawning thread for more interactive display
        if self.needredisp:
            try:
                self.redisplay()
            except:
                print("Received unexpected exception from redisplay() in timeOut(): ")
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_traceback)
                # MRK_NOTE: 2015/12/17 Toshio Moriya
                # Another way to print out exception info...
                # lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
                # print ''.join('!! ' + line for line in lines)
                pass

    def redisplay(self):
        self.needredisp = False
        self.busy = True

        if self.cter_entry_list != None:
            is_child_shown = False
            if not self.wfft.isVisible() and self.is_wfft_minimized:
                self.wfft.show()
                is_child_shown = True
            if (
                not self.wimgmicthumb.isVisible()
                and self.curimgmicthumbdisplay
                and not self.is_wimgmicthumb_minimized
                and os.path.exists(self.cter_micthumb_file_path)
            ):
                self.wimgmicthumb.show()
                is_child_shown = True
            if (
                not self.wplotrotavgcoarse.isVisible()
                and self.curplotrotavgdisplay
                and not self.is_wplotrotavgcoarse_minimized
                and os.path.exists(self.cter_pwrot_file_path)
            ):
                self.wplotrotavgcoarse.show()
                is_child_shown = True
            if (
                not self.wplotrotavgfine.isVisible()
                and self.curplotrotzoomdisplay
                and not self.is_wplotrotavgfine_minimized
                and os.path.exists(self.cter_pwrot_file_path)
            ):
                self.wplotrotavgfine.show()
                is_child_shown = True
            # if not self.whistparam.isVisible() and not self.curhistdisable and not self.is_whistparam_minimized:
            if (
                not self.whistparam.isVisible()
                and self.curhistogramdisplay
                and not self.is_whistparam_minimized
            ):
                self.whistparam.show()
                is_child_shown = True
            # if not self.wscatterparam.isVisible() and not self.curhistdisable and not self.is_wscatterparam_minimized:
            if (
                not self.wscatterparam.isVisible()
                and self.curscatterdisplay
                and not self.is_wscatterparam_minimized
            ):
                self.wscatterparam.show()
                is_child_shown = True
            if is_child_shown == True:
                self.raise_()
                self.activateWindow()

        self.updateImgMicThumb()
        self.updateHist()
        self.updatePlotParam()
        self.updatePlotCurves()
        self.updateFFT()

        self.busy = False

    def eventFilter(self, source, event):
        if event.type() == PyQt5.QtCore.QEvent.WindowStateChange:
            if self.windowState() & PyQt5.QtCore.Qt.WindowMinimized:
                assert self.isMinimized() == True, "MRK_DEBUG"
                #
                # NOTE: 2016/03/09 Toshio Moriya
                # Minimize icon button of child window should be disabled
                #
                if self.cter_entry_list != None:
                    if self.wfft.isVisible() and not self.is_wfft_minimized:
                        self.wfft.hide()
                        self.is_wfft_minimized = True
                    if (
                        self.wimgmicthumb.isVisible()
                        and not self.is_wimgmicthumb_minimized
                    ):
                        assert self.curimgmicthumbdisplay, "MRK_DEBUG"
                        self.wimgmicthumb.hide()
                        self.is_wimgmicthumb_minimized = True
                    if (
                        self.wplotrotavgcoarse.isVisible() == True
                        and not self.is_wplotrotavgcoarse_minimized
                    ):
                        assert self.curplotrotavgdisplay, "MRK_DEBUG"
                        self.wplotrotavgcoarse.hide()
                        self.is_wplotrotavgcoarse_minimized = True
                    if (
                        self.wplotrotavgfine.isVisible()
                        and not self.is_wplotrotavgfine_minimized
                    ):
                        assert self.curplotrotzoomdisplay, "MRK_DEBUG"
                        self.wplotrotavgfine.hide()
                        self.is_wplotrotavgfine_minimized = True
                    if (
                        self.whistparam.isVisible() == True
                        and self.is_whistparam_minimized == False
                    ):
                        # assert self.curhistdisable == False, "MRK_DEBUG"
                        assert self.curhistogramdisplay == True, "MRK_DEBUG"
                        self.whistparam.hide()
                        self.is_whistparam_minimized = True
                    if (
                        self.wscatterparam.isVisible() == True
                        and self.is_wscatterparam_minimized == False
                    ):
                        # assert self.curhistdisable == False, "MRK_DEBUG"
                        assert self.curhistogramdisplay == True, "MRK_DEBUG"
                        self.wscatterparam.hide()
                        self.is_wscatterparam_minimized = True
            else:
                assert self.isMinimized() == False, "MRK_DEBUG"
                #
                # NOTE: 2016/03/09 Toshio Moriya
                # Minimize icon button of child window should be disabled
                if self.cter_entry_list != None:
                    if self.is_wfft_minimized == True:
                        assert (
                            self.wfft.isVisible() == False
                            and self.curfftdisplay == True
                        ), "MRK_DEBUG"
                        self.wfft.show()
                        self.is_wfft_minimized = False
                    if self.is_wimgmicthumb_minimized:
                        assert (
                            not self.wimgmicthumb.isVisible()
                            and self.curimgmicthumbdisplay
                        ), "MRK_DEBUG"
                        self.wimgmicthumb.show()
                        self.is_wimgmicthumb_minimized = False
                    if self.is_wplotrotavgcoarse_minimized:
                        assert (
                            not self.wplotrotavgcoarse.isVisible()
                            and self.curplotrotavgdisplay
                        ), "MRK_DEBUG"
                        self.wplotrotavgcoarse.show()
                        self.is_wplotrotavgcoarse_minimized = False
                    if self.is_wplotrotavgfine_minimized:
                        assert (
                            not self.wplotrotavgfine.isVisible()
                            and self.curplotrotzoomdisplay
                        ), "MRK_DEBUG"
                        self.wplotrotavgfine.show()
                        self.is_wplotrotavgfine_minimized = False
                    if self.is_whistparam_minimized:
                        # assert not self.whistparam.isVisible() and not self.curhistdisable, "MRK_DEBUG"
                        assert (
                            not self.whistparam.isVisible() and self.curhistogramdisplay
                        ), "MRK_DEBUG"
                        self.whistparam.show()
                        self.is_whistparam_minimized = False
                    if self.is_wscatterparam_minimized:
                        # assert not self.wscatterparam.isVisible() and not self.curhistdisable, "MRK_DEBUG"
                        assert (
                            not self.wscatterparam.isVisible()
                            and self.curhistogramdisplay
                        ), "MRK_DEBUG"
                        self.wscatterparam.show()
                        self.is_wscatterparam_minimized = False
                if self.isVisible():  # Depends on timing at startup, this can happen?!!
                    self.raise_()
                    self.activateWindow()
        # If focus-follow-mouse set in window manager, the pop-ups flicker if they're raised automatically. --Tapu, 2019/04/03
        # elif event.type() == QtCore.QEvent.WindowActivate:
        ## print "MRK_DEBUG: sxgui main window has gained focus (become active)"
        ##
        ## NOTE: 2016/03/08 Toshio Moriya
        ## To raise EMGLWidget (SXPlot2DWidget and EMImage2DWidget) window,
        ## we have to go through its qt_parent attribute to call raise_()...
        ##
        # if self.cter_entry_list != None:
        # if self.wfft.isVisible() == True:
        # self.wfft.qt_parent.raise_()
        # if self.wimgmicthumb.isVisible() == True:
        # self.wimgmicthumb.qt_parent.raise_()
        # if self.wplotrotavgcoarse.isVisible() == True:
        # self.wplotrotavgcoarse.qt_parent.raise_()
        # if self.wplotrotavgfine.isVisible() == True:
        # self.wplotrotavgfine.qt_parent.raise_()
        # if self.whistparam.isVisible() == True:
        # self.whistparam.qt_parent.raise_()
        # if self.wscatterparam.isVisible() == True:
        # self.wscatterparam.qt_parent.raise_()
        # assert self.isVisible(), "MRK_DEBUG"
        # self.raise_()

        return super(SXGuiCter, self).eventFilter(source, event)

    def closeEvent(self, event):
        # Save current window layout
        EMAN2.E2saveappwin("sp_gui_cter", "main", self)
        if self.cter_entry_list != None:
            EMAN2.E2saveappwin("sp_gui_cter", "fft", self.wfft.qt_parent)
            EMAN2.E2saveappwin(
                "sp_gui_cter", "imgmicthumb", self.wimgmicthumb.qt_parent
            )
            EMAN2.E2saveappwin(
                "sp_gui_cter", "plotcoarse", self.wplotrotavgcoarse.qt_parent
            )
            EMAN2.E2saveappwin(
                "sp_gui_cter", "plotfine", self.wplotrotavgfine.qt_parent
            )
            EMAN2.E2saveappwin("sp_gui_cter", "plotparam", self.wscatterparam.qt_parent)
            EMAN2.E2saveappwin("sp_gui_cter", "histparam", self.whistparam.qt_parent)

        # close all child windows
        if self.wfft:
            self.wfft.close()
        if self.wimgmicthumb:
            self.wimgmicthumb.close()
        if self.wplotrotavgcoarse:
            self.wplotrotavgcoarse.close()
        if self.wplotrotavgfine:
            self.wplotrotavgfine.close()
        if self.whistparam:
            self.whistparam.close()
        if self.wscatterparam:
            self.wscatterparam.close()

        event.accept()
        # NOTE: Toshio Moriya 2017/11/23
        # The following code seems to be printing "None" to standard output
        # upon exiting the application. But why???
        PyQt5.QtWidgets.qApp.exit(0)

    def updatePlotVisibility(self, val=None):
        if self.wplotrotavgcoarse == None:
            return  # it's closed/not visible
        if self.wplotrotavgfine == None:
            return  # it's closed/not visible
        if self.cter_pwrot_file_path == None:
            return  # no cter entry is selected

        for item in self.graph_map:
            item_widget = self.graph_map[item].item_widget
            name = self.graph_map[item].item_name

            if self.wplotrotavgfine.visibility[name] != item_widget.getValue():
                self.wplotrotavgfine.visibility[name] = item_widget.getValue()
                self.wplotrotavgfine.full_refresh()
                self.wplotrotavgfine.updateGL()

            if self.wplotrotavgcoarse.visibility[name] != item_widget.getValue():
                self.wplotrotavgcoarse.visibility[name] = item_widget.getValue()
                self.wplotrotavgcoarse.full_refresh()
                self.wplotrotavgcoarse.updateGL()

    def newPlotFixScale(self, curplotfixscale):
        if self.curplotfixscale == curplotfixscale:
            return

        # now set the new entry per bin
        self.curplotfixscale = curplotfixscale

        if self.cter_partres_file_path == None:
            return  # no cter ctf file is selected
        if self.cter_entry_list == None:
            return  # no cter ctf file is selected

        self.needredisp = True

    def refreshGraphs(self, item):
        self.needredisp = True

    def plotparammouseup(self, event):
        if self.curthresholdcontrol == self.threshold_control_map["edit_only"].index:
            return

        # Swap control if shift button is pressed
        is_not_reverse_control = True
        modifiers = event.modifiers()
        if modifiers & PyQt5.QtCore.Qt.ShiftModifier:
            is_not_reverse_control = False

        plot_x, plot_y = self.wscatterparam.scr2plot(event.x(), event.y())
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key

        if self.curthresholdcontrol == self.threshold_control_map["lower"].index:
            if is_not_reverse_control:
                selected_hist_param.unapply_widget_lower.setValue(plot_y)
                self.newThresholdLower()
            else:
                selected_hist_param.unapply_widget_upper.setValue(plot_y)
                self.newThresholdUpper()
        else:
            assert (
                self.curthresholdcontrol == self.threshold_control_map["upper"].index
            ), "MRK_DEBUG"
            if is_not_reverse_control:
                selected_hist_param.unapply_widget_upper.setValue(plot_y)
                self.newThresholdUpper()
            else:
                selected_hist_param.unapply_widget_lower.setValue(plot_y)
                self.newThresholdLower()

    def histparammouseup(self, event):
        if self.curthresholdcontrol == self.threshold_control_map["edit_only"].index:
            return

        # Swap control if shift button is pressed
        is_not_reverse_control = True
        modifiers = event.modifiers()
        if modifiers & PyQt5.QtCore.Qt.ShiftModifier:
            is_not_reverse_control = False

        hist_x, hist_y = self.whistparam.scr2plot(event.x(), event.y())
        selected_hist_param = list(self.hist_params.items())[self.curhistidx][
            1
        ]  # hist_params is an OrderedDict, [0] is key

        if self.curthresholdcontrol == self.threshold_control_map["lower"].index:
            if is_not_reverse_control:
                selected_hist_param.unapply_widget_lower.setValue(hist_x)
                self.newThresholdLower()
            else:
                selected_hist_param.unapply_widget_upper.setValue(hist_x)
                self.newThresholdUpper()

        else:
            assert (
                self.curthresholdcontrol == self.threshold_control_map["upper"].index
            ), "MRK_DEBUG"
            if is_not_reverse_control:
                selected_hist_param.unapply_widget_upper.setValue(hist_x)
                self.newThresholdUpper()
            else:
                selected_hist_param.unapply_widget_lower.setValue(hist_x)
                self.newThresholdLower()

def main():
    sp_global_def.print_timestamp("Start")
    sp_global_def.write_command()
    run()
    sp_global_def.print_timestamp("Finish")

if __name__ == "__main__":
    main()
