#!/usr/bin/env python
# sxgui_unblur for analyzing drift parameters made by Unblur and MotionCor2
# Copyright (C) 2016  Markus Stabrin (markus.stabrin@mpi-dortmund.mpg.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from matplotlib import pylab
try:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvasQTAgg
except ImportError:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQT as FigureCanvasQTAgg
try:
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar2QTAgg
except ImportError:
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar2QTAgg
import os
import sys
import glob
import numpy

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MSMainWidget(object):
    def setupUi(self, MSMainWidget):
        widget = QtGui.QWidget(self)
        MSMainWidget.setCentralWidget(widget)

        self.global_layout = QtGui.QVBoxLayout(widget)
        self.h_layout = QtGui.QHBoxLayout()
        self.v_layout_1 = QtGui.QVBoxLayout()
        self.v_layout_2 = QtGui.QVBoxLayout()
        self.v_layout_3 = QtGui.QVBoxLayout()

        self.global_layout.addLayout(self.h_layout, stretch=1)
        self.h_layout.addLayout(self.v_layout_1, stretch=1)
        self.h_layout.addLayout(self.v_layout_2, stretch=0)
        self.h_layout.addLayout(self.v_layout_3, stretch=0)

        self.v_layout_1.addWidget(self.load_section())
        self.v_layout_1.addWidget(self.list_section())
        self.v_layout_2.addWidget(self.status_section())
        self.v_layout_2.addStretch(1)
        self.v_layout_2.addWidget(self.sort_section())
        self.v_layout_3.addWidget(self.plot_section())
        self.v_layout_3.addStretch(1)
        self.v_layout_3.addWidget(self.threshold_section())
        self.global_layout.addWidget(self.save_section())

    def load_section(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 5, 5, 5)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        layout_v.addLayout(layout_h)
        
        label = QtGui.QLabel('Pattern:', self)
        self.lePattern = QtGui.QLineEdit(self)
        self.pbSelectTxt = QtGui.QPushButton('Find pattern', self)
        layout_h.addWidget(label)
        layout_h.addWidget(self.lePattern)
        layout_h.addWidget(self.pbSelectTxt)

        self.pbImportPattern = QtGui.QPushButton('Load files by pattern', self)
        layout_v.addWidget(self.pbImportPattern)
        self.pbImportList = QtGui.QPushButton('Load files by list', self)
        layout_v.addWidget(self.pbImportList)
        return widget

    def list_section(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 5, 5, 5)

        self.lsFiles = QtGui.QListWidget(self)
        self.lsFiles.setEnabled(False)
        layout_v.addWidget(self.lsFiles)
        return widget

    def status_section(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 5, 5, 5)

        label = QtGui.QLabel('Info of current entry:', self)
        layout_v.addWidget(label)
        label = QtGui.QLabel('', self)
        layout_v.addWidget(label)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Micrograph name', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leCurrentMicName = QtGui.QLineEdit(self)
        self.leCurrentMicName.setMaximumWidth(300)
        self.leCurrentMicName.setMinimumWidth(300)
        self.leCurrentMicName.setEnabled(False)
        layout_h.addWidget(self.leCurrentMicName)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Overall drift [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leCurrentOverallDrift = QtGui.QLineEdit(self)
        self.leCurrentOverallDrift.setMaximumWidth(100)
        self.leCurrentOverallDrift.setEnabled(False)
        layout_h.addWidget(self.leCurrentOverallDrift)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Drift per frame [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leCurrentFrameDrift = QtGui.QLineEdit(self)
        self.leCurrentFrameDrift.setMaximumWidth(100)
        self.leCurrentFrameDrift.setEnabled(False)
        layout_h.addWidget(self.leCurrentFrameDrift)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('End to end length [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leCurrentEndToEndDrift = QtGui.QLineEdit(self)
        self.leCurrentEndToEndDrift.setMaximumWidth(100)
        self.leCurrentEndToEndDrift.setEnabled(False)
        layout_h.addWidget(self.leCurrentEndToEndDrift)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Maximum distance between frames [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leCurrentMaxDistance = QtGui.QLineEdit(self)
        self.leCurrentMaxDistance.setMaximumWidth(100)
        self.leCurrentMaxDistance.setEnabled(False)
        layout_h.addWidget(self.leCurrentMaxDistance)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Maximum distance from start frame [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leCurrentMaxDistanceZero = QtGui.QLineEdit(self)
        self.leCurrentMaxDistanceZero.setMaximumWidth(100)
        self.leCurrentMaxDistanceZero.setEnabled(False)
        layout_h.addWidget(self.leCurrentMaxDistanceZero)
        layout_v.addLayout(layout_h)

        label = QtGui.QLabel('', self)
        layout_v.addWidget(label)
        label = QtGui.QLabel('Drift info for selected micrographs:', self)
        layout_v.addWidget(label)
        label = QtGui.QLabel('', self)
        layout_v.addWidget(label)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Nr. of micrographs', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leAllMicNumber = QtGui.QLineEdit(self)
        self.leAllMicNumber.setMaximumWidth(100)
        self.leAllMicNumber.setEnabled(False)
        layout_h.addWidget(self.leAllMicNumber)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Checked micrographs', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leMicChecked = QtGui.QLineEdit(self)
        self.leMicChecked.setMaximumWidth(100)
        self.leMicChecked.setEnabled(False)
        layout_h.addWidget(self.leMicChecked)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Average overall drift [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leAllOverallDrift = QtGui.QLineEdit(self)
        self.leAllOverallDrift.setMaximumWidth(100)
        self.leAllOverallDrift.setEnabled(False)
        layout_h.addWidget(self.leAllOverallDrift)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Average drift per frame [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leAllFrameDrift = QtGui.QLineEdit(self)
        self.leAllFrameDrift.setMaximumWidth(100)
        self.leAllFrameDrift.setEnabled(False)
        layout_h.addWidget(self.leAllFrameDrift)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Average end to end length [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leAllEndToEndDrift = QtGui.QLineEdit(self)
        self.leAllEndToEndDrift.setMaximumWidth(100)
        self.leAllEndToEndDrift.setEnabled(False)
        layout_h.addWidget(self.leAllEndToEndDrift)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Average maximum distance between frames [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leAllMaxDistance = QtGui.QLineEdit(self)
        self.leAllMaxDistance.setMaximumWidth(100)
        self.leAllMaxDistance.setEnabled(False)
        layout_h.addWidget(self.leAllMaxDistance)
        layout_v.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        label = QtGui.QLabel('Average maximum distance from start frame [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leAllMaxDistanceZero = QtGui.QLineEdit(self)
        self.leAllMaxDistanceZero.setMaximumWidth(100)
        self.leAllMaxDistanceZero.setEnabled(False)
        layout_h.addWidget(self.leAllMaxDistanceZero)
        layout_v.addLayout(layout_h)
        
        return widget
    
    def sort_section(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 20, 5, 5)

        label = QtGui.QLabel('Sort entries:', self)
        layout_v.addWidget(label)

        entries = [
            'File name',
            'Overall drift',
            'Drift per frame',
            'End to end length',
            'Maximum distance between frames',
            'Maximum distance from start frame'
            ]
        self.cbSort = QtGui.QComboBox(self)
        self.cbSort.setEnabled(False)
        self.cbSort.addItems(entries)
        layout_v.addWidget(self.cbSort)

        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 0, 0, 0)
        layout_v.addLayout(layout_h)
        self.chDescending = QtGui.QCheckBox('Descending', self)
        self.chDescending.setEnabled(False)
        layout_h.addWidget(self.chDescending)
        self.chSortSelected = QtGui.QCheckBox('Sort selected', self)
        self.chSortSelected.setEnabled(False)
        layout_h.addWidget(self.chSortSelected)

        return widget

    def plot_section(self):
        widget = QtGui.QWidget(self)
        layout_h = QtGui.QHBoxLayout(widget)
        layout_h.setContentsMargins(5, 5, 5, 5)

        layout_v1 = QtGui.QVBoxLayout()
        layout_v1.setContentsMargins(0, 0, 5, 0)
        layout_v2 = QtGui.QVBoxLayout()
        layout_v2.setContentsMargins(5, 0, 0, 0)
        layout_h.addLayout(layout_v1)
        layout_h.addLayout(layout_v2)

        label = QtGui.QLabel('Show plots of current entry:', self)
        layout_v1.addWidget(label)

        self.chPlotDriftMic = QtGui.QCheckBox('Drift', self)
        self.chPlotDriftMic.setEnabled(False)
        layout_v1.addWidget(self.chPlotDriftMic)

        self.chPlotFrameMic = QtGui.QCheckBox('Drift per frame', self)
        self.chPlotFrameMic.setEnabled(False)
        layout_v1.addWidget(self.chPlotFrameMic)

        self.chPlotAngleMic = QtGui.QCheckBox('Angle per frame', self)
        self.chPlotAngleMic.setEnabled(False)
        layout_v1.addWidget(self.chPlotAngleMic)

        label = QtGui.QLabel('Show plots of all micrographs:', self)
        label.setEnabled(False)
        layout_v2.addWidget(label)

        self.chAverageDriftPerFrame = QtGui.QCheckBox('Average drift per frame', self)
        self.chAverageDriftPerFrame.setEnabled(False)
        layout_v2.addWidget(self.chAverageDriftPerFrame)

        self.chPlotDrift = QtGui.QCheckBox('Overall drift histogram', self)
        self.chPlotDrift.setEnabled(False)
        layout_v2.addWidget(self.chPlotDrift)

        self.chPlotFrame = QtGui.QCheckBox('Drift per frame histogram', self)
        self.chPlotFrame.setEnabled(False)
        layout_v2.addWidget(self.chPlotFrame)

        self.chPlotAngle = QtGui.QCheckBox('Angle per frame histogram', self)
        self.chPlotAngle.setEnabled(False)
        layout_v2.addWidget(self.chPlotAngle)

        self.chPlotPerMic = QtGui.QCheckBox('Overall drift per micrograph', self)
        self.chPlotPerMic.setEnabled(False)
        layout_v2.addWidget(self.chPlotPerMic)

        layout_v1.addStretch(1)
        layout_v2.addStretch(1)

        return widget

    def threshold_section(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 20, 5, 5)
        
        layout_h = QtGui.QHBoxLayout()
        layout_h.setContentsMargins(0, 5, 0, 5)
        label = QtGui.QLabel('Start frame:', self)
        layout_h.addWidget(label)
        self.leFrameStart = QtGui.QLineEdit(self)
        self.leFrameStart.setEnabled(False)
        self.leFrameStart.setMaximumWidth(50)
        layout_h.addWidget(self.leFrameStart)
        layout_h.addStretch(1)
        label = QtGui.QLabel('End frame:', self)
        layout_h.addWidget(label)
        self.leFrameStop = QtGui.QLineEdit(self)
        self.leFrameStop.setEnabled(False)
        self.leFrameStop.setMaximumWidth(50)
        layout_h.addWidget(self.leFrameStop)
        layout_v.addLayout(layout_h)
        
        layout_v.addWidget(self.threshold_section_overall())
        layout_v.addWidget(self.threshold_section_frame())
        layout_v.addWidget(self.threshold_section_angle())

        self.pbApply = QtGui.QPushButton('Apply settings marked as criterion', self)
        self.pbApply.setEnabled(False)
        layout_v.addWidget(self.pbApply)
        return widget

    def threshold_section_overall(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 5, 5, 0)

        layout_h1 = QtGui.QHBoxLayout()
        layout_v1 = QtGui.QVBoxLayout()
        layout_v2 = QtGui.QVBoxLayout()

        label = QtGui.QLabel('Threshold overall drift', self)
        layout_v.addWidget(label)
        
        layout_v.addLayout(layout_h1)
        layout_h1.addLayout(layout_v1)
        layout_h1.addLayout(layout_v2)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Start [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartOverall = QtGui.QLineEdit(self)
        self.leStartOverall.setEnabled(False)
        self.leStartOverall.setMaximumWidth(100)
        layout_h.addWidget(self.leStartOverall)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Stop [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopOverall = QtGui.QLineEdit(self)
        self.leStopOverall.setEnabled(False)
        self.leStopOverall.setMaximumWidth(100)
        layout_h.addWidget(self.leStopOverall)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Registered start [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartSaveOverall = QtGui.QLineEdit(self)
        self.leStartSaveOverall.setEnabled(False)
        self.leStartSaveOverall.setMaximumWidth(100)
        layout_h.addWidget(self.leStartSaveOverall)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Registered stop [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopSaveOverall = QtGui.QLineEdit(self)
        self.leStopSaveOverall.setEnabled(False)
        self.leStopSaveOverall.setMaximumWidth(100)
        layout_h.addWidget(self.leStopSaveOverall)
        layout_v1.addLayout(layout_h)

        self.pbSaveOverall = QtGui.QPushButton('Register', self)
        self.pbSaveOverall.setEnabled(False)
        layout_v2.addWidget(self.pbSaveOverall)
        self.chOverallCriterion = QtGui.QCheckBox('Use as criterion', self)
        self.chOverallCriterion.setEnabled(False)
        layout_v2.addWidget(self.chOverallCriterion)
        layout_v2.addStretch(1)

        return widget

    def threshold_section_frame(self):
        tabWidget = QtGui.QTabWidget(self)
        tab_1 = QtGui.QWidget(self)
        tab_2 = QtGui.QWidget(self)

        tabWidget.addTab(tab_1, 'General')
        tabWidget.addTab(tab_2, 'Per frame')

        layout_v = QtGui.QVBoxLayout(tab_1)
        layout_v.setContentsMargins(5, 5, 5, 0)

        layout_h1 = QtGui.QHBoxLayout()
        layout_v1 = QtGui.QVBoxLayout()
        layout_v2 = QtGui.QVBoxLayout()

        label = QtGui.QLabel('Threshold for the drift of every frame', self)
        layout_v.addWidget(label)
        
        layout_v.addLayout(layout_h1)
        layout_h1.addLayout(layout_v1)
        layout_h1.addLayout(layout_v2)
        
        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Start [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartGeneral = QtGui.QLineEdit(self)
        self.leStartGeneral.setEnabled(False)
        self.leStartGeneral.setMaximumWidth(100)
        layout_h.addWidget(self.leStartGeneral)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Stop [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopGeneral = QtGui.QLineEdit(self)
        self.leStopGeneral.setEnabled(False)
        self.leStopGeneral.setMaximumWidth(100)
        layout_h.addWidget(self.leStopGeneral)
        layout_v1.addLayout(layout_h)

        self.pbSaveGeneral = QtGui.QPushButton('Register', self)
        self.pbSaveGeneral.setEnabled(False)
        layout_v2.addWidget(self.pbSaveGeneral)
        self.chGeneralCriterion = QtGui.QCheckBox('Use as criterion', self)
        self.chGeneralCriterion.setEnabled(False)
        layout_v2.addWidget(self.chGeneralCriterion)
        layout_v2.addStretch(1)
        layout_v1.addStretch(1)

        layout_v = QtGui.QVBoxLayout(tab_2)
        layout_v.setContentsMargins(5, 5, 5, 0)

        layout_h1 = QtGui.QHBoxLayout()
        layout_v1 = QtGui.QVBoxLayout()
        layout_v2 = QtGui.QVBoxLayout()

        label = QtGui.QLabel('Threshold drift per frame', self)
        layout_v.addWidget(label)
        
        layout_v.addLayout(layout_h1)
        layout_h1.addLayout(layout_v1)
        layout_h1.addLayout(layout_v2)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Start [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartFrame = QtGui.QLineEdit(self)
        self.leStartFrame.setEnabled(False)
        self.leStartFrame.setMaximumWidth(100)
        layout_h.addWidget(self.leStartFrame)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Stop [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopFrame = QtGui.QLineEdit(self)
        self.leStopFrame.setEnabled(False)
        self.leStopFrame.setMaximumWidth(100)
        layout_h.addWidget(self.leStopFrame)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Registered start [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartSaveFrame = QtGui.QLineEdit(self)
        self.leStartSaveFrame.setEnabled(False)
        self.leStartSaveFrame.setMaximumWidth(100)
        layout_h.addWidget(self.leStartSaveFrame)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Registered stop [A]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopSaveFrame = QtGui.QLineEdit(self)
        self.leStopSaveFrame.setEnabled(False)
        self.leStopSaveFrame.setMaximumWidth(100)
        layout_h.addWidget(self.leStopSaveFrame)
        layout_v1.addLayout(layout_h)

        self.cbFrame = QtGui.QComboBox(self)
        self.cbFrame.setEnabled(False)
        layout_v2.addWidget(self.cbFrame)
        self.pbSaveFrame = QtGui.QPushButton('Register', self)
        self.pbSaveFrame.setEnabled(False)
        layout_v2.addWidget(self.pbSaveFrame)
        self.chFrameCriterion = QtGui.QCheckBox('Use as criterion', self)
        self.chFrameCriterion.setEnabled(False)
        layout_v2.addWidget(self.chFrameCriterion)
        layout_v2.addStretch(1)

        return tabWidget

    def threshold_section_angle(self):
        widget = QtGui.QWidget(self)
        layout_v = QtGui.QVBoxLayout(widget)
        layout_v.setContentsMargins(5, 5, 5, 0)

        layout_h1 = QtGui.QHBoxLayout()
        layout_v1 = QtGui.QVBoxLayout()
        layout_v2 = QtGui.QVBoxLayout()

        label = QtGui.QLabel('Threshold angle', self)
        layout_v.addWidget(label)
        
        layout_v.addLayout(layout_h1)
        layout_h1.addLayout(layout_v1)
        layout_h1.addLayout(layout_v2)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Start [Degree]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartAngle = QtGui.QLineEdit(self)
        self.leStartAngle.setEnabled(False)
        self.leStartAngle.setMaximumWidth(100)
        layout_h.addWidget(self.leStartAngle)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Stop [Degree]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopAngle = QtGui.QLineEdit(self)
        self.leStopAngle.setEnabled(False)
        self.leStopAngle.setMaximumWidth(100)
        layout_h.addWidget(self.leStopAngle)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Registered start [Degree]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStartSaveAngle = QtGui.QLineEdit(self)
        self.leStartSaveAngle.setEnabled(False)
        self.leStartSaveAngle.setMaximumWidth(100)
        layout_h.addWidget(self.leStartSaveAngle)
        layout_v1.addLayout(layout_h)

        layout_h = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Registered stop [Degree]', self)
        layout_h.addWidget(label)
        layout_h.addStretch(1)
        self.leStopSaveAngle = QtGui.QLineEdit(self)
        self.leStopSaveAngle.setEnabled(False)
        self.leStopSaveAngle.setMaximumWidth(100)
        layout_h.addWidget(self.leStopSaveAngle)
        layout_v1.addLayout(layout_h)

        self.cbAngle = QtGui.QComboBox(self)
        self.cbAngle.setEnabled(False)
        layout_v2.addWidget(self.cbAngle)
        self.pbSaveAngle = QtGui.QPushButton('Register', self)
        self.pbSaveAngle.setEnabled(False)
        layout_v2.addWidget(self.pbSaveAngle)
        self.chAngleCriterion = QtGui.QCheckBox('Use as criterion', self)
        self.chAngleCriterion.setEnabled(False)
        layout_v2.addWidget(self.chAngleCriterion)
        self.pbUncheckCriterion = QtGui.QPushButton('Uncheck criteria')
        self.pbUncheckCriterion.setEnabled(False)
        layout_v2.addWidget(self.pbUncheckCriterion)

        return widget

    def save_section(self):
        widget = QtGui.QWidget(self)
        layout_h = QtGui.QHBoxLayout(widget)
        layout_h.setContentsMargins(5, 0, 5, 0)

        label = QtGui.QLabel('Output prefix:', self)
        layout_h.addWidget(label)
        
        self.leOutputName = QtGui.QLineEdit('Trial00', self)
        self.leOutputName.setEnabled(False)
        self.leOutputName.setMinimumWidth(300)
        layout_h.addWidget(self.leOutputName)

        self.pbSaveSelected = QtGui.QPushButton('Select output directory and save selection', self)
        self.pbSaveSelected.setEnabled(False)
        layout_h.addWidget(self.pbSaveSelected)

        layout_h.addStretch(1)

        self.pbSaveSettings = QtGui.QPushButton('Save settings', self)
        self.pbSaveSettings.setEnabled(False)
        layout_h.addWidget(self.pbSaveSettings)

        self.pbLoadSettings = QtGui.QPushButton('Load settings', self)
        layout_h.addWidget(self.pbLoadSettings)

        self.pbAbout = QtGui.QPushButton('About', self)
        layout_h.addWidget(self.pbAbout)

        return widget



class SXUnblurPlot(QtGui.QWidget):

    # Refresh Signal, Frame Changed Signal, Close Signal
    sigRefresh = QtCore.pyqtSignal(list)
    sigFrame = QtCore.pyqtSignal(str)
    sigClose = QtCore.pyqtSignal()

    def __init__(self, title, setframes=False, parent=None):
        super(SXUnblurPlot, self).__init__(parent)

        # Set the Layout of the Widget
        self.setWindowTitle(title)
        self.layout = QtGui.QVBoxLayout(self)

        # Initialise later used Variables
        self.canvas = None
        self.toolbar = None
        self.scrollArea = None

        # Initialise Widgets
        self.widgetPlot = QtGui.QWidget(self)
        self.layoutPlot = QtGui.QVBoxLayout(self.widgetPlot)

        if setframes:
            # Set Variables
            self.dictFrames = {}
            self.widgetFrames = QtGui.QWidget(self)
            self.layoutFrames = QtGui.QVBoxLayout(self.widgetFrames)

            # Add to Layout
            self.layout.addWidget(self.widgetFrames)

        # Add to Layout
        self.layout.addWidget(self.widgetPlot)

        # Remove miminize button
        self.setWindowFlags(QtCore.Qt.WindowCloseButtonHint)

    def _fill_widget_frames_single(self, frame, name):

        # Make space for new frame figures by deleting the old
        # Scroll Area
        if self.scrollArea is not None:
            self.scrollArea.setParent(None)

        # Set Widgets
        self.scrollArea = QtGui.QScrollArea(self.widgetFrames)
        scrollContent = QtGui.QWidget(self.scrollArea)

        # Set Layouts
        layoutContent = QtGui.QVBoxLayout(scrollContent)

        # Fill canvas with the figure and connect it to the click event
        canvas = FigureCanvasQTAgg(frame)
        canvas.mpl_connect('button_press_event', self._change_frame)

        # Canvas Layout settings
        canvas.setMinimumSize(canvas.size())
        layoutContent.addWidget(canvas)

        # Fill dictionary
        self.dictFrames.update({
            '{:s}'.format(canvas): '{:s} {:d}'.format(name, 1)
            })

        # Set Layout
        self.scrollArea.setWidget(scrollContent)
        self.layoutFrames.addWidget(self.scrollArea)
        scrollContent.setLayout(layoutContent)
        self.widgetFrames.setLayout(self.layoutFrames)

    def _fill_widget_frames(self, framesaslist, name, framestart, framestop):
        """Fill the frame widgt with frames"""

        # Make space for new frame figures by deleting the old
        # Scroll Area
        if self.scrollArea is not None:
            self.scrollArea.setParent(None)

        # Set Widgets
        self.scrollArea = QtGui.QScrollArea(self.widgetFrames)
        scrollContent = QtGui.QWidget(self.scrollArea)
        scrollContentUpper = QtGui.QWidget(scrollContent)
        scrollContentLower = QtGui.QWidget(scrollContent)

        # Set Layouts
        layoutScroll = QtGui.QVBoxLayout(scrollContent)
        layoutScrollUpper = QtGui.QHBoxLayout(scrollContentUpper)
        layoutScrollLower = QtGui.QHBoxLayout(scrollContentLower)

        # Calculate how many frames needs to be in the firs line
        intNrFrames = len(framesaslist)
        intLenUpper = intNrFrames // 2
        if intLenUpper % 2 == 0:
            intLenUpper += 1

        # First half of the figure add to the upper widget
        for number in xrange(intLenUpper):
            # Create a Widget for each figure.
            figWidgetUpper = QtGui.QWidget(scrollContentUpper)
            figWidgetLower = QtGui.QWidget(scrollContentLower)
            # Add the figure to the Canvas.
            canvasUpper = FigureCanvasQTAgg(framesaslist[number])
            canvasLower = FigureCanvasQTAgg(
                framesaslist[intNrFrames - 1 - number]
                )
            canvasUpper.setParent(figWidgetLower)
            canvasLower.setParent(figWidgetUpper)
            # Set a pressed Signal to the Canvas.
            canvasUpper.mpl_connect('button_press_event', self._change_frame)
            canvasLower.mpl_connect('button_press_event', self._change_frame)
            # Create a relationship between the name of the Canvas
            # and the Frame number.
            self.dictFrames.update({
                '{:s}'.format(canvasUpper): '{:s} {:d}'.format(
                    name, framestart + number
                    )
                })
            self.dictFrames.update({
                '{:s}'.format(canvasLower): '{:s} {:d}'.format(
                    name, framestop - 1 - number
                    )
                })
            # Create a layout, add the canvas to it and
            # set it to the figWidget.
            figLayoutUpper = QtGui.QVBoxLayout()
            figLayoutLower = QtGui.QVBoxLayout()
            figLayoutUpper.addWidget(canvasUpper)
            figLayoutLower.addWidget(canvasLower)
            figWidgetUpper.setLayout(figLayoutUpper)
            figWidgetLower.setLayout(figLayoutLower)

            # Rescale the Widget and add the Widget to the layoutScrolls
            canvasUpper.setMinimumSize(canvasUpper.size())
            canvasLower.setMinimumSize(canvasLower.size())
            layoutScrollUpper.addWidget(figWidgetUpper)
            layoutScrollLower.addWidget(figWidgetLower)

        # First half of the figure add to the lower widget
        if intNrFrames % 2 != 0:
            # Create a Widget for each figure
            figWidgetUpper = QtGui.QWidget(scrollContentUpper)
            figWidgetLower = QtGui.QWidget(scrollContentLower)
            # Add the figure to the Canvas, funny plot for the missing one
            figFunny = self._funny_plot()
            canvasUpper = FigureCanvasQTAgg(framesaslist[number + 1])
            canvasLower = FigureCanvasQTAgg(
                figFunny
                )
            pylab.close(figFunny)
            canvasUpper.setParent(figWidgetLower)
            canvasLower.setParent(figWidgetUpper)
            # Set a pressed Signal to the Canvas
            canvasUpper.mpl_connect('button_press_event', self._change_frame)
            # Create a relationship between the name of the Canvas
            # and the Frame number.
            self.dictFrames.update({
                '{:s}'.format(canvasUpper): '{:s} {:d}'.format(
                    name, framestart + number + 1
                    )
                })
            # Create a layout, add the canvas to it and
            # set it to the figWidget.
            figLayoutUpper = QtGui.QVBoxLayout()
            figLayoutLower = QtGui.QVBoxLayout()
            figLayoutUpper.addWidget(canvasUpper)
            figLayoutLower.addWidget(canvasLower)
            figWidgetUpper.setLayout(figLayoutUpper)
            figWidgetLower.setLayout(figLayoutLower)

            # Rescale the Widget and add the Widget to the layoutScrolls
            canvasUpper.setMinimumSize(canvasUpper.size())
            canvasLower.setMinimumSize(canvasLower.size())
            layoutScrollUpper.addWidget(figWidgetUpper)
            layoutScrollLower.addWidget(figWidgetLower)

        # Set the Upper and lower scroll content to the scroll content layout
        layoutScroll.addWidget(scrollContentUpper)
        layoutScroll.addWidget(scrollContentLower)

        # Set the Scroll Content to the scroll area and add the
        # Scroll Area to the layoutFrames
        self.scrollArea.setWidget(scrollContent)
        self.layoutFrames.addWidget(self.scrollArea)

        # Set the layoutScrollUpper to the scrollContentUpper and the
        # layoutFrames to the widgetFrames.
        scrollContentUpper.setLayout(layoutScrollUpper)
        scrollContentLower.setLayout(layoutScrollLower)
        scrollContent.setLayout(layoutScroll)
        self.widgetFrames.setLayout(self.layoutFrames)

    def _refresh_plot(self, figure, mode=None):
        """Refresh the plot"""

        # If there is already a plot, delete it
        if self.canvas is not None:
            self.canvas.setParent(None)
            self.toolbar.setParent(None)

        # Create Canvas that displays the figure
        self.canvas = FigureCanvasQTAgg(figure)
        self.canvas.setParent(self.widgetPlot)

        # Prepare for events
        self.canvas.mpl_connect(
            'button_press_event', self._unsaved_settings
            )
        self.mode = mode

        # Create Navigation widget for the canvas
        self.toolbar = NavigationToolbar2QTAgg(self.canvas, self.widgetPlot)
        # Set Layout
        self.layoutPlot.addWidget(self.toolbar)
        self.layoutPlot.addWidget(self.canvas)
        self.widgetPlot.setLayout(self.layoutPlot)

    def _change_frame(self, event):
        """Event when one clicks on one of the frames"""

        # Emit a signal with the current frame
        self.sigFrame.emit(self.dictFrames['{:s}'.format(event.canvas)])

    def _unsaved_settings(self, event):
        """Event when one clicks on the plot"""

        # Event specific variables
        eventMouseLeft = 1
        eventMouseRight = 3
        idxMouse = 0
        idxValue = 1
        idxMode = 2
        listEvent = [idxMouse, idxValue, idxMode]

        # If the click is inaxis and it was a left or right click,
        # emit the button and the xData
        if event.inaxes is not None and (
                event.button == eventMouseLeft or
                event.button == eventMouseRight
                ):
            listEvent[idxMouse] = event.button
            listEvent[idxValue] = event.xdata
            listEvent[idxMode] = self.mode
            self.sigRefresh.emit(listEvent)

    def _funny_plot(self):
        """Plot a list which shows that there is no Plot"""

        figFunny = pylab.figure(figsize=(3, 2), dpi=100)
        pylab.plot([0, 0, 1, 1], [0, 2, 0, 2], 'b')
        pylab.plot([2, 2, 3, 3, 2], [0, 2, 2, 0, 0], 'b')
        pylab.plot([4, 4, 5, 5, 4], [0, 2, 2, 1, 1], 'b')
        pylab.plot([6, 6, 7], [2, 0, 0], 'b')
        pylab.plot([8, 8, 9, 9, 8], [0, 2, 2, 0, 0], 'b')
        pylab.plot([10, 10.5, 10.5, 10.5, 11], [2, 2, 0, 2, 2], 'b')
        listCoordX = []
        listCoordY = []
        arrPhiFace = numpy.linspace(0, 2 * numpy.pi, 100)
        arrPhiMouth = numpy.linspace(0, numpy.pi, 100)
        for i in arrPhiFace:
            listCoordX.append(numpy.cos(i) + 13)
            listCoordY.append(numpy.sin(i) + 1)
        pylab.plot(listCoordX, listCoordY, 'b')
        listCoordX = []
        listCoordY = []
        for i in arrPhiMouth:
            listCoordX.append(0.5 * numpy.cos(i) + 13)
            listCoordY.append(0.5 * numpy.sin(i) + 0.25)
        pylab.plot(listCoordX, listCoordY, 'b')
        pylab.plot([12.75, 12.75], [1.25, 1.5], 'b')
        pylab.plot([13.25, 13.25], [1.25, 1.5], 'b')
        pylab.plot(
            [-0.25, -0.25, 14.25, 14.25, -0.25],
            [-0.25, 2.25, 2.25, -0.25, -0.25],
            'r'
            )

        pylab.xlim([-1, 15])
        pylab.ylim([-0.5, 2.5])
        pylab.grid()
        pylab.tight_layout()

        return figFunny

    def changeEvent(self, event):
        if event.type() == QtCore.QEvent.WindowStateChange:
            if self.windowState() & QtCore.Qt.WindowMinimized:
                self.setWindowState(event.oldState())
                return None

    def closeEvent(self, event):
        """Change the closeEvent to close the application cleanly"""

        # Ignore the incomming event when pressing the "X"
        event.ignore()
        self.sigClose.emit()


class SXDriftUnblur(QtGui.QMainWindow, Ui_MSMainWidget):

    def __init__(self, inputlist=None, inputfile=None, parent=None):
        super(SXDriftUnblur, self).__init__(parent)

        # Setup Gui Elements
        self.setupUi(self)

        # Initialize
        self._set_variables()
        self._connect_events()

        # Check parsed arguments
        if inputlist is not None:
            result = self._check_list_or_file(inputlist)
            if result == 'list':
                try:
                    listOfShiftFiles = numpy.genfromtxt(inputlist, dtype=None, unpack=True)
                except TypeError:
                    message = QtGui.QMessageBox(self)
                    message.setText('Empty File:\n{0}'.format(inputlist))
                    message.exec_()
                except ValueError:
                    message = QtGui.QMessageBox(self)
                    message.setText('File is not valid, only one column allowed:\n{0}'.format(inputlist))
                    message.exec_()
                else:
                    if len(numpy.shape(listOfShiftFiles)) > 1:
                        message = QtGui.QMessageBox(self)
                        message.setText('Too many columns. Expected one column:\n{0}'.format(inputlist))
                        message.exec_()
                    elif len(numpy.shape(listOfShiftFiles)) == 0:
                        self._fill_gui(inputlist=[str(listOfShiftFiles)], inputfile=inputfile)
                    else:
                        listOfShiftFiles = [os.path.relpath(name) for name in listOfShiftFiles]
                        self._fill_gui(inputlist=list(listOfShiftFiles), inputfile=inputfile)
            else:
                message = QtGui.QMessageBox(self)
                message.setText('Not a valid file name. Try again:\n{0}'.format(inputlist))
                message.exec_()
        elif inputfile is not None:
            result = self._check_list_or_file(inputfile)
            if result == 'file':
                self._fill_gui(inputfile=inputfile)
            elif result == 'list':
                try:
                    listOfShiftFiles = numpy.genfromtxt(inputfile, dtype=None, unpack=True)
                except TypeError:
                    message = QtGui.QMessageBox(self)
                    message.setText('Empty File:\n{0}'.format(inputfile))
                    message.exec_()
                except ValueError:
                    message = QtGui.QMessageBox(self)
                    message.setText('File is not valid, only one column allowed:\n{0}'.format(inputfile))
                    message.exec_()
                else:
                    if len(numpy.shape(listOfShiftFiles)) > 1:
                        message = QtGui.QMessageBox(self)
                        message.setText('Too many columns. Expected one column:\n{0}'.format(inputfile))
                        message.exec_()
                    elif len(numpy.shape(listOfShiftFiles)) == 0:
                        self._fill_gui(inputfile=str(listOfShiftFiles))
                    else:
                        self._fill_gui(inputlist=list(listOfShiftFiles))
            else:
                print('Error: {0} not found! Try again'.format(inputfile))

    def _check_list_or_file(self, data):
        """Check if the input name is a file or a list of files"""
        
        if '*' in data:
            return 'file'
        else:
            try:
                with open(data, 'r') as r:
                    first_line = r.readline()
                    if first_line.startswith('# Unblur'):
                        return 'file'
                    elif first_line.startswith('# full frame alignment'):
                        return 'file'
                    elif first_line.startswith('# Patch based alignment'):
                        return 'file'
                    else:
                        return 'list'
            except IOError:
                return None

    def _set_variables(self):
        """Set instance variables"""

        # Minimized
        self.minimized = False

        # Arrays
        self.arrData = None
        self.arrMicNumber = None

        # Variables
        self.strInputDir = None
        self.intNrFrames = None
        self.idxFirstFrame = None
        self.idxLastFrame = None
        self.varAnalyzeOne = None
        self.fileName = None

        # Indices show and hide
        self.idxVisible = 0
        self.idxRect = 1
        self.idxPos = 2

        # Indices Threshold and Widget
        self.idxStart = 0
        self.idxStop = 1
        self.idxStartSave = 2
        self.idxStopSave = 3
        self.idxCriterion = 4
        self.idxSaved = 5
        self.idxName = 6

        # output directory
        self.outputDir = None

        # Modes
        self.modeOverall = 'Overall'
        self.modeFrame = 'Frame'
        self.modeAngle = 'Angle'
        self.modeGeneral = 'General'
        self.modeAverage = 'Average'
        self.modeDrift = 'Drift'
        self.modeDriftPerFrame = 'Drift per Frame'
        self.modeAnglePerFrame = 'Angle per Frame'
        self.modePerMic = 'Drift per Micrograph'

        # DType
        self.dFile = 'fileName'
        self.dFileRaw = 'fileNameRaw'
        self.dMic = 'micName'
        self.dOverall = 'driftOverall'
        self.dFrame = 'driftFrame'
        self.dEnd = 'endToEndDrift'
        self.dMax = 'maxDistance'
        self.dMaxFirst = 'maxDistanceFirst'

        # Sort
        self.sortFile = 'File name'
        self.sortOverall = 'Overall drift'
        self.sortFrame = 'Drift per frame'
        self.sortEnd = 'End to end length'
        self.sortMax = 'Maximum distance between frames'
        self.sortMaxFirst = 'Maximum distance from start frame'

        # Plot Widgets
        self.msPlotDrift = SXUnblurPlot(title='Drift Plot')
        self.msPlotFrame = SXUnblurPlot(title='Drift per Frame Plot')
        self.msPlotAngle = SXUnblurPlot(title='Angle per Frame')
        self.msAllPlotFrameAvg = SXUnblurPlot(title='Average Drift per Frame')
        self.msAllPlotDrift = SXUnblurPlot(
            title='Overall Drift Histogram'
            )
        self.msAllPlotFrame = SXUnblurPlot(
            title='Drift per Frame Histogram', setframes=True
            )
        self.msAllPlotAngle = SXUnblurPlot(
            title='Angle per Frame Histogram', setframes=True
            )
        self.msAllPlotPerMic = SXUnblurPlot(title='Drift per Micrograph')

        # Dictionarys
        self.dictThresh = {}
        self.dictWidgets = {
            self.modeFrame: [
                self.leStartFrame,
                self.leStopFrame,
                self.leStartSaveFrame,
                self.leStopSaveFrame,
                self.chFrameCriterion
                ],
            self.modeAngle: [
                self.leStartAngle,
                self.leStopAngle,
                self.leStartSaveAngle,
                self.leStopSaveAngle,
                self.chAngleCriterion
                ],
            self.modeOverall: [
                self.leStartOverall,
                self.leStopOverall,
                self.leStartSaveOverall,
                self.leStopSaveOverall,
                self.chOverallCriterion
                ]
            }
        self.dictShow = {
            self.chPlotDriftMic: self.msPlotDrift,
            self.chPlotFrameMic: self.msPlotFrame,
            self.chPlotAngleMic: self.msPlotAngle,
            self.chAverageDriftPerFrame: self.msAllPlotFrameAvg,
            self.chPlotDrift: self.msAllPlotDrift,
            self.chPlotFrame: self.msAllPlotFrame,
            self.chPlotAngle: self.msAllPlotAngle,
            self.chPlotPerMic: self.msAllPlotPerMic
            }
        self.dictHide = {
            self.msPlotDrift: self.chPlotDriftMic,
            self.msPlotFrame: self.chPlotFrameMic,
            self.msPlotAngle: self.chPlotAngleMic,
            self.msAllPlotFrameAvg: self.chAverageDriftPerFrame,
            self.msAllPlotDrift: self.chPlotDrift,
            self.msAllPlotFrame: self.chPlotFrame,
            self.msAllPlotAngle: self.chPlotAngle,
            self.msAllPlotPerMic: self.chPlotPerMic
            }
        self.dictSort = {
            self.sortFile: self.dMic,
            self.sortOverall: self.dOverall,
            self.sortFrame: self.dFrame,
            self.sortEnd: self.dEnd,
            self.sortMax: self.dMax,
            self.sortMaxFirst: self.dMaxFirst
            }
        self.dictTranslate = {
            self.modeOverall: self.dOverall
            }
        self.dictButton = {
            self.modeOverall: self.pbSaveOverall,
            self.modeFrame: self.pbSaveFrame,
            self.modeAngle: self.pbSaveAngle,
            self.modeGeneral: self.pbSaveGeneral
            }
        self.dictColor = {
            'modified': "color: rgb(0, 150, 255);",
            'done': "color: rgb(0, 0, 0);"
            }
        self.dictVisible = {
            self.msPlotDrift: [
                self.msPlotDrift.isVisible(),
                self.msPlotDrift.rect(),
                self.msPlotDrift.pos()
                ],
            self.msPlotFrame: [
                self.msPlotFrame.isVisible(),
                self.msPlotFrame.rect(),
                self.msPlotFrame.pos()
                ],
            self.msPlotAngle: [
                self.msPlotAngle.isVisible(),
                self.msPlotAngle.rect(),
                self.msPlotAngle.pos()
                ],
            self.msAllPlotFrameAvg: [
                self.msAllPlotFrameAvg.isVisible(),
                self.msAllPlotFrameAvg.rect(),
                self.msAllPlotFrameAvg.pos()
                ],
            self.msAllPlotDrift: [
                self.msAllPlotDrift.isVisible(),
                self.msAllPlotDrift.rect(),
                self.msAllPlotDrift.pos()
                ],
            self.msAllPlotFrame: [
                self.msAllPlotFrame.isVisible(),
                self.msAllPlotFrame.rect(),
                self.msAllPlotFrame.pos()
                ],
            self.msAllPlotAngle: [
                self.msAllPlotAngle.isVisible(),
                self.msAllPlotAngle.rect(),
                self.msAllPlotAngle.pos()
                ],
            self.msAllPlotPerMic: [
                self.msAllPlotPerMic.isVisible(),
                self.msAllPlotPerMic.rect(),
                self.msAllPlotPerMic.pos()
                ]
            }

        # Lists
        self.listChecked = []
        self.listUnchecked = []
        self.listFile = []
        self.listCoordX = []
        self.listCoordY = []
        self.listFrames = []
        self.listAngles = []
        self.listDType = []
        self.connectedWidgets = []

        # List widget flags
        self.newItemFlags = \
            Qt.ItemFlags(Qt.ItemIsSelectable) | \
            Qt.ItemFlags(Qt.ItemIsEnabled) | \
            Qt.ItemFlags(Qt.ItemIsUserCheckable)

    def _disconnect_events(self):
        """disconnect the widgets to the events"""
        for entry in self.connectedWidgets:
            entry.disconnect()

    def _connect_events(self):
        """Connect the widgets to the events"""

        # Connect buttons
        self.pbSelectTxt.clicked.connect(self._select_txt)
        self.connectedWidgets.append(self.pbSelectTxt.clicked)
        self.pbImportList.clicked.connect(self._load_from_list)
        self.connectedWidgets.append(self.pbImportList.clicked)
        self.pbImportPattern.clicked.connect(self._load_from_pattern)
        self.connectedWidgets.append(self.pbImportPattern.clicked)
        #self.pbSelectAll.clicked.connect(self._select_all)
        #self.connectedWidgets.append(self.pbSelectAll.clicked)
        #self.pbInvert.clicked.connect(self._invert_selection)
        #self.connectedWidgets.append(self.pbInvert.clicked)
        self.pbAbout.clicked.connect(self._show_about)
        self.connectedWidgets.append(self.pbAbout.clicked)
        self.pbSaveSelected.clicked.connect(self._write_selection)
        self.connectedWidgets.append(self.pbSaveSelected.clicked)
        self.pbSaveOverall.clicked.connect(
            lambda: self._save_settings(mode=self.modeOverall)
            )
        self.connectedWidgets.append(self.pbSaveOverall.clicked)
        self.pbSaveGeneral.clicked.connect(
            lambda: self._save_settings(mode=self.modeGeneral)
            )
        self.connectedWidgets.append(self.pbSaveGeneral.clicked)
        self.pbSaveFrame.clicked.connect(
            lambda: self._save_settings(mode=self.modeFrame)
            )
        self.connectedWidgets.append(self.pbSaveFrame.clicked)
        self.pbSaveAngle.clicked.connect(
            lambda: self._save_settings(mode=self.modeAngle)
            )
        self.connectedWidgets.append(self.pbSaveAngle.clicked)
        self.pbApply.clicked.connect(self._apply_setting)
        self.connectedWidgets.append(self.pbApply.clicked)
        self.pbUncheckCriterion.clicked.connect(self._uncheck_angle_criterion)
        self.connectedWidgets.append(self.pbUncheckCriterion.clicked)
        self.pbSaveSettings.clicked.connect(self._write_settings)
        self.connectedWidgets.append(self.pbSaveSettings.clicked)
        self.pbLoadSettings.clicked.connect(self._load_settings)
        self.connectedWidgets.append(self.pbLoadSettings.clicked)

        # Connect list widget
        self.lsFiles.itemClicked.connect(self._current_info)
        self.connectedWidgets.append(self.lsFiles.itemClicked)

        # Connect entry widgets to calculation
        self.leFrameStart.returnPressed.connect(self._refresh_calculations)
        self.connectedWidgets.append(self.leFrameStart.returnPressed)
        self.leFrameStop.returnPressed.connect(self._refresh_calculations)
        self.connectedWidgets.append(self.leFrameStop.returnPressed)
        self.leStartOverall.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeOverall
                )
            )
        self.connectedWidgets.append(self.leStartOverall.editingFinished)
        self.leStopOverall.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeOverall
                )
            )
        self.connectedWidgets.append(self.leStopOverall.editingFinished)
        self.leStartGeneral.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeGeneral
                )
            )
        self.connectedWidgets.append(self.leStartGeneral.editingFinished)
        self.leStopGeneral.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeGeneral
                )
            )
        self.connectedWidgets.append(self.leStopGeneral.editingFinished)
        self.leStartFrame.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeFrame
                )
            )
        self.connectedWidgets.append(self.leStartFrame.editingFinished)
        self.leStopFrame.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeFrame
                )
            )
        self.connectedWidgets.append(self.leStopFrame.editingFinished)
        self.leStartAngle.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeAngle
                )
            )
        self.connectedWidgets.append(self.leStartAngle.editingFinished)
        self.leStopAngle.editingFinished.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeAngle
                )
            )
        self.connectedWidgets.append(self.leStopAngle.editingFinished)

        # Connect entry widgets to change color
        self.leFrameStart.textChanged.connect(
            lambda: self._change_color(
                widget=self.leFrameStart
                )
            )
        self.connectedWidgets.append(self.leFrameStart.textChanged)
        self.leFrameStop.textChanged.connect(
            lambda: self._change_color(
                widget=self.leFrameStop
                )
            )
        self.connectedWidgets.append(self.leFrameStop.textChanged)
        self.leStartOverall.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartOverall
                )
            )
        self.connectedWidgets.append(self.leStartOverall.textChanged)
        self.leStopOverall.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopOverall
                )
            )
        self.connectedWidgets.append(self.leStopOverall.textChanged)
        self.leStartGeneral.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartGeneral
                )
            )
        self.connectedWidgets.append(self.leStartGeneral.textChanged)
        self.leStopGeneral.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopGeneral
                )
            )
        self.connectedWidgets.append(self.leStopGeneral.textChanged)
        self.leStartFrame.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartFrame
                )
            )
        self.connectedWidgets.append(self.leStartFrame.textChanged)
        self.leStopFrame.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopFrame
                )
            )
        self.connectedWidgets.append(self.leStopFrame.textChanged)
        self.leStartAngle.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartAngle
                )
            )
        self.connectedWidgets.append(self.leStartAngle.textChanged)
        self.leStopAngle.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopAngle
                )
            )
        self.connectedWidgets.append(self.leStopAngle.textChanged)

        # Connect check boxes
        self.chDescending.setChecked(False)
        self.chDescending.stateChanged.connect(self._apply_sorting)
        self.connectedWidgets.append(self.chDescending.stateChanged)
        self.chSortSelected.setChecked(False)
        self.chSortSelected.stateChanged.connect(self._apply_sorting)
        self.connectedWidgets.append(self.chSortSelected.stateChanged)
        self.chOverallCriterion.setChecked(False)
        self.chOverallCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeOverall)
            )
        self.connectedWidgets.append(self.chOverallCriterion.stateChanged)
        self.chGeneralCriterion.setChecked(False)
        self.chGeneralCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeGeneral)
            )
        self.connectedWidgets.append(self.chGeneralCriterion.stateChanged)
        self.chFrameCriterion.setChecked(False)
        self.chFrameCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeFrame)
            )
        self.connectedWidgets.append(self.chFrameCriterion.stateChanged)
        self.chAngleCriterion.setChecked(False)
        self.chAngleCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeAngle)
            )
        self.connectedWidgets.append(self.chAngleCriterion.stateChanged)
        self.chPlotDriftMic.setChecked(False)
        self.chPlotDriftMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotDriftMic)
            )
        self.connectedWidgets.append(self.chPlotDriftMic.stateChanged)
        self.chPlotFrameMic.setChecked(False)
        self.chPlotFrameMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotFrameMic)
            )
        self.connectedWidgets.append(self.chPlotFrameMic.stateChanged)
        self.chPlotAngleMic.setChecked(False)
        self.chPlotAngleMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotAngleMic)
            )
        self.connectedWidgets.append(self.chPlotAngleMic.stateChanged)
        self.chAverageDriftPerFrame.setChecked(False)
        self.chAverageDriftPerFrame.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chAverageDriftPerFrame)
            )
        self.connectedWidgets.append(self.chAverageDriftPerFrame.stateChanged)
        self.chPlotDrift.setChecked(False)
        self.chPlotDrift.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotDrift)
            )
        self.connectedWidgets.append(self.chPlotDrift.stateChanged)
        self.chPlotFrame.setChecked(False)
        self.chPlotFrame.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotFrame)
            )
        self.connectedWidgets.append(self.chPlotFrame.stateChanged)
        self.chPlotAngle.setChecked(False)
        self.chPlotAngle.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotAngle)
            )
        self.connectedWidgets.append(self.chPlotAngle.stateChanged)
        self.chPlotPerMic.setChecked(False)
        self.chPlotPerMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotPerMic)
            )
        self.connectedWidgets.append(self.chPlotPerMic.stateChanged)

        # Connect combo boxes
        self.cbSort.activated.connect(self._apply_sorting)
        self.connectedWidgets.append(self.cbSort.activated)
        self.cbFrame.activated.connect(
            lambda: self._plot_threshold(mode=self.modeFrame, fill=True)
            )
        self.connectedWidgets.append(self.cbFrame.activated)
        self.cbAngle.activated.connect(
            lambda: self._plot_threshold(mode=self.modeAngle, fill=True)
            )
        self.connectedWidgets.append(self.cbAngle.activated)

        # Connect canvas close events
        self.msPlotDrift.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msPlotDrift)
            )
        self.connectedWidgets.append(self.msPlotDrift.sigClose)
        self.msPlotFrame.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msPlotFrame)
            )
        self.connectedWidgets.append(self.msPlotFrame.sigClose)
        self.msPlotAngle.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msPlotAngle)
            )
        self.connectedWidgets.append(self.msPlotAngle.sigClose)
        self.msAllPlotFrameAvg.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotFrameAvg)
            )
        self.connectedWidgets.append(self.msAllPlotFrameAvg.sigClose)
        self.msAllPlotDrift.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotDrift)
            )
        self.connectedWidgets.append(self.msAllPlotDrift.sigClose)
        self.msAllPlotFrame.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotFrame)
            )
        self.connectedWidgets.append(self.msAllPlotFrame.sigClose)
        self.msAllPlotAngle.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotAngle)
            )
        self.connectedWidgets.append(self.msAllPlotAngle.sigClose)
        self.msAllPlotPerMic.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotPerMic)
            )
        self.connectedWidgets.append(self.msAllPlotPerMic.sigClose)

        # Connect canvas refresh frame event
        self.msAllPlotFrame.sigFrame.connect(self._refresh_frame)
        self.connectedWidgets.append(self.msAllPlotFrame.sigFrame)
        self.msAllPlotAngle.sigFrame.connect(self._refresh_frame)
        self.connectedWidgets.append(self.msAllPlotAngle.sigFrame)

        # Connect canvas refresh plot event
        self.msAllPlotDrift.sigRefresh.connect(
            self._refresh_histogram_mouse
            )
        self.connectedWidgets.append(self.msAllPlotDrift.sigRefresh)
        self.msAllPlotFrame.sigRefresh.connect(
            self._refresh_histogram_mouse
            )
        self.connectedWidgets.append(self.msAllPlotFrame.sigRefresh)
        self.msAllPlotAngle.sigRefresh.connect(
            self._refresh_histogram_mouse
            )
        self.connectedWidgets.append(self.msAllPlotAngle.sigRefresh)

    def _enable_all(self):
        """Enable all widgets"""

        self.pbSelectTxt.setEnabled(True)
        #self.pbSelectAll.setEnabled(True)
        #self.pbInvert.setEnabled(True)
        self.pbAbout.setEnabled(True)
        self.pbSaveSelected.setEnabled(True)
        self.pbSaveGeneral.setEnabled(True)
        self.pbApply.setEnabled(True)
        self.pbSaveSettings.setEnabled(True)
        self.lsFiles.setEnabled(True)
        self.leOutputName.setEnabled(True)
        self.leFrameStart.setEnabled(True)
        self.leFrameStop.setEnabled(True)
        self.leStartOverall.setEnabled(True)
        self.leStopOverall.setEnabled(True)
        self.leStartGeneral.setEnabled(True)
        self.leStopGeneral.setEnabled(True)
        self.leStartFrame.setEnabled(True)
        self.leStopFrame.setEnabled(True)
        self.leStartAngle.setEnabled(True)
        self.leStopAngle.setEnabled(True)
        self.chDescending.setEnabled(True)
        self.chSortSelected.setEnabled(True)
        self.chOverallCriterion.setEnabled(True)
        self.chGeneralCriterion.setEnabled(True)
        self.chFrameCriterion.setEnabled(True)
        self.chAngleCriterion.setEnabled(True)
        self.chPlotDriftMic.setEnabled(True)
        self.chPlotFrameMic.setEnabled(True)
        self.chPlotAngleMic.setEnabled(True)
        self.chAverageDriftPerFrame.setEnabled(True)
        self.chPlotDrift.setEnabled(True)
        self.chPlotFrame.setEnabled(True)
        self.chPlotAngle.setEnabled(True)
        self.chPlotPerMic.setEnabled(True)
        self.cbSort.setEnabled(True)
        self.cbFrame.setEnabled(True)
        self.cbAngle.setEnabled(True)

    def _default_color(self):
        """Set all colors to black"""

        self.leFrameStart.setStyleSheet(self.dictColor['done'])
        self.leFrameStop.setStyleSheet(self.dictColor['done'])
        self.leStartOverall.setStyleSheet(self.dictColor['done'])
        self.leStopOverall.setStyleSheet(self.dictColor['done'])
        self.leStartGeneral.setStyleSheet(self.dictColor['done'])
        self.leStopGeneral.setStyleSheet(self.dictColor['done'])
        self.leStartFrame.setStyleSheet(self.dictColor['done'])
        self.leStopFrame.setStyleSheet(self.dictColor['done'])
        self.leStartAngle.setStyleSheet(self.dictColor['done'])
        self.leStopAngle.setStyleSheet(self.dictColor['done'])

    def _change_color(self, widget):
        """Changes the color of the text to red"""

        widget.setStyleSheet(self.dictColor['modified'])

    def _select_txt(self):
        """Open the drift files and do the drift calculations"""

        # Find Directory
        strInputFile = str(QtGui.QFileDialog.getOpenFileName(
            directory=os.getcwd(),
            options=QtGui.QFileDialog.DontUseNativeDialog,
            filter='Unblur (*.txt);;MotionCor2 (*.log);;All (*)'
            ))

        # If the return value is not empty, fill the line edit
        if self.strInputDir != '':
            self.lePattern.clear()
            self.lePattern.insert(strInputFile)

    def _apply_sorting(self):
        """Sort the list widget entrys"""

        # Check if the sort selected is checked
        if self.chSortSelected.isChecked():

            # Create an empty array for entrys
            arrChecked = numpy.empty(
                len(self.listChecked), dtype=self.listDType
                )
            arrUnchecked = numpy.empty(
                len(self.listUnchecked), dtype=self.listDType
                )

            # Fill the arrays with the entrys from the lists
            for number, entry in enumerate(self.listChecked):
                arrChecked[number] = self.arrData[
                    self.arrData[self.dFile] == entry
                    ]
            for number, entry in enumerate(self.listUnchecked):
                arrUnchecked[number] = self.arrData[
                    self.arrData[self.dFile] == entry
                    ]

            # Sort the arrays in ascending or descending order
            arrSortedChecked = numpy.sort(
                arrChecked, order=self.dictSort[str(self.cbSort.currentText())]
                )
            if self.chDescending.isChecked():
                arrSortedChecked = arrSortedChecked[::-1]

            arrSortedUnchecked = numpy.sort(
                arrUnchecked, order=self.dictSort[str(self.cbSort.currentText())]
                )
            if self.chDescending.isChecked():
                arrSortedUnchecked = arrSortedUnchecked[::-1]

            # Fill the list widget
            self.lsFiles.clear()
            for file in arrSortedUnchecked[self.dFile]:
                newItem = QtGui.QListWidgetItem(file)
                newItem.setFlags(self.newItemFlags)
                newItem.setCheckState(Qt.Unchecked)
                self.lsFiles.addItem(newItem)
            for file in arrSortedChecked[self.dFile]:
                newItem = QtGui.QListWidgetItem(file)
                newItem.setFlags(self.newItemFlags)
                newItem.setCheckState(Qt.Checked)
                self.lsFiles.addItem(newItem)

        # If sort selected is not checked
        else:

            # Create a set for faster search and sort the array
            # in ascending or descending order.
            setChecked = set(self.listChecked)
            arrSorted = numpy.sort(
                self.arrData, order=self.dictSort[str(self.cbSort.currentText())]
                )
            if self.chDescending.isChecked():
                arrSorted = arrSorted[::-1]

            # Fill the list widget, but leave the selection untouched
            self.lsFiles.clear()
            for file in arrSorted[self.dFile]:
                newItem = QtGui.QListWidgetItem(file)
                newItem.setFlags(self.newItemFlags)
                if file in setChecked:
                    newItem.setCheckState(Qt.Checked)
                else:
                    newItem.setCheckState(Qt.Unchecked)
                self.lsFiles.addItem(newItem)

    def _fill_gui(self, inputlist=None, inputfile=None):
        """Do the Calculations and fill the GUI"""

        # Reset variables
        self._disconnect_events()
        self._set_variables()
        self._connect_events()

        # Reset lists
        self.listChecked = []
        self.listUnchecked = []
        self.listFile = []
        self.listCoordX = []
        self.listCoordY = []
        self.listFrames = []
        self.listAngles = []
        self.listDType = []

        # Message Box!
        messageBox = QtGui.QMessageBox()
        messageBox.setText('Do drift calculations...')
        messageBox.setStandardButtons(QtGui.QMessageBox().NoButton)
        messageBox.show()
        print('Do drift calculations...')

        # Clear the list widget and load the shift files as list
        if inputfile is not None:
            strDirectory = os.getcwd()
            strSuffix = os.path.relpath(inputfile)
            path = '{:s}/{:s}'.format(strDirectory, strSuffix)
            self.outputDir = os.path.relpath(path[:-len(path.split('/')[-1])])

            self.listFile = glob.glob(os.path.relpath(path))

            # Try to find micrograph numbers
            if len(strSuffix.split('/')) == 1:
                listSplit = strSuffix.split('*')
            else:
                listSplit = strSuffix.split('/')[-1].split('*')
            if len(listSplit) != 2:
                self.arrMicNumber = numpy.arange(len(self.listFile))
                print(
                    'Warning: Could not identify micrograph serial number.\n' +
                    'X-axis of the plots do not represent' +
                    'the serial micrograph number.'
                    )
            else:
                listMicNumber = []
                varBreak = False
                for file in self.listFile:
                    name = file.split('/')[-1]
                    prefix = name.split(listSplit[-1])
                    try:
                        number = int(prefix[0].split(listSplit[0])[-1])
                    except:
                        self.arrMicNumber = numpy.arange(len(self.listFile))
                        varBreak = True
                        print(
                            'Warning: Could not identify micrograph serial number.\n' +
                            'X-axis of the plots do not represent' +
                            'the serial micrograph number.'
                            )
                        break
                    else:
                        listMicNumber.append(number)

                if not varBreak:
                    self.arrMicNumber = numpy.array(listMicNumber)

        if inputlist is not None:
            if inputfile is not None:
                self.listFile = [name for name in inputlist if os.path.relpath(name) in self.listFile]
            else:
                self.listFile = inputlist
            self.arrMicNumber = numpy.arange(len(self.listFile))

        # If no or few files were found
        if not self.listFile:
            messageBox2 = QtGui.QMessageBox()
            messageBox2.setText(
                'Error: No matching drift files found or pattern' + 
                ' and list entries does not match\n{:s}'.format(inputfile)
                )
            messageBox2.exec_()
            messageBox.hide()
            return None
        elif len(self.listFile) <= 5:
            self.varAnalyzeOne = True
            print(
                '\nWarning: !!!! Only few shift files selected, ' +
                'so plots of all micrographs could not work as expected. !!!!\n'
                )

        # Load and calculate data for the first time
        value = self._first_time_calculations()

        if len(self.arrMicNumber) != len(self.listFile):
            self.arrMicNumber = numpy.arange(len(self.listFile))
            print(
                'Warning: Files were corrupt, lost track of micrograph serial number.\n' +
                'X-axis of the plots do not represent' +
                'the serial micrograph number.'
                )

        if value is not None:
            # Fill list widget
            self.lsFiles.clear()
            for file in self.listFile:
                if os.path.exists(file) and file in self.arrData[self.dFileRaw]:
                    file = file.split('/')[-1]
                    newItem = QtGui.QListWidgetItem(file)
                    newItem.setFlags(self.newItemFlags)
                    newItem.setCheckState(Qt.Checked)
                    self.lsFiles.addItem(newItem)

            # Fill the frame widgets
            self.leFrameStart.clear()
            self.leFrameStop.clear()
            self.leFrameStart.insert('1')
            self.leFrameStop.insert('{:d}'.format(self.intFrames))

            # Do Calculations
            self._calculations(oldfirst=1, oldlast=self.intFrames)

            # Fill the GUI
            self._refresh_gui()
            
            # Enable all widges
            self._enable_all()

            # Go use the GUI
            print('Done')
        messageBox.hide()

    def _refresh_gui(self):
        """Refresh the GUI"""

        # Fill the Dictionarys
        self._fill_dictionary()

        # Select the first item
        self.lsFiles.setItemSelected(self.lsFiles.item(0), True)
        self._current_info(item=self.lsFiles.item(0))

        # Fill the General widgets with maximum and minimum values
        arrFrames = numpy.array([
            item for entry in self.arrData[self.listFrames] for item in entry
            ])
        self.leStartGeneral.setText(
            '{:f}'.format(numpy.amin(arrFrames))
            )
        self.leStopGeneral.setText(
            '{:f}'.format(numpy.amax(arrFrames))
            )

        self.varOldStartGeneral = numpy.amin(arrFrames)
        self.varOldStopGeneral = numpy.amax(arrFrames)

        # plot the frame dependant plots
        self._plot_scroll(mode=self.modeFrame)
        self._plot_scroll(mode=self.modeAngle)
        self._plot_single(mode=self.modeAverage)
        self._plot_single(mode=self.modePerMic)
        self._plot_threshold(mode=self.modeOverall)
        self._plot_threshold(mode=self.modeFrame)
        self._plot_threshold(mode=self.modeAngle)

        # Fill the threshold widgets
        self._fill_widgets(mode=self.modeFrame)
        self._fill_widgets(mode=self.modeAngle)
        self._fill_widgets(mode=self.modeOverall)

        # Enable Save general button
        self.pbSaveGeneral.setEnabled(True)

        # Save selected micrographs to list
        self._save_selection()

        # Print it black
        self._default_color()

        # Check Overall and general criterion
        self.chOverallCriterion.setChecked(True)
        self.chGeneralCriterion.setChecked(True)

    def _first_time_calculations(self):
        """Calculate the Drift"""

        # Check how many frames are there
        self.intFrames = 0
        input_typ = None
        for files in self.listFile:
            try:
                with open(files, 'r') as f:
                    first_line = f.readline()
                    if first_line.startswith('# Unblur'):
                        input_typ = 'Unblur'
                        for linenumber, line in enumerate(f):
                            if linenumber == 2:
                                self.intFrames = int(line.split()[-1])
                                break
                    elif first_line.startswith('# full frame alignment'):
                        input_typ = 'MotionCor2'
                        for linenumber, line in enumerate(f):
                            pass
                        self.intFrames = linenumber + 1
                    elif first_line.startswith('# Patch based alignment'):
                        input_typ = 'MotionCor2'
                        for linenumber, line in enumerate(f):
                            pass
                        self.intFrames = linenumber - 1
                    else:
                        raise IOError
                break
            except IOError:
                continue
        else:
            message = QtGui.QMessageBox(self)
            message.setText(
                    'No files in given file list available:\n{0}'.format(
                        self.fileName
                        )
                    )
            message.exec_()
            return None

        # Fill the listDType with the right number of frames and angles.
        # Add entrys to the lists and fill the translation Dictionary.
        self.listDType = [
            (self.dFile, '|S200'),
            (self.dFileRaw, '|S200'),
            (self.dMic, '|S200'),
            (self.dOverall, '<f8'),
            (self.dFrame, '<f8'),
            (self.dEnd, '<f8'),
            (self.dMax, '<f8'),
            (self.dMaxFirst, '<f8')
            ]
        for index in xrange(1, self.intFrames + 1):
            self.listDType.append(('x{:d}'.format(index), '<f8'))
            self.listDType.append(('y{:d}'.format(index), '<f8'))
            self.listCoordX.append('x{:d}'.format(index))
            self.listCoordY.append('y{:d}'.format(index))

            if index <= self.intFrames - 1:
                self.listDType.append(('frame{:d}'.format(index), '<f8'))
                self.listFrames.append('frame{:d}'.format(index))
                self.dictTranslate.update({
                    'Frame {:d}'.format(index): 'frame{:d}'.format(index)
                    })

            if index <= self.intFrames - 2:
                self.listDType.append(('angle{:d}'.format(index), '<f8'))
                self.listAngles.append('angle{:d}'.format(index))
                self.dictTranslate.update({
                    'Angle {:d}'.format(index): 'angle{:d}'.format(index)
                    })

        validFiles = []
        for file in self.listFile:

            # Import the data
            try:
                if input_typ == 'Unblur':
                    arrCoord = numpy.genfromtxt(file, unpack=True)
                elif input_typ == 'MotionCor2':
                    arrCoord = numpy.genfromtxt(file, unpack=True)[1:]
                    # Transpose the array
                    arrCoord = numpy.transpose(arrCoord)
            except ValueError:
                print('Warning: File corrupt, skip:', file)
                continue

            if len(arrCoord) != self.intFrames:
                print('Warning: File does not have {0} Frames, skip:'.format(self.intFrames), file)
                continue

            # Get the micrograph name
            try:
                with open(file, 'r') as f:
                    pass
            except IOError:
                continue
            validFiles.append(file)

        self.listFile = validFiles

        # Create an empty array
        self.arrData = numpy.empty(len(self.listFile), dtype=self.listDType)

        # Load the data
        idxX = 0
        idxY = 1
        for number, file in enumerate(self.listFile):

            # Import the data
            try:
                if input_typ == 'Unblur':
                    arrCoord = numpy.genfromtxt(file, unpack=True)
                elif input_typ == 'MotionCor2':
                    arrCoord = numpy.genfromtxt(file, unpack=True)[1:]
                    # Transpose the array
                    arrCoord = numpy.transpose(arrCoord)
            except ValueError:
                print('Warning: File corrupt, skip:', file)
                continue

            if len(arrCoord) != self.intFrames:
                print('Warning: File does not have {0} Frames, skip:'.format(self.intFrames), file)
                continue

            # Get the micrograph name
            try:
                if input_typ == 'Unblur':
                    with open(file, 'r') as f:
                        self.arrData[self.dMic][number] = f.readline().split()[-1]
                elif input_typ == 'MotionCor2':
                    self.arrData[self.dMic][number] = file.split('/')[-1].replace('0-Patch-Full.log', '').replace('0-Full.log', '')
            except IOError:
                continue

            # Get the file name
            self.arrData[self.dFile][number] = file.split('/')[-1]
            self.arrData[self.dFileRaw][number] = file

            for linenumber, coord in enumerate(arrCoord):
                self.arrData['x{:d}'.format(linenumber + 1)][number] = \
                    round(coord[idxX], 6)
                self.arrData['y{:d}'.format(linenumber + 1)][number] = \
                    round(coord[idxY], 6)

            # Calculate per frame drift
            for index in xrange(1, self.intFrames):

                fltDistanceX = \
                    self.arrData['x{:d}'.format(index)][number] - \
                    self.arrData['x{:d}'.format(index + 1)][number]

                fltDistanceY = \
                    self.arrData['y{:d}'.format(index)][number] - \
                    self.arrData['y{:d}'.format(index + 1)][number]

                fltDistance = numpy.sqrt(fltDistanceX**2 + fltDistanceY**2)

                self.arrData['frame{:d}'.format(index)][number] = \
                    round(fltDistance, 6)
        return 'Done!'

    def _calculations(self, oldfirst, oldlast):
        """Calculate End-to-End drift, Overall drift, Maximum Distance
        Maximum Distance from first frame, Drift per Frame
        """

        # Get the first and last frame variables
        try:
            self.idxFirstFrame = int(self.leFrameStart.text())
            self.idxLastFrame = int(self.leFrameStop.text())
        except ValueError:
            messageBox = QtGui.QMessageBox(self)
            messageBox.setText('Frame must be integer!')
            messageBox.exec_()
            self.leFrameStart.setText('{:d}'.format(oldfirst))
            self.leFrameStop.setText('{:d}'.format(oldlast))
            self.leFrameStart.setStyleSheet(self.dictColor['done'])
            self.leFrameStop.setStyleSheet(self.dictColor['done'])
            self.idxFirstFrame = oldfirst
            self.idxLastFrame = oldlast
            return False

        # Check some abort situations
        if self.idxFirstFrame < 1:
            # Warning box when refreshing frames
            warningBox = QtGui.QMessageBox(self)
            warningBox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
            warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
            warningBox.setText(
                'Start frame too small (minimum 1)!\n' +
                'Continue with minimum value?'
                )
            warningBox.exec_()
            if warningBox.result() == QtGui.QMessageBox.No:
                self.leFrameStart.setText('{:d}'.format(oldfirst))
                self.leFrameStart.setStyleSheet(self.dictColor['done'])
                self.idxFirstFrame = oldfirst
                return False
            elif warningBox.result() == QtGui.QMessageBox.Yes:
                self.leFrameStart.setText('{:d}'.format(1))
                self.leFrameStart.setStyleSheet(self.dictColor['done'])
                self.idxFirstFrame = 1

        if self.idxLastFrame > self.intFrames:
            warningBox = QtGui.QMessageBox(self)
            warningBox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
            warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
            warningBox.setText(
                'Stop frame too large (maximum {:d})!\n'.format(self.intFrames) +
                'Continue with maximum value?'
                )
            warningBox.exec_()
            if warningBox.result() == QtGui.QMessageBox.No:
                self.leFrameStart.setText('{:d}'.format(oldlast))
                self.leFrameStart.setStyleSheet(self.dictColor['done'])
                self.idxLastFrame = oldlast
                return False
            elif warningBox.result() == QtGui.QMessageBox.Yes:
                self.leFrameStop.setText('{:d}'.format(self.intFrames))
                self.leFrameStop.setStyleSheet(self.dictColor['done'])
                self.idxLastFrame = self.intFrames

        if self.idxLastFrame <= self.idxFirstFrame:
            messageBox = QtGui.QMessageBox(self)
            messageBox.setText(
                'Start frame must be smaller stop frame!'
                )
            messageBox.exec_()
            self.leFrameStart.setText('{:d}'.format(oldfirst))
            self.leFrameStop.setText('{:d}'.format(oldlast))
            self.leFrameStart.setStyleSheet(self.dictColor['done'])
            self.leFrameStop.setStyleSheet(self.dictColor['done'])
            self.idxFirstFrame = oldfirst
            self.idxLastFrame = oldlast
            return False

        # Do the calculations for each array entry
        for entry in self.arrData:

            # Set Variables
            varMaximum = 0
            varMaximumFirst = 0
            varOverallDrift = 0

            # Use all wanted frames
            for index in xrange(self.idxFirstFrame, self.idxLastFrame):

                # Calculate Angles
                if index <= self.idxLastFrame - 2:
                    fltDistanceX1 = \
                        entry['x{:d}'.format(index)] - \
                        entry['x{:d}'.format(index + 1)]

                    fltDistanceY1 = \
                        entry['y{:d}'.format(index)] - \
                        entry['y{:d}'.format(index + 1)]

                    fltDistanceX2 = \
                        entry['x{:d}'.format(index + 1)] - \
                        entry['x{:d}'.format(index + 2)]

                    fltDistanceY2 = \
                        entry['y{:d}'.format(index + 1)] - \
                        entry['y{:d}'.format(index + 2)]

                    fltPointProduct = \
                        fltDistanceX1 * fltDistanceX2 + \
                        fltDistanceY1 * fltDistanceY2

                    fltAbsOne = numpy.sqrt(
                        fltDistanceX1**2 +
                        fltDistanceY1**2
                        )

                    fltAbsTwo = numpy.sqrt(
                        fltDistanceX2**2 +
                        fltDistanceY2**2
                        )

                    fltAbsProduct = fltAbsOne * fltAbsTwo

                    if 0 <= fltAbsProduct < 0.000001:
                        angleCos = 1
                    else:
                        angleCos = fltPointProduct / fltAbsProduct
                    if angleCos > 1:
                        angleCos = 1
                    elif angleCos < -1:
                        angleCos = -1
                    angleRad = numpy.arccos(angleCos)

                    angleDeg = angleRad * 180 / numpy.pi

                    # Save calculations of angles
                    entry['angle{:d}'.format(index)] = round(angleDeg, 6)

                # Calculate Overall drift and maximum drift
                if index <= self.idxLastFrame - 1:

                    # Maximum Drift between two frames
                    if varMaximum < entry['frame{:d}'.format(index)]:
                        varMaximum = entry['frame{:d}'.format(index)]

                    varOverallDrift += entry['frame{:d}'.format(index)]

                # Calculate end to end distance and maximum distance from
                # first frame
                fltDistanceX = \
                    entry['x{:d}'.format(self.idxFirstFrame)] - \
                    entry['x{:d}'.format(index)]

                fltDistanceY = \
                    entry['y{:d}'.format(self.idxFirstFrame)] - \
                    entry['y{:d}'.format(index)]

                fltDistance = numpy.sqrt(fltDistanceX**2 + fltDistanceY**2)

                if fltDistance > varMaximumFirst:
                    varMaximumFirst = fltDistance

            # Save calculations
            entry[self.dMaxFirst] = round(varMaximumFirst, 6)
            entry[self.dEnd] = round(fltDistance, 6)
            entry[self.dOverall] = round(varOverallDrift, 6)
            entry[self.dMax] = round(varMaximum, 6)
            entry[self.dFrame] = round(
                varOverallDrift / (self.idxLastFrame - self.idxFirstFrame), 6
                )

        return True

    def _fill_dictionary(self):
        """Fill the Dictionarys"""

        # Clear the combo box
        self.cbAngle.clear()
        self.cbFrame.clear()

        # Fill the dictionary with overall information
        self.dictThresh.update({
            self.modeOverall: [
                self.idxStart,
                self.idxStop,
                self.idxStartSave,
                self.idxStopSave,
                self.idxCriterion,
                self.idxSaved,
                None
                ]
            })

        # Fill the entry
        self.dictThresh[self.modeOverall][
            self.idxStart
            ] = numpy.min(self.arrData[self.dOverall])
        self.dictThresh[self.modeOverall][
            self.idxStop
            ] = numpy.max(self.arrData[self.dOverall])
        self.dictThresh[self.modeOverall][
            self.idxStartSave
            ] = numpy.min(self.arrData[self.dOverall])
        self.dictThresh[self.modeOverall][
            self.idxStopSave
            ] = numpy.max(self.arrData[self.dOverall])
        self.dictThresh[self.modeOverall][
            self.idxCriterion
            ] = False
        self.dictThresh[self.modeOverall][
            self.idxSaved
            ] = True

        # Fill the dictionary with frame and angle
        for index in xrange(self.idxFirstFrame, self.idxLastFrame):

            # With angles
            if index <= self.idxLastFrame - 1 \
                    and index <= self.intFrames - 2:
                # Update the Entry
                self.dictThresh.update({
                    'Angle {:d}'.format(index): [
                        self.idxStart,
                        self.idxStop,
                        self.idxStartSave,
                        self.idxStopSave,
                        self.idxCriterion,
                        self.idxSaved,
                        self.idxName
                        ]
                    })

                # Fill the entry
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxName
                    ] = 'angle{:d}'.format(index)
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxStart
                    ] = numpy.min(self.arrData['angle{:d}'.format(index)])
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxStop
                    ] = numpy.max(self.arrData['angle{:d}'.format(index)])
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxStartSave
                    ] = numpy.min(self.arrData['angle{:d}'.format(index)])
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxStopSave
                    ] = numpy.max(self.arrData['angle{:d}'.format(index)])
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxCriterion
                    ] = False
                self.dictThresh['Angle {:d}'.format(index)][
                    self.idxSaved
                    ] = True

                # Add to combo box

                self.cbAngle.addItem('Angle {:d}'.format(index))

            # With frames, Update the entry
            self.dictThresh.update({
                'Frame {:d}'.format(index): [
                    self.idxStart,
                    self.idxStop,
                    self.idxStartSave,
                    self.idxStopSave,
                    self.idxCriterion,
                    self.idxSaved,
                    self.idxName
                    ]
                })

            # Fill the entry
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxName
                ] = 'frame{:d}'.format(index)
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxStart
                ] = numpy.min(self.arrData['frame{:d}'.format(index)])
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxStop
                ] = numpy.max(self.arrData['frame{:d}'.format(index)])
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxStartSave
                ] = numpy.min(self.arrData['frame{:d}'.format(index)])
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxStopSave
                ] = numpy.max(self.arrData['frame{:d}'.format(index)])
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxCriterion
                ] = False
            self.dictThresh['Frame {:d}'.format(index)][
                self.idxSaved
                ] = True

            # Add to combo box
            self.cbFrame.addItem('Frame {:d}'.format(index))

    def _refresh_calculations(self, goon=False):
        ''''Refresh the calculations for the chosen frames'''

        # Skip warning when loading from saved data
        # Safe old state of first and last frame
        intOldFirst = self.idxFirstFrame
        intOldLast = self.idxLastFrame

        # Warning box when refreshing frames
        warningBox = QtGui.QMessageBox(self)
        warningBox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
        warningBox.setText(
            'Threshold settings will be lost when calculating new drift data!\n' +
            'Do you really want to continue?'
            )
        if not goon:
            warningBox.exec_()

        if warningBox.result() == QtGui.QMessageBox.Yes or goon:

            # Message Box!
            messageBox = QtGui.QMessageBox()
            messageBox.setText('Do drift calculations...')
            messageBox.setStandardButtons(QtGui.QMessageBox().NoButton)
            messageBox.show()
            print('Do drift calculations...')

            # Set all criterions to False
            for key in self.dictThresh:
                self.dictThresh[key][self.idxCriterion] = False

            # If the input is correct continue
            varContinue = self._calculations(oldfirst=intOldFirst, oldlast=intOldLast)
            if varContinue:
                # Refresh GUI
                self._refresh_gui()

            # Refresh overall calculations
            self._invert_selection()
            self._invert_selection()
            # Use the GUI!
            messageBox.hide()
            print('Done')

        else:
            self.leFrameStart.setText('{:d}'.format(intOldFirst))
            self.leFrameStop.setText('{:d}'.format(intOldLast))
            self.leFrameStart.setStyleSheet(self.dictColor['done'])
            self.leFrameStop.setStyleSheet(self.dictColor['done'])

    def _fill_widgets(self, mode):
        """Fill the Widget with the current Setting"""

        # Get the right entry
        if mode == self.modeFrame or mode == self.modeGeneral:
            mode = self.modeFrame
            entry = str(self.cbFrame.currentText())
        elif mode == self.modeAngle:
            if str(self.cbAngle.currentText()) != '':
                entry = str(self.cbAngle.currentText())
            else:
                # If there is no Angle entry
                self.dictWidgets[mode][self.idxStart].setText(
                    '{:s}'.format('No Entry')
                    )
                self.dictWidgets[mode][self.idxStop].setText(
                    '{:s}'.format('No Entry')
                    )
                self.dictWidgets[mode][self.idxStartSave].setText(
                    '{:s}'.format('No Entry')
                    )
                self.dictWidgets[mode][self.idxStopSave].setText(
                    '{:s}'.format('No Entry')
                    )
                self.dictWidgets[mode][self.idxCriterion].setChecked(
                    False
                    )
                return None
        elif mode == self.modeOverall:
            entry = self.modeOverall

        # Fill the widgets
        self.dictWidgets[mode][self.idxStart].setText(
            '{:f}'.format(self.dictThresh[entry][self.idxStart])
            )
        self.dictWidgets[mode][self.idxStop].setText(
            '{:f}'.format(self.dictThresh[entry][self.idxStop])
            )
        self.dictWidgets[mode][self.idxStartSave].setText(
            '{:f}'.format(self.dictThresh[entry][self.idxStartSave])
            )
        self.dictWidgets[mode][self.idxStopSave].setText(
            '{:f}'.format(self.dictThresh[entry][self.idxStopSave])
            )
        self.dictWidgets[mode][self.idxCriterion].setChecked(
            self.dictThresh[entry][self.idxCriterion]
            )

        # Enable or disable the save button, depending on the
        # save index.
        if self.dictThresh[entry][self.idxSaved]:
            self.dictButton[mode].setEnabled(False)
        else:
            self.dictButton[mode].setEnabled(True)

        # Black is the color
        self.dictWidgets[mode][self.idxStart].setStyleSheet(self.dictColor['done'])
        self.dictWidgets[mode][self.idxStop].setStyleSheet(self.dictColor['done'])
        self.dictWidgets[mode][self.idxStartSave].setStyleSheet(self.dictColor['done'])
        self.dictWidgets[mode][self.idxStopSave].setStyleSheet(self.dictColor['done'])

    def _plot_scroll(self, mode):
        """Plot the plots for the Scroll Widget"""

        # List for figures
        listFigure = []

        # Special case for angle: Just the last frame is selected
        if self.idxFirstFrame == self.idxLastFrame - 1 and \
                self.idxLastFrame >= self.intFrames and \
                mode == self.modeAngle:

            # Create a funny plot
            figFunny = pylab.figure(figsize=(3, 2), dpi=100)
            pylab.plot([0, 0, 1, 1], [0, 2, 0, 2], 'b')
            pylab.plot([2, 2, 3, 3, 2], [0, 2, 2, 0, 0], 'b')
            pylab.plot([4, 4, 5, 5, 4], [0, 2, 2, 1, 1], 'b')
            pylab.plot([6, 6, 7], [2, 0, 0], 'b')
            pylab.plot([8, 8, 9, 9, 8], [0, 2, 2, 0, 0], 'b')
            pylab.plot([10, 10.5, 10.5, 10.5, 11], [2, 2, 0, 2, 2], 'b')
            listCoordX = []
            listCoordY = []
            arrPhiFace = numpy.linspace(0, 2 * numpy.pi, 100)
            arrPhiMouth = numpy.linspace(0, numpy.pi, 100)
            for i in arrPhiFace:
                listCoordX.append(numpy.cos(i) + 13)
                listCoordY.append(numpy.sin(i) + 1)
            pylab.plot(listCoordX, listCoordY, 'b')
            listCoordX = []
            listCoordY = []
            for i in arrPhiMouth:
                listCoordX.append(0.5 * numpy.cos(i) + 13)
                listCoordY.append(0.5 * numpy.sin(i) + 0.25)
            pylab.plot(listCoordX, listCoordY, 'b')
            pylab.plot([12.75, 12.75], [1.25, 1.5], 'b')
            pylab.plot([13.25, 13.25], [1.25, 1.5], 'b')
            pylab.plot(
                [-0.25, -0.25, 14.25, 14.25, -0.25],
                [-0.25, 2.25, 2.25, -0.25, -0.25],
                'r'
                )

            pylab.xlim([-1, 15])
            pylab.ylim([-0.5, 2.5])
            pylab.grid()
            pylab.tight_layout()
            self.msAllPlotAngle._fill_widget_frames_single(
                frame=figFunny, name=mode
                )

            pylab.close(figFunny)

            # Make the angle Threshold widgets unenabled
            self.chAngleCriterion.setChecked(False)
            self.chAngleCriterion.setEnabled(False)
            self.leStartAngle.setEnabled(False)
            self.leStopAngle.setEnabled(False)
            self.cbAngle.setEnabled(False)
            self.pbSaveAngle.setEnabled(False)
            self.pbUncheckCriterion.setEnabled(False)

        # Special case does not occure
        else:

            # Make the angle Threshold widgets unenabled
            self.chAngleCriterion.setEnabled(True)
            self.leStartAngle.setEnabled(True)
            self.leStopAngle.setEnabled(True)
            self.cbAngle.setEnabled(True)

            for number in xrange(self.idxFirstFrame, self.idxLastFrame):

                # Special case for angle: break if the number is too huge
                if mode == self.modeAngle and number > self.intFrames - 2:
                    break

                # Find name and plot histogram
                strName = self.dictThresh['{:s} {:d}'.format(mode, number)][
                    self.idxName
                    ]
                figFrames = pylab.figure(figsize=(3, 2), dpi=100)
                if self.varAnalyzeOne:
                    intBins = 3
                else:
                    intBins = int(self.lsFiles.count() / 3)
                arrBins = numpy.linspace(
                    numpy.min(self.arrData[strName]),
                    numpy.max(self.arrData[strName]) + 0.0001,
                    intBins
                    )
                hist = numpy.histogram(
                    self.arrData[strName],
                    bins=arrBins
                    )
                pylab.plot(hist[1][:-1], hist[0], 'k.')
                pylab.grid()
                pylab.xlim([hist[1][0] - hist[1][-1] * 0.1, hist[1][-1] * 1.1])
                pylab.ylim([0, numpy.max(hist[0]) + 1])
                if mode == self.modeFrame:
                    pylab.xlabel(r'Drift / Angstrom')
                if mode == self.modeAngle:
                    pylab.xlabel(r'Angle / Degree')
                pylab.ylabel(r'Nr. of Micrographs')
                pylab.title('{:s} {:d}'.format(mode, number))
                pylab.tight_layout()
                listFigure.append(figFrames)
                pylab.close(figFrames)

            # Special case: just one figure
            if len(listFigure) == 1:
                if mode == self.modeFrame:
                    self.msAllPlotFrame._fill_widget_frames_single(
                        frame=listFigure[0], name=mode
                        )
                elif mode == self.modeAngle:
                    self.msAllPlotAngle._fill_widget_frames_single(
                        frame=listFigure[0], name=mode
                        )
            # More than one figure
            else:
                if mode == self.modeFrame:
                    self.msAllPlotFrame._fill_widget_frames(
                        framesaslist=listFigure, name=mode,
                        framestart=self.idxFirstFrame,
                        framestop=self.idxLastFrame
                        )
                elif mode == self.modeAngle:
                    self.msAllPlotAngle._fill_widget_frames(
                        framesaslist=listFigure, name=mode,
                        framestart=self.idxFirstFrame,
                        framestop=self.idxLastFrame - 1
                        )

    def _plot_threshold(self, mode, fill=False):
        """Plot the right values for selecting a threshold"""

        # Figure and bin size for all plots
        figThresh = pylab.figure(figsize=(5, 4), dpi=100)
        if self.varAnalyzeOne:
            intBins = 3
        else:
            intBins = int(self.lsFiles.count() / 3)

        # Special case, if there is no angle available
        if self.idxFirstFrame == self.idxLastFrame - 1 and \
                self.idxLastFrame >= self.intFrames and \
                mode == self.modeAngle:

            # Plot a list which shows that there is no Plot
            figFrames = pylab.figure(figsize=(3, 2), dpi=100)
            pylab.plot([0, 0, 1, 1], [0, 2, 0, 2], 'b')
            pylab.plot([2, 2, 3, 3, 2], [0, 2, 2, 0, 0], 'b')
            pylab.plot([4, 4, 5, 5, 4], [0, 2, 2, 1, 1], 'b')
            pylab.plot([6, 6, 7], [2, 0, 0], 'b')
            pylab.plot([8, 8, 9, 9, 8], [0, 2, 2, 0, 0], 'b')
            pylab.plot([10, 10.5, 10.5, 10.5, 11], [2, 2, 0, 2, 2], 'b')
            listCoordX = []
            listCoordY = []
            arrPhiFace = numpy.linspace(0, 2 * numpy.pi, 100)
            arrPhiMouth = numpy.linspace(0, numpy.pi, 100)
            for i in arrPhiFace:
                listCoordX.append(numpy.cos(i) + 13)
                listCoordY.append(numpy.sin(i) + 1)
            pylab.plot(listCoordX, listCoordY, 'b')
            listCoordX = []
            listCoordY = []
            for i in arrPhiMouth:
                listCoordX.append(0.5 * numpy.cos(i) + 13)
                listCoordY.append(0.5 * numpy.sin(i) + 0.25)
            pylab.plot(listCoordX, listCoordY, 'b')
            pylab.plot([12.75, 12.75], [1.25, 1.5], 'b')
            pylab.plot([13.25, 13.25], [1.25, 1.5], 'b')
            pylab.plot(
                [-0.25, -0.25, 14.25, 14.25, -0.25],
                [-0.25, 2.25, 2.25, -0.25, -0.25],
                'r'
                )

            pylab.xlim([-1, 15])
            pylab.ylim([-0.5, 2.5])
            pylab.grid()
            pylab.tight_layout()
            self.msAllPlotAngle._refresh_plot(
                figure=figFrames
                )
            pylab.close(figFrames)

        # Normal case
        else:

            # Set settings depending on mode
            if mode == self.modeOverall:
                arrBins = numpy.linspace(
                    numpy.min(self.arrData[self.dOverall]),
                    numpy.max(self.arrData[self.dOverall]) + 0.0001,
                    intBins
                    )
                arrInput = self.arrData[self.dOverall]
                fltStart = self.dictThresh[mode][self.idxStart]
                fltStop = self.dictThresh[mode][self.idxStop]
                fltStartSave = self.dictThresh[mode][self.idxStartSave]
                fltStopSave = self.dictThresh[mode][self.idxStopSave]
                strTitle = r'Overall Drift Histogram'
                strXLabel = r'Drift / Angstrom'
                varOutput = self.msAllPlotDrift

            if mode == self.modeFrame:
                strCurrent = str(self.cbFrame.currentText())
                # Abort if there is no current entry
                if strCurrent == '':
                    return None
                idxCurrent = self.dictTranslate[strCurrent]
                arrBins = numpy.linspace(
                    numpy.min(self.arrData[idxCurrent]),
                    numpy.max(self.arrData[idxCurrent]) + 0.0001,
                    intBins
                    )
                arrInput = self.arrData[idxCurrent]
                fltStart = self.dictThresh[strCurrent][self.idxStart]
                fltStop = self.dictThresh[strCurrent][self.idxStop]
                fltStartSave = self.dictThresh[strCurrent][self.idxStartSave]
                fltStopSave = self.dictThresh[strCurrent][self.idxStopSave]
                strTitle = r'{:s} Histogram'.format(strCurrent)
                strXLabel = r'Drift / Angstrom'
                varOutput = self.msAllPlotFrame

            if mode == self.modeAngle:
                strCurrent = str(self.cbAngle.currentText())
                # Abort if there is no current entry
                if strCurrent == '':
                    return None
                idxCurrent = self.dictTranslate[strCurrent]
                arrBins = numpy.linspace(
                    numpy.min(self.arrData[idxCurrent]),
                    numpy.max(self.arrData[idxCurrent]) + 0.0001,
                    intBins
                    )
                arrInput = self.arrData[idxCurrent]
                fltStart = self.dictThresh[strCurrent][self.idxStart]
                fltStop = self.dictThresh[strCurrent][self.idxStop]
                fltStartSave = self.dictThresh[strCurrent][self.idxStartSave]
                fltStopSave = self.dictThresh[strCurrent][self.idxStopSave]
                strTitle = r'{:s} Histogram'.format(strCurrent)
                strXLabel = r'Angle / Degree'
                varOutput = self.msAllPlotAngle

            # Calculate and plot the histogram
            hist = numpy.histogram(
                arrInput,
                bins=arrBins
                )
            pylab.plot(hist[1][:-1], hist[0], 'k.')
            pylab.vlines(
                fltStart, 0, numpy.max(hist[0]) + 1,
                color='b', linestyle='dashed'
                )
            pylab.vlines(
                fltStop, 0, numpy.max(hist[0]) + 1,
                color='r', linestyle='dashed'
                )
            pylab.vlines(
                fltStartSave, 0, numpy.max(hist[0]) + 1,
                color='b', linestyle='solid'
                )
            pylab.vlines(
                fltStopSave, 0, numpy.max(hist[0]) + 1,
                color='r', linestyle='solid'
                )
            pylab.vlines(
                hist[1][0], 0, 3,
                color='k', linestyle='solid'
                )
            pylab.vlines(
                hist[1][-1], 0, 3,
                color='k', linestyle='solid'
                )
            pylab.grid()
            pylab.xlim([hist[1][0] - hist[1][-1] * 0.1, hist[1][-1] * 1.1])
            pylab.ylim([0, numpy.max(hist[0]) + 1])
            pylab.xlabel(strXLabel)
            pylab.ylabel(r'Nr. of Micrographs')
            pylab.title(strTitle)
            pylab.tight_layout()
            varOutput._refresh_plot(figure=figThresh, mode=mode)
            pylab.close(figThresh)

        # # Only do this when the combo box changes
        if fill:
            # Refresh the current widget
            self._fill_widgets(mode=mode)

    def _save_settings(self, mode):
        """Save the settings"""

        # Overall mode
        if mode == self.modeOverall:

            # Equalize saved and not saved and fill the widgets
            self.dictThresh[mode][self.idxStartSave] = \
                self.dictThresh[mode][self.idxStart]
            self.dictThresh[mode][self.idxStopSave] = \
                self.dictThresh[mode][self.idxStop]
            self.dictThresh[mode][self.idxSaved] = True

        # Frame mode
        if mode == self.modeFrame:

            # Equalize saved and not saved and fill the widgets
            strCurrent = str(self.cbFrame.currentText())
            self.dictThresh[strCurrent][self.idxStartSave] = \
                self.dictThresh[strCurrent][self.idxStart]
            self.dictThresh[strCurrent][self.idxStopSave] = \
                self.dictThresh[strCurrent][self.idxStop]
            self.dictThresh[strCurrent][self.idxSaved] = True

            # Enable the general save button
            self.dictButton[self.modeGeneral].setEnabled(True)

        # Angle mode
        if mode == self.modeAngle:

            # Equalize saved and not saved and fill the widgets
            strCurrent = str(self.cbAngle.currentText())
            self.dictThresh[strCurrent][self.idxStartSave] = \
                self.dictThresh[strCurrent][self.idxStart]
            self.dictThresh[strCurrent][self.idxStopSave] = \
                self.dictThresh[strCurrent][self.idxStop]
            self.dictThresh[strCurrent][self.idxSaved] = True

        # General mode
        if mode == self.modeGeneral:

            # Try to cast as float, if it fails show an error
            try:
                fltStart = float(self.leStartGeneral.text())
                fltStop = float(self.leStopGeneral.text())
            except ValueError:
                messageBox = QtGui.QMessageBox()
                messageBox.setText('General input needs to be a float!')
                messageBox.exec_()

            # If it not fails equalize saved and not saved and fill the widgets
            else:

                # Set all frames to the general settings
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxStartSave
                        ] = fltStart
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxStopSave
                        ] = fltStop
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxStart
                        ] = fltStart
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxStop
                        ] = fltStop

                # Disable the general and frame save button
                # and set the mode to frames.
                self.dictButton[mode].setEnabled(False)
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    self.dictThresh['Frame {:d}'.format(index)][self.idxSaved] = True
                mode = self.modeFrame

        # Refresh Plot and widgets, disable the save button
        self.dictButton[mode].setEnabled(False)
        self._fill_widgets(mode=mode)
        self._plot_threshold(mode=mode)

    def _mark_as_criterion(self, mode):
        """Mark the current frame or option as criterion"""

        # Overall mode
        if mode == self.modeOverall:
            if self.chOverallCriterion.isChecked():
                self.dictThresh[mode][self.idxCriterion] = True
            else:
                self.dictThresh[mode][self.idxCriterion] = False

        # Frame mode
        if mode == self.modeFrame:
            strCurrent = str(self.cbFrame.currentText())

            # If the current check box is now checked, set it true
            if self.chFrameCriterion.isChecked():
                self.dictThresh[strCurrent][self.idxCriterion] = True

                # If all check boxes are checked,
                # mark the general check box as checked.
                # Else partially checked.
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    if not self.dictThresh['Frame {:d}'.format(index)][
                            self.idxCriterion
                            ]:
                        self.chGeneralCriterion.setCheckState(Qt.PartiallyChecked)
                        break
                    elif index == self.idxLastFrame - 1:
                        self.chGeneralCriterion.setCheckState(Qt.Checked)

            # If the current check box is now unchecked, set it false
            else:
                self.dictThresh[strCurrent][self.idxCriterion] = False

                # If all check boxes are unchecked,
                # mark the general check box as unchecked.
                # Else partially checked.
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    if self.dictThresh['Frame {:d}'.format(index)][
                            self.idxCriterion
                            ]:
                        self.chGeneralCriterion.setCheckState(Qt.PartiallyChecked)
                        break
                    elif index == self.idxLastFrame - 1:
                        self.chGeneralCriterion.setCheckState(Qt.Unchecked)

        # Angle mode
        if mode == self.modeAngle:
            strCurrent = str(self.cbAngle.currentText())

            # If the current check box is now checked, set it True
            if self.chAngleCriterion.isChecked():
                self.dictThresh[strCurrent][self.idxCriterion] = True
                # Enable the uncheck all button
                self.pbUncheckCriterion.setEnabled(True)

            # If the current check box is now unchecked, set it False
            else:
                self.dictThresh[strCurrent][self.idxCriterion] = False

                # If all check boxes are unchecked, disable the button
                for index in xrange(
                        self.idxFirstFrame,
                        self.idxFirstFrame + self.cbAngle.count()
                        ):
                    if self.dictThresh['Angle {:d}'.format(index)][self.idxCriterion]:
                        break
                    if index == self.idxFirstFrame + self.cbAngle.count() - 1:
                        self.pbUncheckCriterion.setEnabled(False)

        # General mode
        if mode == self.modeGeneral:

            # If its checked, check all frames
            if self.chGeneralCriterion.checkState() == Qt.Checked:
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxCriterion
                        ] = True

                # Set the current frame check box as checked
                self.chFrameCriterion.setChecked(True)

            # If its unchecked, uncheck all frames
            elif self.chGeneralCriterion.checkState() == Qt.Unchecked:
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxCriterion
                        ] = False

                # Set the current frame check box unchecked.
                # Disable the Tristate mode of the general check box.
                self.chFrameCriterion.setChecked(False)
                self.chGeneralCriterion.setTristate(False)

        # # Refresh widgets
        # self._fill_widgets(mode=mode)

    def _uncheck_angle_criterion(self):
        """Uncheck all angle criterions"""

        # Set all criterions to False
        for index in xrange(self.idxFirstFrame, self.idxLastFrame):
            if index >= self.intFrames - 1:
                break
            self.dictThresh['Angle {:d}'.format(index)][
                self.idxCriterion
                ] = False

        self.chAngleCriterion.setChecked(False)
        self.pbUncheckCriterion.setEnabled(False)

    def _apply_setting(self):
        """Apply the saved settings"""

        # Ask the user if he really wants to continue
        warningBox = QtGui.QMessageBox(self)
        warningBox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
        warningBox.setText(
                'Do you really want to apply criteria?\n' +
                'All selections will be lost.'
            )
        warningBox.exec_()

        if warningBox.result() == QtGui.QMessageBox.No:
            return None

        # Copy of the data
        arrBetweenThres = numpy.copy(self.arrData)

        # Check all dictionary entrys and shrink it depending on threshold
        for key in self.dictThresh:
            if self.dictThresh[key][self.idxCriterion]:
                arrBetweenThres = arrBetweenThres[
                    arrBetweenThres[self.dictTranslate[key]] >=
                    self.dictThresh[key][self.idxStartSave]
                    ]
                arrBetweenThres = arrBetweenThres[
                    arrBetweenThres[self.dictTranslate[key]] <=
                    self.dictThresh[key][self.idxStopSave]
                    ]

        # Set the file names for faster search and select this files
        setFileNames = set(arrBetweenThres[self.dFile])
        for index in xrange(int(self.lsFiles.count())):
            if str(self.lsFiles.item(index).text()) in setFileNames:
                self.lsFiles.item(index).setCheckState(Qt.Checked)
            else:
                self.lsFiles.item(index).setCheckState(Qt.Unchecked)

        # Save check selection
        self._save_selection()

    def _refresh_histogram_mouse(self, eventaslist):
        """Refresh the Histograms in the Widgets with mouse event"""

        # Indices
        eventMouseLeft = 1
        eventMouseRight = 3
        idxMouse = 0
        idxValue = 1
        idxMode = 2

        # Set the varInput to the right index
        if eventaslist[idxMode] == self.modeAngle:
            varInput = str(self.cbAngle.currentText())
        if eventaslist[idxMode] == self.modeFrame:
            varInput = str(self.cbFrame.currentText())
        if eventaslist[idxMode] == self.modeOverall:
            varInput = self.modeOverall

        # Enable the related save button
        self.dictThresh[varInput][self.idxSaved] = False

        # Choose the right widget, depending on the mouse event
        if eventaslist[idxMouse] == eventMouseLeft:
            fltStart = eventaslist[idxValue]
            fltStop = self.dictThresh[varInput][
                self.idxStop
                ]
        if eventaslist[idxMouse] == eventMouseRight:
            fltStart = self.dictThresh[varInput][
                self.idxStart
                ]
            fltStop = eventaslist[idxValue]

        # Start value must be smaller then stop value
        if fltStart > fltStop:

            # If event left, both values are the start value
            if eventaslist[idxMouse] == eventMouseLeft:
                self.dictThresh[varInput][
                    self.idxStop
                    ] = fltStart
                self.dictThresh[varInput][
                    self.idxStart
                    ] = fltStart

            # If event left, both values are the stop value
            if eventaslist[idxMouse] == eventMouseRight:
                self.dictThresh[varInput][
                    self.idxStop
                    ] = fltStop
                self.dictThresh[varInput][
                    self.idxStart
                    ] = fltStop

        # If everything is alright
        else:
            self.dictThresh[varInput][
                self.idxStart
                ] = fltStart
            self.dictThresh[varInput][
                self.idxStop
                ] = fltStop

        # Fill the widgets and refresh plot
        self._fill_widgets(mode=eventaslist[idxMode])
        self._plot_threshold(mode=eventaslist[idxMode])

    def _refresh_histogram_widget(self, event, mode):
        """Refresh the Histograms in the Widgets with widget event"""

        # Set the settings for different modes
        if mode == self.modeAngle:
            varStart = self.leStartAngle
            varStop = self.leStopAngle
            varInput = str(self.cbAngle.currentText())
            varStartOld = self.dictThresh[varInput][self.idxStart]
            varStopOld = self.dictThresh[varInput][self.idxStop]
        if mode == self.modeFrame:
            varStart = self.leStartFrame
            varStop = self.leStopFrame
            varInput = str(self.cbFrame.currentText())
            varStartOld = self.dictThresh[varInput][self.idxStart]
            varStopOld = self.dictThresh[varInput][self.idxStop]
        if mode == self.modeOverall:
            varStart = self.leStartOverall
            varStop = self.leStopOverall
            varInput = self.modeOverall
            varStartOld = self.dictThresh[varInput][self.idxStart]
            varStopOld = self.dictThresh[varInput][self.idxStop]
        if mode == self.modeGeneral:
            varStart = self.leStartGeneral
            varStop = self.leStopGeneral
            varStartOld = self.varOldStartGeneral
            varStopOld = self.varOldStopGeneral

        # Try to cast this as float
        try:
            fltStart = float(varStart.text())
            fltStop = float(varStop.text())
        except ValueError:
            messageBox = QtGui.QMessageBox()
            messageBox.setText(
                'Error with {:s}! Input must be float!'.format(
                    mode
                    )
                )
            messageBox.exec_()

            varStart.setText('{:f}'.format(varStartOld))
            varStop.setText('{:f}'.format(varStopOld))
            varStart.setStyleSheet(self.dictColor['done'])
            varStop.setStyleSheet(self.dictColor['done'])

        # Everything is alright
        else:

            # Settings for general mode
            if mode == self.modeGeneral:

                # Break with wrong input
                if fltStart > fltStop:
                    messageBox = QtGui.QMessageBox()
                    messageBox.setText(
                        'Error with {:s}! Larger must be smaller then Smaller! ;)'.format(
                            mode
                            )
                        )
                    messageBox.exec_()
                    varStart.setText('{:f}'.format(varStartOld))
                    varStop.setText('{:f}'.format(varStopOld))
                    varStart.setStyleSheet(self.dictColor['done'])
                    varStop.setStyleSheet(self.dictColor['done'])
                    return None

                self.varOldStartGeneral = fltStart
                self.varOldStopGeneral = fltStop

                # Set the settings to all frames
                for index in xrange(self.idxFirstFrame, self.idxLastFrame - 1):
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxStart
                        ] = fltStart
                    self.dictThresh['Frame {:d}'.format(index)][
                        self.idxStop
                        ] = fltStop

                # Enable the save button for general and all frames
                self.dictButton[mode].setEnabled(True)
                for index in xrange(self.idxFirstFrame, self.idxLastFrame):
                    self.dictThresh['Frame {:d}'.format(index)][self.idxSaved] = False

                # Set the Text to the general widgets and turn it black
                varStart.setText('{:f}'.format(fltStart))
                varStop.setText('{:f}'.format(fltStop))
                varStart.setStyleSheet(self.dictColor['done'])
                varStop.setStyleSheet(self.dictColor['done'])

                # Set the mode to Frame
                mode = self.modeFrame

            # Setting for modes except General
            else:

                # Check, if start is smaller then stop
                if fltStart > fltStop:
                    if event == 'start':
                        # Set both to start
                        self.dictThresh[varInput][
                            self.idxStop
                            ] = fltStart
                        self.dictThresh[varInput][
                            self.idxStart
                            ] = fltStart
                    if event == 'stop':
                        # Set both to stop
                        self.dictThresh[varInput][
                            self.idxStop
                            ] = fltStop
                        self.dictThresh[varInput][
                            self.idxStart
                            ] = fltStop
                # If everything is finde
                else:
                    # Set the things where they belong
                    self.dictThresh[varInput][
                        self.idxStart
                        ] = fltStart
                    self.dictThresh[varInput][
                        self.idxStop
                        ] = fltStop

                # Set it unsaved
                self.dictThresh[varInput][self.idxSaved] = False

            # Refresh widgets and plots
            self._fill_widgets(mode=mode)
            self._plot_threshold(mode=mode)

    def _refresh_frame(self, event):
        """Refresh, which Plot will be shown"""

        # Refresh the current plot to the chosen frame

        # Frame
        if str(event[:5]) == self.modeFrame:
            intIndex = self.cbFrame.findText(event)
            self.cbFrame.setCurrentIndex(intIndex)
            self._fill_widgets(mode=self.modeFrame)
            self._plot_threshold(mode=self.modeFrame)

        # Angle
        if str(event[:5]) == self.modeAngle:
            intIndex = self.cbAngle.findText(event)
            self.cbAngle.setCurrentIndex(intIndex)
            self._fill_widgets(mode=self.modeAngle)
            self._plot_threshold(mode=self.modeAngle)

    def _hide_plot(self, msplot):
        """Hide the Plot Widget when its closed"""

        # Uncheck the related check box and hide the widget
        msplot.hide()
        self.dictHide[msplot].setCheckState(Qt.Unchecked)
        self.dictVisible.update({
            msplot: [
                msplot.isVisible(),
                msplot.rect(),
                msplot.pos()
                ]
            })

    def _show_plot(self, checkbox):
        """Show the Plot Widget"""
        msplot = self.dictShow[checkbox]
        # If the checkbox is checked, show the related widget
        if checkbox.checkState() == Qt.Checked:
            msplot.setGeometry(self.dictVisible[msplot][self.idxRect])
            msplot.move(self.dictVisible[msplot][self.idxPos])
            msplot.activateWindow()
            msplot.show()
            self.dictVisible.update({
                msplot: [
                    msplot.isVisible(),
                    msplot.rect(),
                    msplot.pos()
                    ]
                })
        # Otherwise hide it
        elif checkbox.checkState() == Qt.Unchecked:
            msplot.hide()
            self.dictVisible.update({
            msplot: [
                    msplot.isVisible(),
                    msplot.rect(),
                    msplot.pos()
                    ]
                })

#    def _select_all(self):
#        """Select all entrys"""
#
#        # Set all items to checked and save the current selection state
#        for index in xrange(self.lsFiles.count()):
#            self.lsFiles.item(index).setCheckState(Qt.Checked)
#
#        # Save new selection in cache
#        self._save_selection()

    def _invert_selection(self):
        """Invert Selection"""

        # Invert the selection and save the current selection state
        for index in xrange(int(self.lsFiles.count())):
            if self.lsFiles.item(index).checkState() == Qt.Checked:
                self.lsFiles.item(index).setCheckState(Qt.Unchecked)
            else:
                self.lsFiles.item(index).setCheckState(Qt.Checked)

        # Save new selection in cache
        self._save_selection()

    def _current_info(self, item):
        """Refresh the Plots and Entrys of the current Info"""

        # get the current entry
        arrCurrentEntry = self.arrData[
            self.arrData[self.dFile] == str(item.text())
            ]

        # Clear all widgets
        self.leCurrentMicName.clear()
        self.leCurrentOverallDrift.clear()
        self.leCurrentFrameDrift.clear()
        self.leCurrentEndToEndDrift.clear()
        self.leCurrentMaxDistance.clear()
        self.leCurrentMaxDistanceZero.clear()

        # Fill all entrys
        self.leCurrentMicName.insert(
            '{:s}'.format(arrCurrentEntry[self.dMic][0].split('/')[-1])
            )
        self.leCurrentOverallDrift.insert(
            '{:f}'.format(arrCurrentEntry[self.dOverall][0])
            )
        self.leCurrentFrameDrift.insert(
            '{:f}'.format(arrCurrentEntry[self.dFrame][0])
            )
        self.leCurrentEndToEndDrift.insert(
            '{:f}'.format(arrCurrentEntry[self.dEnd][0])
            )
        self.leCurrentMaxDistance.insert(
            '{:f}'.format(arrCurrentEntry[self.dMax][0])
            )
        self.leCurrentMaxDistanceZero.insert(
            '{:f}'.format(arrCurrentEntry[self.dMaxFirst][0])
            )

        # Plot info
        self._plot_single(mode=self.modeDrift, item=item)
        self._plot_single(mode=self.modeDriftPerFrame, item=item)
        self._plot_single(mode=self.modeAnglePerFrame, item=item)

        # Refresh the checked info
        self._save_selection()

    def _plot_single(self, mode, item=None):
        """Plots without Threshold"""

        # Create a figure to plot on
        figSingle = pylab.figure()

        # If the mode is Average
        if mode == self.modeAverage:

            # Define plot settings
            arrX = numpy.linspace(1, self.intFrames - 1, self.intFrames - 1)
            arrY = numpy.zeros(self.intFrames - 1)
            for number in arrX:
                arrY[int(number - 1)] = \
                    numpy.sum(self.arrData['frame{:d}'.format(int(number))]) / \
                    int(self.lsFiles.count())
            strTitle = r'Average drift per Frame'
            strXLabel = r'Frame'
            strYLabel = r'Average Drift / Angstrom'
            varOutput = self.msAllPlotFrameAvg
            pylab.xlim([arrX[0] - arrX[-1] * 0.1, arrX[-1] * 1.1])

        elif mode == self.modePerMic:

            # Sort array by number
            listDType = []
            listDType.append(('micDrift', '<f8'))
            listDType.append(('micNum', '<i8'))
            arrXY = numpy.empty(len(self.arrMicNumber), dtype=listDType)
            arrXY['micNum'] = self.arrMicNumber
            arrXY['micDrift'] = self.arrData[self.dOverall]
            arrXY = numpy.sort(arrXY, order='micNum')

            # Define plot settings
            arrX = arrXY['micNum']
            arrY = arrXY['micDrift']
            strTitle = r'Drift per Micrograph'
            strXLabel = r'Micrograph'
            strYLabel = r'Overall Drift / Angstrom'
            varOutput = self.msAllPlotPerMic
            pylab.xlim([arrX[0] - arrX[-1] * 0.01, arrX[-1] * 1.01])

        # If the mode is not Average
        else:

            # Get the array entry of the current selected list item
            arrCurrentEntry = self.arrData[
                self.arrData[self.dFile] == str(item.text())
                ]

            # Plot settings depending on the mode
            if mode == self.modeDrift:
                arrX = list(arrCurrentEntry[self.listCoordX][0])
                arrY = list(arrCurrentEntry[self.listCoordY][0])
                strTitle = r'Drift of file {:s}'.format(arrCurrentEntry[self.dFile][0])
                strXLabel = r'Drift X / Angstrom'
                strYLabel = r'Drift Y / Angstrom'
                varOutput = self.msPlotDrift
                pylab.plot(
                    arrCurrentEntry[self.listCoordX[0]],
                    arrCurrentEntry[self.listCoordY[0]],
                    'dg', label=r'Frame 1'
                    )
                pylab.plot(
                    arrCurrentEntry[self.listCoordX[-1]],
                    arrCurrentEntry[self.listCoordY[-1]],
                    'dr', label=r'Frame {:d}'.format(self.intFrames)
                    )
                pylab.legend(loc='best')

            if mode == self.modeDriftPerFrame:
                arrX = numpy.linspace(1, self.intFrames - 1, self.intFrames - 1)
                arrY = list(arrCurrentEntry[self.listFrames][0])
                strTitle = r'Frame to frame drift of file {:s}'.format(
                    arrCurrentEntry[self.dFile][0]
                    )
                strXLabel = r'Frame'
                strYLabel = r'Drift / Angstrom'
                varOutput = self.msPlotFrame

            if mode == self.modeAnglePerFrame:
                arrX = numpy.linspace(1, self.intFrames - 2, self.intFrames - 2)
                arrY = list(arrCurrentEntry[self.listAngles][0])
                strTitle = 'Frame to frame angle of file {:s}'.format(
                    arrCurrentEntry[self.dFile][0]
                    )
                strXLabel = r'Frame'
                strYLabel = r'Angle / Degree'
                varOutput = self.msPlotAngle

        # Plot it!
        if not mode == self.modePerMic:
            pylab.plot(
                arrX,
                arrY,
                'kx'
                )
        pylab.plot(
            arrX,
            arrY,
            'k-'
            )
        pylab.grid()
        pylab.title(strTitle)
        pylab.xlabel(strXLabel)
        pylab.ylabel(strYLabel)
        varOutput._refresh_plot(figure=figSingle)
        pylab.close(figSingle)

    def _save_selection(self):
        """Save the selected micrographs"""

        # Fill the lists with current entrys
        listChecked = []
        listUnchecked = []
        for index in xrange(int(self.lsFiles.count())):
            if self.lsFiles.item(index).checkState() == Qt.Checked:
                listChecked.append(str(self.lsFiles.item(index).text()))
            else:
                listUnchecked.append(str(self.lsFiles.item(index).text()))

        # If the selection has changed, calculate the selection info
        if set(listChecked) != set(self.listChecked):
            self.listChecked = listChecked
            self.listUnchecked = listUnchecked

            # Uncheck sorted button
            self.chSortSelected.setChecked(False)

            # Create an array for the checked entrys and add all checked data to it
            arrChecked = numpy.empty(
                len(self.listChecked), dtype=self.listDType
                )
            for number, entry in enumerate(self.listChecked):
                arrChecked[number] = self.arrData[
                    self.arrData[self.dFile] == entry
                    ]

            # If no item is selected, set everything to 0
            if len(self.listChecked) == 0:
                fltOverallDrift = 0
                fltFrameDrift = 0
                fltEndToEndDrift = 0
                fltMaxDistance = 0
                fltMaxDistanceZero = 0

            # Else calculate the entrys
            else:
                fltOverallDrift = \
                    numpy.sum(arrChecked[self.dOverall]) / \
                    len(self.listChecked)
                fltFrameDrift = \
                    numpy.sum(arrChecked[self.dFrame]) / \
                    len(self.listChecked)
                fltEndToEndDrift = \
                    numpy.sum(arrChecked[self.dEnd]) / \
                    len(self.listChecked)
                fltMaxDistance = \
                    numpy.sum(arrChecked[self.dMax]) / \
                    len(self.listChecked)
                fltMaxDistanceZero = \
                    numpy.sum(arrChecked[self.dMaxFirst]) / \
                    len(self.listChecked)

            # Fill the widgets
            self.leAllMicNumber.setText(
                '{:d}'.format(self.lsFiles.count())
                )
            self.leMicChecked.setText(
                '{:d}'.format(len(self.listChecked))
                )
            self.leAllOverallDrift.setText(
                '{:f}'.format(fltOverallDrift)
                )
            self.leAllFrameDrift.setText(
                '{:f}'.format(fltFrameDrift)
                )
            self.leAllEndToEndDrift.setText(
                '{:f}'.format(fltEndToEndDrift)
                )
            self.leAllMaxDistance.setText(
                '{:f}'.format(fltMaxDistance)
                )
            self.leAllMaxDistanceZero.setText(
                '{:f}'.format(fltMaxDistanceZero)
                )

    def _load_from_list(self):
        """Load shift files from list"""

        self.fileName = str(QtGui.QFileDialog.getOpenFileName(
            directory=os.getcwd(),
            options=QtGui.QFileDialog.DontUseNativeDialog,
            filter='Text files (*.txt)'
            ))
        # Abort if empty
        if self.fileName == '':
            return None
        else:
            assert(self.fileName != '')
        if self._check_list_or_file(self.fileName) == 'list':
            try:
                listOfShiftFiles = numpy.genfromtxt(self.fileName, dtype=None, unpack=True)
            except TypeError:
                message = QtGui.QMessageBox(self)
                message.setText('Empty File:\n{0}'.format(self.fileName))
                message.exec_()
            except ValueError:
                message = QtGui.QMessageBox(self)
                message.setText('File is not valid, only one column allowed:\n{0}'.format(self.fileName))
                message.exec_()
            else:
                if len(numpy.shape(listOfShiftFiles)) > 1:
                    message = QtGui.QMessageBox(self)
                    message.setText('Too many columns. Expected one column:\n{0}'.format(self.fileName))
                    message.exec_()
                else:
                    if numpy.ndim(listOfShiftFiles) == 0:
                        listOfShiftFiles = [listOfShiftFiles[()]]
                    else:
                        listOfShiftFiles = list(listOfShiftFiles)
                    self._fill_gui(inputlist=listOfShiftFiles)
        else:
            print('Error: Input is no valid list')

    def _load_from_pattern(self):
        """Load files from pattern"""

        filePattern = str(self.lePattern.text()).replace('\n','')
        if self._check_list_or_file(filePattern) == 'file':
            self._fill_gui(inputfile=filePattern)
        else:
            message = QtGui.QMessageBox(self)
            message.setText('Not valid unblur shift files:\n{0}'.format(filePattern))
            message.exec_()

    def _write_selection(self):
        """Write the selected micrographs to a file"""

        # Get output path and file name
        outputPath = QtGui.QFileDialog.getExistingDirectory(
            directory=os.getcwd(),
            options=QtGui.QFileDialog.DontUseNativeDialog
            )
        if not outputPath:
            return None

        # Output names
        outputPrefix = str(self.leOutputName.text()).replace('.txt', '')
        outputSelected = '{0}/{1}_selected.txt'.format(outputPath, outputPrefix)
        outputDiscarded = '{0}/{1}_discarded.txt'.format(outputPath, outputPrefix)
        shiftSelected = '{0}/{1}_shift_selected.txt'.format(outputPath, outputPrefix)
        shiftDiscarded = '{0}/{1}_shift_discarded.txt'.format(outputPath, outputPrefix)

        # Check if output files already exists
        if os.path.exists(outputSelected) or \
                os.path.exists(outputDiscarded) or\
                os.path.exists(shiftSelected) or \
                os.path.exists(shiftDiscarded):
            # Ask the user if he really wants to continue
            warningBox = QtGui.QMessageBox(self)
            warningBox.setStandardButtons(
                QtGui.QMessageBox.No | QtGui.QMessageBox.Yes
                )
            warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
            warningBox.setText(
                    'Do you really want to save?\n\n' +
                    '{0}\nor\n{1}\nor\n{2}\nor\n{3}\nalready exists.'.format(
                        outputSelected,
                        outputDiscarded,
                        shiftSelected,
                        shiftDiscarded
                        )
                )
            warningBox.exec_()

            if warningBox.result() == QtGui.QMessageBox.No:
                return None

        with open(outputSelected, 'w') as w:
            with open(shiftSelected, 'w') as f:
                for name in sorted(self.listChecked):
                    arrCurrentEntry = self.arrData[self.arrData[self.dFile] == name]
                    micName = arrCurrentEntry[self.dMic][0]\
                        .split('/')[-1]\
                        .replace('_temp.mrc', '_sum.mrc')\
                        .replace('_frames.mrc', '_sum.mrc')
                    w.write(
                        '{:s}\n'.format(
                            micName
                            )
                        )
                    f.write(
                        '{:s}\n'.format(
                            os.path.relpath(arrCurrentEntry[self.dFileRaw][0])
                            )
                        )

        with open(outputDiscarded, 'w') as w:
            with open(shiftDiscarded, 'w') as f:
                for name in sorted(self.listUnchecked):
                    arrCurrentEntry = self.arrData[
                            self.arrData[self.dFile] == name
                            ]
                    micName = arrCurrentEntry[self.dMic][0]\
                        .split('/')[-1]\
                        .replace('_temp.mrc', '_sum.mrc')\
                        .replace('_frames.mrc', '_sum.mrc')
                    w.write(
                        '{:s}\n'.format(
                            micName
                            )
                        )
                    f.write(
                        '{:s}\n'.format(
                            os.path.relpath(arrCurrentEntry[self.dFileRaw][0])
                            )
                        )

        # Ask the user if he really wants to continue
        warningBox = QtGui.QMessageBox(self)
        warningBox.setStandardButtons(
            QtGui.QMessageBox.Yes
            )
        warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
        warningBox.setText(
                'Selection saved to:\n\n' +
                '{0}\n{1}\n{2}\n{3}\n\n'.format(
                    outputSelected,
                    outputDiscarded,
                    shiftSelected,
                    shiftDiscarded
                    ) +
                'You selected {0} of {1} micrographs ({2}%).\n'.format(
                    len(self.listChecked),
                    len(self.listChecked) + len(self.listUnchecked),
                    100 * len(self.listChecked) / (len(self.listChecked) + len(self.listUnchecked))
                    ) +
                'You discarded {0} of {1} micrographs ({2}%).'.format(
                    len(self.listUnchecked),
                    len(self.listChecked) + len(self.listUnchecked),
                    100 * len(self.listUnchecked) / (len(self.listChecked) + len(self.listUnchecked))
                    )
            )
        warningBox.exec_()

    @QtCore.pyqtSlot()
    def _write_settings(self):
        """Save settings"""

        # Output file name
        strSaveName = str(QtGui.QFileDialog.getSaveFileName(
            directory=os.getcwd(),
            options=QtGui.QFileDialog.DontUseNativeDialog
            ))

        # Save data to file
        if strSaveName != '':
            with open(strSaveName, 'w') as f:
                f.write('# Frames\n')
                f.write('{:d} {:d} {:d}\n'.format(
                    self.idxFirstFrame, self.idxLastFrame, self.intFrames
                    ))
                f.write('# DType\n')
                for row in self.listDType:
                    for element in row:
                        f.write('{0} '.format(element))
                    f.write('\n')
                f.write('# Array\n')
                for row in self.arrData:
                    for element in row:
                        f.write('{0} '.format(element))
                    f.write('\n')
                f.write('# Checked\n')
                if self.listChecked:
                    for element in self.listChecked:
                        f.write('{0} '.format(element))
                else:
                    f.write('No checked files')
                f.write('\n')
                f.write('# Unchecked\n')
                if self.listUnchecked:
                    for element in self.listUnchecked:
                        f.write('{0} '.format(element))
                else:
                    f.write('No unchecked files')
                f.write('\n')
                f.write('# listCoordX\n')
                for element in self.listCoordX:
                    f.write('{0} '.format(element))
                f.write('\n')
                f.write('# listCoordY\n')
                for element in self.listCoordY:
                    f.write('{0} '.format(element))
                f.write('\n')
                f.write('# listFrames\n')
                for element in self.listFrames:
                    f.write('{0} '.format(element))
                f.write('\n')
                f.write('# listAngles\n')
                for element in self.listAngles:
                    f.write('{0} '.format(element))
                f.write('\n')
                f.write('# dictTranslate\n')
                for key in self.dictTranslate:
                    if key == self.modeOverall:
                        f.write('{0} {1} {2}\n'.format(key, 'None', self.dictTranslate[key]))
                    else:
                        f.write('{0} {1}\n'.format(key, self.dictTranslate[key]))
                f.write('# dictThresh\n')
                for key in self.dictTranslate:
                    if key == self.modeOverall:
                        f.write('{0} {1} '.format(key, 'None'))
                    else:
                        f.write('{0} '.format(key))
                    for element in self.dictThresh[key]:
                        f.write('{0} '.format(element))
                    f.write('\n')
                f.write('# General\n')
                f.write('{0} {1}\n'.format(self.varOldStartGeneral, self.varOldStopGeneral))
                f.write('# MicNumber\n')
                for entry in self.arrMicNumber:
                    f.write('{0}\n'.format(entry))
                f.write('# End\n')

    @QtCore.pyqtSlot()
    def _load_settings(self):
        """Load settings"""

        # Warning box when refreshing frames
        warningBox = QtGui.QMessageBox(self)
        warningBox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        warningBox.setDefaultButton(QtGui.QMessageBox.Yes)
        warningBox.setText(
            'Not saved changes will be lost when loading drift data!\n' +
            'Do you really want to continue?'
            )
        warningBox.exec_()

        if warningBox.result() == QtGui.QMessageBox.No:
            return None

        # Input file
        strLoadName = str(QtGui.QFileDialog.getOpenFileName(
            directory=os.getcwd(),
            options=QtGui.QFileDialog.DontUseNativeDialog
            ))

        # If not cancel
        if strLoadName == '':
            return None
        else:
            assert(strLoadName != '')

        # Reset variables
        self._disconnect_events()
        self._set_variables()
        self._connect_events()

        try:
            arrThresh, arrGeneral = self._load_settings_new(strLoadName=strLoadName)
        except TypeError:
            arrThresh, arrGeneral = self._load_settings_old_v2(strLoadName=strLoadName)
        except IndexError:
            print('Detected old settings file! Save the settings now again to convert it to the latest version.')
            arrThresh, arrGeneral = self._load_settings_old_v1(strLoadName=strLoadName)

        # Refresh everything
        setChecked = set(self.listChecked)

        # Fill the list widget, and return the selection
        self.lsFiles.clear()
        for file in self.arrData[self.dFile]:
            newItem = QtGui.QListWidgetItem(file)
            newItem.setFlags(self.newItemFlags)
            if file in setChecked:
                newItem.setCheckState(Qt.Checked)
            else:
                newItem.setCheckState(Qt.Unchecked)
            self.lsFiles.addItem(newItem)

        # Refresh GUI representation
        self._refresh_calculations(goon=True)

        for row in arrThresh:
            listElement = []
            for index in xrange(2, len(row)):
                listElement.append(row[index])
            if row[0] != self.modeOverall:
                self.dictThresh.update({'{0} {1}'.format(row[0], row[1]): listElement})
            else:
                self.dictThresh.update({'{0}'.format(row[0]): listElement})

        # Refresh histograms
        self._fill_widgets(mode=self.modeOverall)
        self._refresh_histogram_widget(
            event='start', mode=self.modeOverall
            )
        self._refresh_histogram_widget(
            event='stop', mode=self.modeOverall
            )
        for index in range(self.cbFrame.count()):
            self.cbFrame.setCurrentIndex(index)
            self._fill_widgets(mode=self.modeFrame)
            self._refresh_histogram_widget(
                event='start', mode=self.modeFrame
                )
            self._refresh_histogram_widget(
                event='stop', mode=self.modeFrame
                )
        self.cbFrame.setCurrentIndex(0)
        for index in range(self.cbAngle.count()):
            self.cbAngle.setCurrentIndex(index)
            self._fill_widgets(mode=self.modeAngle)
            self._refresh_histogram_widget(
                event='start', mode=self.modeAngle
                )
            self._refresh_histogram_widget(
                event='stop', mode=self.modeAngle
                )
        self.cbAngle.setCurrentIndex(0)

        self.leStartGeneral.setText(arrGeneral[0])
        self.leStopGeneral.setText(arrGeneral[1])
        self._save_selection()

        # Enable everything and color all black
        self._default_color()
        self._enable_all()


    def _load_settings_new(self, strLoadName):
        """Load the settings in the new way"""

        separate_list = []
        idxFrames = 0
        idxDType = 1
        idxArray = 2
        idxChecked = 3
        idxUnchecked = 4
        idxCoordX = 5
        idxCoordY = 6
        idxFrame = 7
        idxAngles = 8
        idxTranslate = 9
        idxThresh = 10
        idxGeneral = 11
        idxMicNumber = 12
        idxEnd = 13
        with open(strLoadName, 'r') as f:
            # Append line of keywords
            for linenumber, line in enumerate(f):
                if line == '# Frames\n' or line == '# DType\n' or \
                        line == '# Array\n' or line == '# Checked\n' or  \
                        line == '# Unchecked\n' or line == '# listCoordX\n' or \
                        line == '# listCoordY\n' or line == '# listFrames\n' or \
                        line == '# listAngles\n' or line == '# dictTranslate\n' or \
                        line == '# dictThresh\n' or line == '# General\n' or \
                        line == '# MicNumber\n' or line == '# End\n':
                    separate_list.append(linenumber)

        # Fill frame widgets
        end = separate_list[idxEnd] - separate_list[idxDType] + 1
        self.idxFirstFrame, self.idxLastFrame, self.intFrames = numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxFrames] + 1, skip_footer=end,
            comments='$', dtype='<i8'
            )
        self.leFrameStart.setText('{:d}'.format(self.idxFirstFrame))
        self.leFrameStop.setText('{:d}'.format(self.idxLastFrame))

        # Get Dtype
        end = separate_list[idxEnd] - separate_list[idxArray] + 1
        arrName, arrType = numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxDType] + 1,
            skip_footer=end, comments='$',
            dtype=None, unpack=True
            )

        self.listDType = []
        for name, dtype in zip(arrName, arrType):
            self.listDType.append((name, dtype))

        # Load data array
        end = separate_list[idxEnd] - separate_list[idxChecked] + 1
        self.arrData = numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxArray] + 1,
            skip_footer=end, comments='$',
            dtype=self.listDType
            )
        self.arrData = numpy.atleast_1d(self.arrData)

        # Load checked list
        end = separate_list[idxEnd] - separate_list[idxUnchecked] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxChecked] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.listChecked = [str(loaded_data)]
        else:
            self.listChecked = list(loaded_data)

        if self.listChecked[0] == 'No':
            self.listChecked = []

        # Load unchecked list
        end = separate_list[idxEnd] - separate_list[idxCoordX] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxUnchecked] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.listUnchecked = [str(loaded_data)]
        else:
            self.listUnchecked = list(loaded_data)

        if self.listUnchecked[0] == 'No':
            self.listUnchecked = []

        # Load coord x list
        end = separate_list[idxEnd] - separate_list[idxCoordY] + 1
        self.listCoordX = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxCoordX] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load coord y list
        end = separate_list[idxEnd] - separate_list[idxFrame] + 1
        self.listCoordY = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxCoordY] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load frame list
        end = separate_list[idxEnd] - separate_list[idxAngles] + 1
        self.listFrames = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxFrame] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load angle list
        end = separate_list[idxEnd] - separate_list[idxTranslate] + 1
        self.listAngles = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxAngles] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load Translate dictionary
        end = separate_list[idxEnd] - separate_list[idxThresh] + 1
        arrNames, arrNumber, arrTrans = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxTranslate] + 1,
            skip_footer=end, comments='$',
            dtype=None, unpack=True
            ))

        for name, number, trans in zip(arrNames, arrNumber, arrTrans):
            if name != self.modeOverall:
                self.dictTranslate.update({'{0} {1}'.format(name, number): trans})
            else:
                self.dictTranslate.update({name: trans})

        # Load Thresh dictionary
        end = separate_list[idxEnd] - separate_list[idxGeneral] + 1
        arrThresh = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxThresh] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load general settings
        end = separate_list[idxEnd] - separate_list[idxMicNumber] + 1
        arrGeneral = list(numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxGeneral] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load MicNumber settings
        end = separate_list[idxEnd] - separate_list[idxEnd] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skip_header=separate_list[idxMicNumber] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.arrMicNumber = [str(loaded_data)]
        else:
            self.arrMicNumber = list(loaded_data)

        # Check if there are less than 5
        if len(self.arrMicNumber) <= 5:
            self.varAnalyzeOne = True
            print(
                '\nWarning: !!!! Only few shift files selected, ' +
                'so plots of all micrographs could not work as expected. !!!!\n'
                )

        return arrThresh, arrGeneral


    def _load_settings_old_v2(self, strLoadName):
        """Load the settings in the new way"""

        separate_list = []
        idxFrames = 0
        idxDType = 1
        idxArray = 2
        idxChecked = 3
        idxUnchecked = 4
        idxCoordX = 5
        idxCoordY = 6
        idxFrame = 7
        idxAngles = 8
        idxTranslate = 9
        idxThresh = 10
        idxGeneral = 11
        idxMicNumber = 12
        idxEnd = 13
        with open(strLoadName, 'r') as f:
            # Append line of keywords
            for linenumber, line in enumerate(f):
                if line == '# Frames\n' or line == '# DType\n' or \
                        line == '# Array\n' or line == '# Checked\n' or  \
                        line == '# Unchecked\n' or line == '# listCoordX\n' or \
                        line == '# listCoordY\n' or line == '# listFrames\n' or \
                        line == '# listAngles\n' or line == '# dictTranslate\n' or \
                        line == '# dictThresh\n' or line == '# General\n' or \
                        line == '# MicNumber\n' or line == '# End\n':
                    separate_list.append(linenumber)

        # Fill frame widgets
        end = separate_list[idxEnd] - separate_list[idxDType] + 1
        self.idxFirstFrame, self.idxLastFrame, self.intFrames = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxFrames] + 1, skip_footer=end,
            comments='$', dtype='<i8'
            )
        self.leFrameStart.setText('{:d}'.format(self.idxFirstFrame))
        self.leFrameStop.setText('{:d}'.format(self.idxLastFrame))

        # Get Dtype
        end = separate_list[idxEnd] - separate_list[idxArray] + 1
        arrName, arrType = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxDType] + 1,
            skip_footer=end, comments='$',
            dtype=None, unpack=True
            )

        self.listDType = []
        for name, dtype in zip(arrName, arrType):
            self.listDType.append((name, dtype))

        # Load data array
        end = separate_list[idxEnd] - separate_list[idxChecked] + 1
        self.arrData = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxArray] + 1,
            skip_footer=end, comments='$',
            dtype=self.listDType
            )
        self.arrData = numpy.atleast_1d(self.arrData)

        # Load checked list
        end = separate_list[idxEnd] - separate_list[idxUnchecked] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxChecked] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.listChecked = [str(loaded_data)]
        else:
            self.listChecked = list(loaded_data)

        if self.listChecked[0] == 'No':
            self.listChecked = []

        # Load unchecked list
        end = separate_list[idxEnd] - separate_list[idxCoordX] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxUnchecked] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.listUnchecked = [str(loaded_data)]
        else:
            self.listUnchecked = list(loaded_data)

        if self.listUnchecked[0] == 'No':
            self.listUnchecked = []

        # Load coord x list
        end = separate_list[idxEnd] - separate_list[idxCoordY] + 1
        self.listCoordX = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxCoordX] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load coord y list
        end = separate_list[idxEnd] - separate_list[idxFrame] + 1
        self.listCoordY = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxCoordY] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load frame list
        end = separate_list[idxEnd] - separate_list[idxAngles] + 1
        self.listFrames = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxFrame] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load angle list
        end = separate_list[idxEnd] - separate_list[idxTranslate] + 1
        self.listAngles = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxAngles] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load Translate dictionary
        end = separate_list[idxEnd] - separate_list[idxThresh] + 1
        arrNames, arrNumber, arrTrans = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxTranslate] + 1,
            skip_footer=end, comments='$',
            dtype=None, unpack=True
            ))

        for name, number, trans in zip(arrNames, arrNumber, arrTrans):
            if name != self.modeOverall:
                self.dictTranslate.update({'{0} {1}'.format(name, number): trans})
            else:
                self.dictTranslate.update({name: trans})

        # Load Thresh dictionary
        end = separate_list[idxEnd] - separate_list[idxGeneral] + 1
        arrThresh = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxThresh] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load general settings
        end = separate_list[idxEnd] - separate_list[idxMicNumber] + 1
        arrGeneral = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxGeneral] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load MicNumber settings
        end = separate_list[idxEnd] - separate_list[idxEnd] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxMicNumber] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.arrMicNumber = [str(loaded_data)]
        else:
            self.arrMicNumber = list(loaded_data)

        # Check if there are less than 5
        if len(self.arrMicNumber) <= 5:
            self.varAnalyzeOne = True
            print(
                '\nWarning: !!!! Only few shift files selected, ' +
                'so plots of all micrographs could not work as expected. !!!!\n'
                )

        return arrThresh, arrGeneral


    def _load_settings_old_v1(self, strLoadName):
        """Load the settings in the new way"""

        separate_list = []
        idxFrames = 0
        idxDType = 1
        idxArray = 2
        idxChecked = 3
        idxUnchecked = 4
        idxCoordX = 5
        idxCoordY = 6
        idxFrame = 7
        idxAngles = 8
        idxTranslate = 9
        idxThresh = 10
        idxGeneral = 11
        idxEnd = 12

        # Check for the keywords
        with open(strLoadName, 'r') as f:
            # Append line of keywords
            for linenumber, line in enumerate(f):
                if line == '# Frames\n' or line == '# DType\n' or \
                        line == '# Array\n' or line == '# Checked\n' or  \
                        line == '# Unchecked\n' or line == '# listCoordX\n' or \
                        line == '# listCoordY\n' or line == '# listFrames\n' or \
                        line == '# listAngles\n' or line == '# dictTranslate\n' or \
                        line == '# dictThresh\n' or line == '# General\n' or \
                        line == '# End\n':
                    separate_list.append(linenumber)

        # Fill frame widgets
        end = separate_list[idxEnd] - separate_list[idxDType] + 1
        self.idxFirstFrame, self.idxLastFrame, self.intFrames = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxFrames] + 1, skip_footer=end,
            comments='$', dtype='<i8'
            )
        self.leFrameStart.setText('{:d}'.format(self.idxFirstFrame))
        self.leFrameStop.setText('{:d}'.format(self.idxLastFrame))

        # Get Dtype
        end = separate_list[idxEnd] - separate_list[idxArray] + 1
        arrName, arrType = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxDType] + 1,
            skip_footer=end, comments='$',
            dtype=None, unpack=True
            )

        self.listDType = []
        for name, dtype in zip(arrName, arrType):
            self.listDType.append((name, dtype))

        # Load data array
        end = separate_list[idxEnd] - separate_list[idxChecked] + 1
        self.arrData = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxArray] + 1,
            skip_footer=end, comments='$',
            dtype=self.listDType
            )
        self.arrData = numpy.atleast_1d(self.arrData)

        # Load checked list
        end = separate_list[idxEnd] - separate_list[idxUnchecked] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxChecked] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.listChecked = [str(loaded_data)]
        else:
            self.listChecked = list(loaded_data)

        if self.listChecked[0] == 'No':
            self.listChecked = []

        # Load unchecked list
        end = separate_list[idxEnd] - separate_list[idxCoordX] + 1
        loaded_data = numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxUnchecked] + 1,
            skip_footer=end, comments='$',
            dtype=None
            )
        if len(numpy.shape(loaded_data)) == 0:
            self.listUnchecked = [str(loaded_data)]
        else:
            self.listUnchecked = list(loaded_data)

        if self.listUnchecked[0] == 'No':
            self.listUnchecked = []

        # Load coord x list
        end = separate_list[idxEnd] - separate_list[idxCoordY] + 1
        self.listCoordX = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxCoordX] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load coord y list
        end = separate_list[idxEnd] - separate_list[idxFrame] + 1
        self.listCoordY = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxCoordY] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load frame list
        end = separate_list[idxEnd] - separate_list[idxAngles] + 1
        self.listFrames = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxFrame] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load angle list
        end = separate_list[idxEnd] - separate_list[idxTranslate] + 1
        self.listAngles = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxAngles] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load Translate dictionary
        end = separate_list[idxEnd] - separate_list[idxThresh] + 1
        arrNames, arrNumber, arrTrans = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxTranslate] + 1,
            skip_footer=end, comments='$',
            dtype=None, unpack=True
            ))

        for name, number, trans in zip(arrNames, arrNumber, arrTrans):
            if name != self.modeOverall:
                self.dictTranslate.update({'{0} {1}'.format(name, number): trans})
            else:
                self.dictTranslate.update({name: trans})

        # Load Thresh dictionary
        end = separate_list[idxEnd] - separate_list[idxGeneral] + 1
        arrThresh = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxThresh] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        # Load general settings
        end = separate_list[idxEnd] - separate_list[idxEnd] + 1
        arrGeneral = list(numpy.genfromtxt(
            strLoadName,
            skiprows=separate_list[idxGeneral] + 1,
            skip_footer=end, comments='$',
            dtype=None
            ))

        print(
            'Warning: Could not identify micrograph serial number due to an old settings file.\n' +
            'X-axis of the plots do not represent' +
            'the serial micrograph number.'
            )
        self.arrMicNumber = numpy.arange(len(self.arrData))

        return arrThresh, arrGeneral


    def _show_about(self):
        """Show the about info"""

        # Generate a about message box
        about = QtGui.QMessageBox()
        about.setText(
            """
            sxgui_drift for analyzing drift parameters
            made by Unblur
            Copyright (C) 2016  Markus Stabrin
            (markus.stabrin@mpi-dortmund.mpg.de)

            This software is issued under a joint
            BSD/GNU license. You may use the
            source code in this file under either
            license.
            However, note that the
            complete EMAN2 and SPARX software
            packages have some GPL dependencies,
            so you are responsible for compliance
            with the licenses of these packages
            if you opt to use BSD licensing.
            The warranty disclaimer below holds
            in either instance.

            This complete copyright notice must
            be included in any revised version of
            the source code.
            Additional authorship citations may
            be added, but existing author
            citations must be preserved.

            This program is free software;
            you can redistribute it and/or modify
            it under the terms of the
            GNU General Public License as
            published by
            the Free Software Foundation;
             either version 2 of the License, or
            (at your option) any later version.

            This program is distributed in
            the hope that it will be useful,
            but WITHOUT ANY WARRANTY;
            without even the implied warranty of
            MERCHANTABILITY or
            FITNESS FOR A PARTICULAR PURPOSE.
            See the GNU General Public License
            for more details.
            """
            )
        about.exec_()

    def changeEvent(self, event):
        if event.type() == QtCore.QEvent.WindowStateChange:
            if self.isMinimized():
                self.dictVisible.update({
                    self.msPlotDrift: [
                        self.msPlotDrift.isVisible(),
                        self.msPlotDrift.rect(),
                        self.msPlotDrift.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msPlotFrame: [
                        self.msPlotFrame.isVisible(),
                        self.msPlotFrame.rect(),
                        self.msPlotFrame.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msPlotAngle: [
                        self.msPlotAngle.isVisible(),
                        self.msPlotAngle.rect(),
                        self.msPlotAngle.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msAllPlotFrameAvg: [
                        self.msAllPlotFrameAvg.isVisible(),
                        self.msAllPlotFrameAvg.rect(),
                        self.msAllPlotFrameAvg.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msAllPlotDrift: [
                        self.msAllPlotDrift.isVisible(),
                        self.msAllPlotDrift.rect(),
                        self.msAllPlotDrift.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msAllPlotFrame: [
                        self.msAllPlotFrame.isVisible(),
                        self.msAllPlotFrame.rect(),
                        self.msAllPlotFrame.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msAllPlotAngle: [
                        self.msAllPlotAngle.isVisible(),
                        self.msAllPlotAngle.rect(),
                        self.msAllPlotAngle.pos()
                        ]
                    })
                self.dictVisible.update({
                    self.msAllPlotPerMic: [
                        self.msAllPlotPerMic.isVisible(),
                        self.msAllPlotPerMic.rect(),
                        self.msAllPlotPerMic.pos()
                        ]
                    })

                for key in self.dictVisible:
                    if self.dictVisible[key]:
                        key.hide()

                self.minimized = True
            elif self.minimized:
                self.minimized = False
                for key in self.dictVisible:
                    if self.dictVisible[key][self.idxVisible]:
                        key.setGeometry(self.dictVisible[key][self.idxRect])
                        key.move(self.dictVisible[key][self.idxPos])
                        key.activateWindow()
                        key.show()

                if self.isVisible():
                    self.raise_()
                    self.activateWindow()

        elif event.type() == QtCore.QEvent.ActivationChange:

            if self.isActiveWindow():
                for key in self.dictVisible:
                    if self.dictVisible[key][self.idxVisible]:
                        key.raise_()
                self.raise_()
                self.activateWindow()

    def closeEvent(self, event):
        """Change the closeEvent to close the application cleanly"""

        # Ignore the incomming event when pressing the "X" and
        # quit the application instead
        event.ignore()
        print('Bye Bye!')
        QtCore.QCoreApplication.instance().quit()


def _main():
    app = QtGui.QApplication(sys.argv)

    if len(sys.argv) > 3:
        print(
            'Too many arguments!\n'
            'Please start the GUI without argument or use \n'
            'quote marks around the wildcard "example-*_shift.txt".\n'
            'Usage: sxgui_unblur.py "shift_files" "list_file"\n'
            )
        return None
    elif len(sys.argv) == 3:
        shiftFiles, shiftList = sys.argv[-2:]
        msApp = SXDriftUnblur(inputfile=shiftFiles, inputlist=shiftList)
    elif len(sys.argv) == 2:
        shiftFiles = sys.argv[-1]
        msApp = SXDriftUnblur(inputfile=shiftFiles)
    else:
        arg = None
        msApp = SXDriftUnblur(arg)

    msApp.show()

    app.exec_()


if __name__ == '__main__':
    _main()
