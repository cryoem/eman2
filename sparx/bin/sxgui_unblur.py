#!/usr/bin/env python
# sxgui_drift for analyzing drift parameters made by Unblur
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
import os
import sys
import glob
import numpy
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg

# The name changed in newer versions of matplotlib
try:
    from matplotlib.backends.backend_qt4agg \
        import NavigationToolbar2QTAgg as NavigationToolbar2QT
except:
    from matplotlib.backends.backend_qt4agg \
        import NavigationToolbar2QT as NavigationToolbar2QT

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
        MSMainWidget.setObjectName(_fromUtf8("MSMainWidget"))
        MSMainWidget.resize(1316, 773)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MSMainWidget.sizePolicy().hasHeightForWidth())
        MSMainWidget.setSizePolicy(sizePolicy)
        self.layoutWidget = QtGui.QWidget(MSMainWidget)
        self.layoutWidget.setGeometry(QtCore.QRect(-20, 0, 1367, 771))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_2.setMargin(0)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.layoutCentralWidget = QtGui.QHBoxLayout()
        self.layoutCentralWidget.setObjectName(_fromUtf8("layoutCentralWidget"))
        self.layoutInput = QtGui.QVBoxLayout()
        self.layoutInput.setObjectName(_fromUtf8("layoutInput"))
        spacerItem = QtGui.QSpacerItem(20, 4, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.layoutInput.addItem(spacerItem)
        self.horizontalLayout_12 = QtGui.QHBoxLayout()
        self.horizontalLayout_12.setObjectName(_fromUtf8("horizontalLayout_12"))
        self.pbFindDir = QtGui.QPushButton(self.layoutWidget)
        self.pbFindDir.setEnabled(True)
        self.pbFindDir.setCheckable(False)
        self.pbFindDir.setObjectName(_fromUtf8("pbFindDir"))
        self.horizontalLayout_12.addWidget(self.pbFindDir)
        self.label_7 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.horizontalLayout_12.addWidget(self.label_7)
        self.leSuffix = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leSuffix.sizePolicy().hasHeightForWidth())
        self.leSuffix.setSizePolicy(sizePolicy)
        self.leSuffix.setObjectName(_fromUtf8("leSuffix"))
        self.horizontalLayout_12.addWidget(self.leSuffix)
        self.layoutInput.addLayout(self.horizontalLayout_12)
        self.lsFiles = QtGui.QListWidget(self.layoutWidget)
        self.lsFiles.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lsFiles.sizePolicy().hasHeightForWidth())
        self.lsFiles.setSizePolicy(sizePolicy)
        self.lsFiles.setMinimumSize(QtCore.QSize(200, 0))
        self.lsFiles.setObjectName(_fromUtf8("lsFiles"))
        self.layoutInput.addWidget(self.lsFiles)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.pbSelectAll = QtGui.QPushButton(self.layoutWidget)
        self.pbSelectAll.setEnabled(False)
        self.pbSelectAll.setObjectName(_fromUtf8("pbSelectAll"))
        self.horizontalLayout_5.addWidget(self.pbSelectAll)
        self.pbInvert = QtGui.QPushButton(self.layoutWidget)
        self.pbInvert.setEnabled(False)
        self.pbInvert.setObjectName(_fromUtf8("pbInvert"))
        self.horizontalLayout_5.addWidget(self.pbInvert)
        self.layoutInput.addLayout(self.horizontalLayout_5)
        self.layoutCentralWidget.addLayout(self.layoutInput)
        self.layoutInfo = QtGui.QVBoxLayout()
        self.layoutInfo.setObjectName(_fromUtf8("layoutInfo"))
        self.label_6 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.layoutInfo.addWidget(self.label_6)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.label_4 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.horizontalLayout_4.addWidget(self.label_4)
        self.leCurrentMicName = QtGui.QLineEdit(self.layoutWidget)
        self.leCurrentMicName.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(22)
        sizePolicy.setHeightForWidth(self.leCurrentMicName.sizePolicy().hasHeightForWidth())
        self.leCurrentMicName.setSizePolicy(sizePolicy)
        self.leCurrentMicName.setMinimumSize(QtCore.QSize(300, 0))
        self.leCurrentMicName.setText(_fromUtf8(""))
        self.leCurrentMicName.setReadOnly(True)
        self.leCurrentMicName.setObjectName(_fromUtf8("leCurrentMicName"))
        self.horizontalLayout_4.addWidget(self.leCurrentMicName)
        self.layoutInfo.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.label_3 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout_3.addWidget(self.label_3)
        self.leCurrentOverallDrift = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leCurrentOverallDrift.sizePolicy().hasHeightForWidth())
        self.leCurrentOverallDrift.setSizePolicy(sizePolicy)
        self.leCurrentOverallDrift.setReadOnly(True)
        self.leCurrentOverallDrift.setObjectName(_fromUtf8("leCurrentOverallDrift"))
        self.horizontalLayout_3.addWidget(self.leCurrentOverallDrift)
        self.layoutInfo.addLayout(self.horizontalLayout_3)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label = QtGui.QLabel(self.layoutWidget)
        self.label.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.leCurrentFrameDrift = QtGui.QLineEdit(self.layoutWidget)
        self.leCurrentFrameDrift.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leCurrentFrameDrift.sizePolicy().hasHeightForWidth())
        self.leCurrentFrameDrift.setSizePolicy(sizePolicy)
        self.leCurrentFrameDrift.setReadOnly(True)
        self.leCurrentFrameDrift.setObjectName(_fromUtf8("leCurrentFrameDrift"))
        self.horizontalLayout.addWidget(self.leCurrentFrameDrift)
        self.layoutInfo.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label_2 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_2.addWidget(self.label_2)
        self.leCurrentEndToEndDrift = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leCurrentEndToEndDrift.sizePolicy().hasHeightForWidth())
        self.leCurrentEndToEndDrift.setSizePolicy(sizePolicy)
        self.leCurrentEndToEndDrift.setReadOnly(True)
        self.leCurrentEndToEndDrift.setObjectName(_fromUtf8("leCurrentEndToEndDrift"))
        self.horizontalLayout_2.addWidget(self.leCurrentEndToEndDrift)
        self.layoutInfo.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_27 = QtGui.QHBoxLayout()
        self.horizontalLayout_27.setObjectName(_fromUtf8("horizontalLayout_27"))
        self.label_30 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_30.sizePolicy().hasHeightForWidth())
        self.label_30.setSizePolicy(sizePolicy)
        self.label_30.setObjectName(_fromUtf8("label_30"))
        self.horizontalLayout_27.addWidget(self.label_30)
        self.leCurrentMaxDistance = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leCurrentMaxDistance.sizePolicy().hasHeightForWidth())
        self.leCurrentMaxDistance.setSizePolicy(sizePolicy)
        self.leCurrentMaxDistance.setReadOnly(True)
        self.leCurrentMaxDistance.setObjectName(_fromUtf8("leCurrentMaxDistance"))
        self.horizontalLayout_27.addWidget(self.leCurrentMaxDistance)
        self.layoutInfo.addLayout(self.horizontalLayout_27)
        self.horizontalLayout_29 = QtGui.QHBoxLayout()
        self.horizontalLayout_29.setObjectName(_fromUtf8("horizontalLayout_29"))
        self.label_32 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_32.sizePolicy().hasHeightForWidth())
        self.label_32.setSizePolicy(sizePolicy)
        self.label_32.setObjectName(_fromUtf8("label_32"))
        self.horizontalLayout_29.addWidget(self.label_32)
        self.leCurrentMaxDistanceZero = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leCurrentMaxDistanceZero.sizePolicy().hasHeightForWidth())
        self.leCurrentMaxDistanceZero.setSizePolicy(sizePolicy)
        self.leCurrentMaxDistanceZero.setReadOnly(True)
        self.leCurrentMaxDistanceZero.setObjectName(_fromUtf8("leCurrentMaxDistanceZero"))
        self.horizontalLayout_29.addWidget(self.leCurrentMaxDistanceZero)
        self.layoutInfo.addLayout(self.horizontalLayout_29)
        self.label_8 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_8.sizePolicy().hasHeightForWidth())
        self.label_8.setSizePolicy(sizePolicy)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.layoutInfo.addWidget(self.label_8)
        self.horizontalLayout_11 = QtGui.QHBoxLayout()
        self.horizontalLayout_11.setObjectName(_fromUtf8("horizontalLayout_11"))
        self.label_13 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_13.sizePolicy().hasHeightForWidth())
        self.label_13.setSizePolicy(sizePolicy)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.horizontalLayout_11.addWidget(self.label_13)
        self.leAllMicNumber = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leAllMicNumber.sizePolicy().hasHeightForWidth())
        self.leAllMicNumber.setSizePolicy(sizePolicy)
        self.leAllMicNumber.setReadOnly(True)
        self.leAllMicNumber.setObjectName(_fromUtf8("leAllMicNumber"))
        self.horizontalLayout_11.addWidget(self.leAllMicNumber)
        self.layoutInfo.addLayout(self.horizontalLayout_11)
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName(_fromUtf8("horizontalLayout_20"))
        self.label_29 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_29.sizePolicy().hasHeightForWidth())
        self.label_29.setSizePolicy(sizePolicy)
        self.label_29.setObjectName(_fromUtf8("label_29"))
        self.horizontalLayout_20.addWidget(self.label_29)
        self.leMicChecked = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leMicChecked.sizePolicy().hasHeightForWidth())
        self.leMicChecked.setSizePolicy(sizePolicy)
        self.leMicChecked.setReadOnly(True)
        self.leMicChecked.setObjectName(_fromUtf8("leMicChecked"))
        self.horizontalLayout_20.addWidget(self.leMicChecked)
        self.layoutInfo.addLayout(self.horizontalLayout_20)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.label_9 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.horizontalLayout_7.addWidget(self.label_9)
        self.leAllOverallDrift = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leAllOverallDrift.sizePolicy().hasHeightForWidth())
        self.leAllOverallDrift.setSizePolicy(sizePolicy)
        self.leAllOverallDrift.setReadOnly(True)
        self.leAllOverallDrift.setObjectName(_fromUtf8("leAllOverallDrift"))
        self.horizontalLayout_7.addWidget(self.leAllOverallDrift)
        self.layoutInfo.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.label_12 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.horizontalLayout_10.addWidget(self.label_12)
        self.leAllFrameDrift = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leAllFrameDrift.sizePolicy().hasHeightForWidth())
        self.leAllFrameDrift.setSizePolicy(sizePolicy)
        self.leAllFrameDrift.setReadOnly(True)
        self.leAllFrameDrift.setObjectName(_fromUtf8("leAllFrameDrift"))
        self.horizontalLayout_10.addWidget(self.leAllFrameDrift)
        self.layoutInfo.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.label_10 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_10.sizePolicy().hasHeightForWidth())
        self.label_10.setSizePolicy(sizePolicy)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.horizontalLayout_8.addWidget(self.label_10)
        self.leAllEndToEndDrift = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leAllEndToEndDrift.sizePolicy().hasHeightForWidth())
        self.leAllEndToEndDrift.setSizePolicy(sizePolicy)
        self.leAllEndToEndDrift.setReadOnly(True)
        self.leAllEndToEndDrift.setObjectName(_fromUtf8("leAllEndToEndDrift"))
        self.horizontalLayout_8.addWidget(self.leAllEndToEndDrift)
        self.layoutInfo.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_26 = QtGui.QHBoxLayout()
        self.horizontalLayout_26.setObjectName(_fromUtf8("horizontalLayout_26"))
        self.label_28 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_28.sizePolicy().hasHeightForWidth())
        self.label_28.setSizePolicy(sizePolicy)
        self.label_28.setObjectName(_fromUtf8("label_28"))
        self.horizontalLayout_26.addWidget(self.label_28)
        self.leAllMaxDistance = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leAllMaxDistance.sizePolicy().hasHeightForWidth())
        self.leAllMaxDistance.setSizePolicy(sizePolicy)
        self.leAllMaxDistance.setReadOnly(True)
        self.leAllMaxDistance.setObjectName(_fromUtf8("leAllMaxDistance"))
        self.horizontalLayout_26.addWidget(self.leAllMaxDistance)
        self.layoutInfo.addLayout(self.horizontalLayout_26)
        self.horizontalLayout_28 = QtGui.QHBoxLayout()
        self.horizontalLayout_28.setObjectName(_fromUtf8("horizontalLayout_28"))
        self.label_31 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_31.sizePolicy().hasHeightForWidth())
        self.label_31.setSizePolicy(sizePolicy)
        self.label_31.setObjectName(_fromUtf8("label_31"))
        self.horizontalLayout_28.addWidget(self.label_31)
        self.leAllMaxDistanceZero = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leAllMaxDistanceZero.sizePolicy().hasHeightForWidth())
        self.leAllMaxDistanceZero.setSizePolicy(sizePolicy)
        self.leAllMaxDistanceZero.setReadOnly(True)
        self.leAllMaxDistanceZero.setObjectName(_fromUtf8("leAllMaxDistanceZero"))
        self.horizontalLayout_28.addWidget(self.leAllMaxDistanceZero)
        self.layoutInfo.addLayout(self.horizontalLayout_28)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_16 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.verticalLayout.addWidget(self.label_16)
        self.cbSort = QtGui.QComboBox(self.layoutWidget)
        self.cbSort.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cbSort.sizePolicy().hasHeightForWidth())
        self.cbSort.setSizePolicy(sizePolicy)
        self.cbSort.setObjectName(_fromUtf8("cbSort"))
        self.cbSort.addItem(_fromUtf8(""))
        self.cbSort.addItem(_fromUtf8(""))
        self.cbSort.addItem(_fromUtf8(""))
        self.cbSort.addItem(_fromUtf8(""))
        self.cbSort.addItem(_fromUtf8(""))
        self.cbSort.addItem(_fromUtf8(""))
        self.verticalLayout.addWidget(self.cbSort)
        self.horizontalLayout_21 = QtGui.QHBoxLayout()
        self.horizontalLayout_21.setObjectName(_fromUtf8("horizontalLayout_21"))
        self.chDescending = QtGui.QCheckBox(self.layoutWidget)
        self.chDescending.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chDescending.sizePolicy().hasHeightForWidth())
        self.chDescending.setSizePolicy(sizePolicy)
        self.chDescending.setObjectName(_fromUtf8("chDescending"))
        self.horizontalLayout_21.addWidget(self.chDescending)
        self.chSortSelected = QtGui.QCheckBox(self.layoutWidget)
        self.chSortSelected.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chSortSelected.sizePolicy().hasHeightForWidth())
        self.chSortSelected.setSizePolicy(sizePolicy)
        self.chSortSelected.setObjectName(_fromUtf8("chSortSelected"))
        self.horizontalLayout_21.addWidget(self.chSortSelected)
        self.verticalLayout.addLayout(self.horizontalLayout_21)
        self.layoutInfo.addLayout(self.verticalLayout)
        self.layoutCentralWidget.addLayout(self.layoutInfo)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        spacerItem1 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_3.addItem(spacerItem1)
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.label_14 = QtGui.QLabel(self.layoutWidget)
        self.label_14.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_14.sizePolicy().hasHeightForWidth())
        self.label_14.setSizePolicy(sizePolicy)
        self.label_14.setBaseSize(QtCore.QSize(0, 0))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.verticalLayout_4.addWidget(self.label_14)
        self.chPlotDriftMic = QtGui.QCheckBox(self.layoutWidget)
        self.chPlotDriftMic.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chPlotDriftMic.sizePolicy().hasHeightForWidth())
        self.chPlotDriftMic.setSizePolicy(sizePolicy)
        self.chPlotDriftMic.setObjectName(_fromUtf8("chPlotDriftMic"))
        self.verticalLayout_4.addWidget(self.chPlotDriftMic)
        self.chPlotFrameMic = QtGui.QCheckBox(self.layoutWidget)
        self.chPlotFrameMic.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chPlotFrameMic.sizePolicy().hasHeightForWidth())
        self.chPlotFrameMic.setSizePolicy(sizePolicy)
        self.chPlotFrameMic.setObjectName(_fromUtf8("chPlotFrameMic"))
        self.verticalLayout_4.addWidget(self.chPlotFrameMic)
        self.chPlotAngleMic = QtGui.QCheckBox(self.layoutWidget)
        self.chPlotAngleMic.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chPlotAngleMic.sizePolicy().hasHeightForWidth())
        self.chPlotAngleMic.setSizePolicy(sizePolicy)
        self.chPlotAngleMic.setObjectName(_fromUtf8("chPlotAngleMic"))
        self.verticalLayout_4.addWidget(self.chPlotAngleMic)
        spacerItem2 = QtGui.QSpacerItem(20, 19, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_4.addItem(spacerItem2)
        self.horizontalLayout_9.addLayout(self.verticalLayout_4)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.label_15 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.verticalLayout_5.addWidget(self.label_15)
        self.chAverageDriftPerFrame = QtGui.QCheckBox(self.layoutWidget)
        self.chAverageDriftPerFrame.setEnabled(False)
        self.chAverageDriftPerFrame.setObjectName(_fromUtf8("chAverageDriftPerFrame"))
        self.verticalLayout_5.addWidget(self.chAverageDriftPerFrame)
        self.chPlotDrift = QtGui.QCheckBox(self.layoutWidget)
        self.chPlotDrift.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chPlotDrift.sizePolicy().hasHeightForWidth())
        self.chPlotDrift.setSizePolicy(sizePolicy)
        self.chPlotDrift.setObjectName(_fromUtf8("chPlotDrift"))
        self.verticalLayout_5.addWidget(self.chPlotDrift)
        self.chPlotFrame = QtGui.QCheckBox(self.layoutWidget)
        self.chPlotFrame.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chPlotFrame.sizePolicy().hasHeightForWidth())
        self.chPlotFrame.setSizePolicy(sizePolicy)
        self.chPlotFrame.setObjectName(_fromUtf8("chPlotFrame"))
        self.verticalLayout_5.addWidget(self.chPlotFrame)
        self.chPlotAngle = QtGui.QCheckBox(self.layoutWidget)
        self.chPlotAngle.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chPlotAngle.sizePolicy().hasHeightForWidth())
        self.chPlotAngle.setSizePolicy(sizePolicy)
        self.chPlotAngle.setObjectName(_fromUtf8("chPlotAngle"))
        self.verticalLayout_5.addWidget(self.chPlotAngle)
        self.horizontalLayout_9.addLayout(self.verticalLayout_5)
        self.verticalLayout_6.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_32 = QtGui.QHBoxLayout()
        self.horizontalLayout_32.setObjectName(_fromUtf8("horizontalLayout_32"))
        self.label_5 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.horizontalLayout_32.addWidget(self.label_5)
        self.leFrameStart = QtGui.QLineEdit(self.layoutWidget)
        self.leFrameStart.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leFrameStart.sizePolicy().hasHeightForWidth())
        self.leFrameStart.setSizePolicy(sizePolicy)
        self.leFrameStart.setObjectName(_fromUtf8("leFrameStart"))
        self.horizontalLayout_32.addWidget(self.leFrameStart)
        self.label_11 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.horizontalLayout_32.addWidget(self.label_11)
        self.leFrameStop = QtGui.QLineEdit(self.layoutWidget)
        self.leFrameStop.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leFrameStop.sizePolicy().hasHeightForWidth())
        self.leFrameStop.setSizePolicy(sizePolicy)
        self.leFrameStop.setObjectName(_fromUtf8("leFrameStop"))
        self.horizontalLayout_32.addWidget(self.leFrameStop)
        self.verticalLayout_6.addLayout(self.horizontalLayout_32)
        self.verticalLayout_3.addLayout(self.verticalLayout_6)
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.label_21 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_21.sizePolicy().hasHeightForWidth())
        self.label_21.setSizePolicy(sizePolicy)
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.verticalLayout_7.addWidget(self.label_21)
        self.horizontalLayout_13 = QtGui.QHBoxLayout()
        self.horizontalLayout_13.setObjectName(_fromUtf8("horizontalLayout_13"))
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setObjectName(_fromUtf8("verticalLayout_8"))
        self.label_22 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_22.sizePolicy().hasHeightForWidth())
        self.label_22.setSizePolicy(sizePolicy)
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.verticalLayout_8.addWidget(self.label_22)
        self.label_23 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_23.sizePolicy().hasHeightForWidth())
        self.label_23.setSizePolicy(sizePolicy)
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.verticalLayout_8.addWidget(self.label_23)
        self.label_24 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_24.sizePolicy().hasHeightForWidth())
        self.label_24.setSizePolicy(sizePolicy)
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.verticalLayout_8.addWidget(self.label_24)
        self.label_25 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_25.sizePolicy().hasHeightForWidth())
        self.label_25.setSizePolicy(sizePolicy)
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.verticalLayout_8.addWidget(self.label_25)
        self.horizontalLayout_13.addLayout(self.verticalLayout_8)
        self.verticalLayout_9 = QtGui.QVBoxLayout()
        self.verticalLayout_9.setObjectName(_fromUtf8("verticalLayout_9"))
        self.leStartOverall = QtGui.QLineEdit(self.layoutWidget)
        self.leStartOverall.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartOverall.sizePolicy().hasHeightForWidth())
        self.leStartOverall.setSizePolicy(sizePolicy)
        self.leStartOverall.setObjectName(_fromUtf8("leStartOverall"))
        self.verticalLayout_9.addWidget(self.leStartOverall)
        self.leStopOverall = QtGui.QLineEdit(self.layoutWidget)
        self.leStopOverall.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopOverall.sizePolicy().hasHeightForWidth())
        self.leStopOverall.setSizePolicy(sizePolicy)
        self.leStopOverall.setObjectName(_fromUtf8("leStopOverall"))
        self.verticalLayout_9.addWidget(self.leStopOverall)
        self.leStartSaveOverall = QtGui.QLineEdit(self.layoutWidget)
        self.leStartSaveOverall.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartSaveOverall.sizePolicy().hasHeightForWidth())
        self.leStartSaveOverall.setSizePolicy(sizePolicy)
        self.leStartSaveOverall.setReadOnly(True)
        self.leStartSaveOverall.setObjectName(_fromUtf8("leStartSaveOverall"))
        self.verticalLayout_9.addWidget(self.leStartSaveOverall)
        self.leStopSaveOverall = QtGui.QLineEdit(self.layoutWidget)
        self.leStopSaveOverall.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopSaveOverall.sizePolicy().hasHeightForWidth())
        self.leStopSaveOverall.setSizePolicy(sizePolicy)
        self.leStopSaveOverall.setReadOnly(True)
        self.leStopSaveOverall.setObjectName(_fromUtf8("leStopSaveOverall"))
        self.verticalLayout_9.addWidget(self.leStopSaveOverall)
        self.horizontalLayout_13.addLayout(self.verticalLayout_9)
        self.verticalLayout_10 = QtGui.QVBoxLayout()
        self.verticalLayout_10.setObjectName(_fromUtf8("verticalLayout_10"))
        spacerItem3 = QtGui.QSpacerItem(20, 26, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_10.addItem(spacerItem3)
        self.pbSaveOverall = QtGui.QPushButton(self.layoutWidget)
        self.pbSaveOverall.setEnabled(False)
        self.pbSaveOverall.setObjectName(_fromUtf8("pbSaveOverall"))
        self.verticalLayout_10.addWidget(self.pbSaveOverall)
        self.chOverallCriterion = QtGui.QCheckBox(self.layoutWidget)
        self.chOverallCriterion.setEnabled(False)
        self.chOverallCriterion.setObjectName(_fromUtf8("chOverallCriterion"))
        self.verticalLayout_10.addWidget(self.chOverallCriterion)
        spacerItem4 = QtGui.QSpacerItem(20, 26, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_10.addItem(spacerItem4)
        self.horizontalLayout_13.addLayout(self.verticalLayout_10)
        self.verticalLayout_7.addLayout(self.horizontalLayout_13)
        self.verticalLayout_3.addLayout(self.verticalLayout_7)
        self.tabWidget = QtGui.QTabWidget(self.layoutWidget)
        self.tabWidget.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tab.sizePolicy().hasHeightForWidth())
        self.tab.setSizePolicy(sizePolicy)
        self.tab.setObjectName(_fromUtf8("tab"))
        self.layoutWidget_4 = QtGui.QWidget(self.tab)
        self.layoutWidget_4.setGeometry(QtCore.QRect(0, 0, 441, 141))
        self.layoutWidget_4.setObjectName(_fromUtf8("layoutWidget_4"))
        self.verticalLayout_15 = QtGui.QVBoxLayout(self.layoutWidget_4)
        self.verticalLayout_15.setMargin(0)
        self.verticalLayout_15.setObjectName(_fromUtf8("verticalLayout_15"))
        self.label_37 = QtGui.QLabel(self.layoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_37.sizePolicy().hasHeightForWidth())
        self.label_37.setSizePolicy(sizePolicy)
        self.label_37.setObjectName(_fromUtf8("label_37"))
        self.verticalLayout_15.addWidget(self.label_37)
        self.horizontalLayout_15 = QtGui.QHBoxLayout()
        self.horizontalLayout_15.setObjectName(_fromUtf8("horizontalLayout_15"))
        self.verticalLayout_16 = QtGui.QVBoxLayout()
        self.verticalLayout_16.setObjectName(_fromUtf8("verticalLayout_16"))
        self.label_38 = QtGui.QLabel(self.layoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_38.sizePolicy().hasHeightForWidth())
        self.label_38.setSizePolicy(sizePolicy)
        self.label_38.setObjectName(_fromUtf8("label_38"))
        self.verticalLayout_16.addWidget(self.label_38)
        self.label_39 = QtGui.QLabel(self.layoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_39.sizePolicy().hasHeightForWidth())
        self.label_39.setSizePolicy(sizePolicy)
        self.label_39.setObjectName(_fromUtf8("label_39"))
        self.verticalLayout_16.addWidget(self.label_39)
        self.horizontalLayout_15.addLayout(self.verticalLayout_16)
        self.verticalLayout_17 = QtGui.QVBoxLayout()
        self.verticalLayout_17.setObjectName(_fromUtf8("verticalLayout_17"))
        self.leStartGeneral = QtGui.QLineEdit(self.layoutWidget_4)
        self.leStartGeneral.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartGeneral.sizePolicy().hasHeightForWidth())
        self.leStartGeneral.setSizePolicy(sizePolicy)
        self.leStartGeneral.setStyleSheet(_fromUtf8(""))
        self.leStartGeneral.setText(_fromUtf8(""))
        self.leStartGeneral.setObjectName(_fromUtf8("leStartGeneral"))
        self.verticalLayout_17.addWidget(self.leStartGeneral)
        self.leStopGeneral = QtGui.QLineEdit(self.layoutWidget_4)
        self.leStopGeneral.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopGeneral.sizePolicy().hasHeightForWidth())
        self.leStopGeneral.setSizePolicy(sizePolicy)
        self.leStopGeneral.setText(_fromUtf8(""))
        self.leStopGeneral.setObjectName(_fromUtf8("leStopGeneral"))
        self.verticalLayout_17.addWidget(self.leStopGeneral)
        self.horizontalLayout_15.addLayout(self.verticalLayout_17)
        self.verticalLayout_18 = QtGui.QVBoxLayout()
        self.verticalLayout_18.setObjectName(_fromUtf8("verticalLayout_18"))
        self.pbSaveGeneral = QtGui.QPushButton(self.layoutWidget_4)
        self.pbSaveGeneral.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbSaveGeneral.sizePolicy().hasHeightForWidth())
        self.pbSaveGeneral.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(False)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(50)
        font.setStrikeOut(False)
        self.pbSaveGeneral.setFont(font)
        self.pbSaveGeneral.setAutoRepeat(False)
        self.pbSaveGeneral.setObjectName(_fromUtf8("pbSaveGeneral"))
        self.verticalLayout_18.addWidget(self.pbSaveGeneral)
        self.chGeneralCriterion = QtGui.QCheckBox(self.layoutWidget_4)
        self.chGeneralCriterion.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chGeneralCriterion.sizePolicy().hasHeightForWidth())
        self.chGeneralCriterion.setSizePolicy(sizePolicy)
        self.chGeneralCriterion.setTristate(False)
        self.chGeneralCriterion.setObjectName(_fromUtf8("chGeneralCriterion"))
        self.verticalLayout_18.addWidget(self.chGeneralCriterion)
        self.horizontalLayout_15.addLayout(self.verticalLayout_18)
        self.verticalLayout_15.addLayout(self.horizontalLayout_15)
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tab_2.sizePolicy().hasHeightForWidth())
        self.tab_2.setSizePolicy(sizePolicy)
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.layoutWidget_3 = QtGui.QWidget(self.tab_2)
        self.layoutWidget_3.setGeometry(QtCore.QRect(0, 10, 441, 161))
        self.layoutWidget_3.setObjectName(_fromUtf8("layoutWidget_3"))
        self.verticalLayout_11 = QtGui.QVBoxLayout(self.layoutWidget_3)
        self.verticalLayout_11.setMargin(0)
        self.verticalLayout_11.setObjectName(_fromUtf8("verticalLayout_11"))
        self.label_26 = QtGui.QLabel(self.layoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_26.sizePolicy().hasHeightForWidth())
        self.label_26.setSizePolicy(sizePolicy)
        self.label_26.setObjectName(_fromUtf8("label_26"))
        self.verticalLayout_11.addWidget(self.label_26)
        self.horizontalLayout_14 = QtGui.QHBoxLayout()
        self.horizontalLayout_14.setObjectName(_fromUtf8("horizontalLayout_14"))
        self.verticalLayout_12 = QtGui.QVBoxLayout()
        self.verticalLayout_12.setObjectName(_fromUtf8("verticalLayout_12"))
        self.label_33 = QtGui.QLabel(self.layoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_33.sizePolicy().hasHeightForWidth())
        self.label_33.setSizePolicy(sizePolicy)
        self.label_33.setObjectName(_fromUtf8("label_33"))
        self.verticalLayout_12.addWidget(self.label_33)
        self.label_34 = QtGui.QLabel(self.layoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_34.sizePolicy().hasHeightForWidth())
        self.label_34.setSizePolicy(sizePolicy)
        self.label_34.setObjectName(_fromUtf8("label_34"))
        self.verticalLayout_12.addWidget(self.label_34)
        self.label_35 = QtGui.QLabel(self.layoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_35.sizePolicy().hasHeightForWidth())
        self.label_35.setSizePolicy(sizePolicy)
        self.label_35.setObjectName(_fromUtf8("label_35"))
        self.verticalLayout_12.addWidget(self.label_35)
        self.label_36 = QtGui.QLabel(self.layoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_36.sizePolicy().hasHeightForWidth())
        self.label_36.setSizePolicy(sizePolicy)
        self.label_36.setObjectName(_fromUtf8("label_36"))
        self.verticalLayout_12.addWidget(self.label_36)
        self.horizontalLayout_14.addLayout(self.verticalLayout_12)
        self.verticalLayout_13 = QtGui.QVBoxLayout()
        self.verticalLayout_13.setObjectName(_fromUtf8("verticalLayout_13"))
        self.leStartFrame = QtGui.QLineEdit(self.layoutWidget_3)
        self.leStartFrame.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartFrame.sizePolicy().hasHeightForWidth())
        self.leStartFrame.setSizePolicy(sizePolicy)
        self.leStartFrame.setObjectName(_fromUtf8("leStartFrame"))
        self.verticalLayout_13.addWidget(self.leStartFrame)
        self.leStopFrame = QtGui.QLineEdit(self.layoutWidget_3)
        self.leStopFrame.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopFrame.sizePolicy().hasHeightForWidth())
        self.leStopFrame.setSizePolicy(sizePolicy)
        self.leStopFrame.setObjectName(_fromUtf8("leStopFrame"))
        self.verticalLayout_13.addWidget(self.leStopFrame)
        self.leStartSaveFrame = QtGui.QLineEdit(self.layoutWidget_3)
        self.leStartSaveFrame.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartSaveFrame.sizePolicy().hasHeightForWidth())
        self.leStartSaveFrame.setSizePolicy(sizePolicy)
        self.leStartSaveFrame.setReadOnly(True)
        self.leStartSaveFrame.setObjectName(_fromUtf8("leStartSaveFrame"))
        self.verticalLayout_13.addWidget(self.leStartSaveFrame)
        self.leStopSaveFrame = QtGui.QLineEdit(self.layoutWidget_3)
        self.leStopSaveFrame.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopSaveFrame.sizePolicy().hasHeightForWidth())
        self.leStopSaveFrame.setSizePolicy(sizePolicy)
        self.leStopSaveFrame.setReadOnly(True)
        self.leStopSaveFrame.setObjectName(_fromUtf8("leStopSaveFrame"))
        self.verticalLayout_13.addWidget(self.leStopSaveFrame)
        self.horizontalLayout_14.addLayout(self.verticalLayout_13)
        self.verticalLayout_14 = QtGui.QVBoxLayout()
        self.verticalLayout_14.setObjectName(_fromUtf8("verticalLayout_14"))
        self.cbFrame = QtGui.QComboBox(self.layoutWidget_3)
        self.cbFrame.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cbFrame.sizePolicy().hasHeightForWidth())
        self.cbFrame.setSizePolicy(sizePolicy)
        self.cbFrame.setObjectName(_fromUtf8("cbFrame"))
        self.verticalLayout_14.addWidget(self.cbFrame)
        self.pbSaveFrame = QtGui.QPushButton(self.layoutWidget_3)
        self.pbSaveFrame.setEnabled(False)
        self.pbSaveFrame.setObjectName(_fromUtf8("pbSaveFrame"))
        self.verticalLayout_14.addWidget(self.pbSaveFrame)
        self.chFrameCriterion = QtGui.QCheckBox(self.layoutWidget_3)
        self.chFrameCriterion.setEnabled(False)
        self.chFrameCriterion.setObjectName(_fromUtf8("chFrameCriterion"))
        self.verticalLayout_14.addWidget(self.chFrameCriterion)
        spacerItem5 = QtGui.QSpacerItem(20, 27, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_14.addItem(spacerItem5)
        self.horizontalLayout_14.addLayout(self.verticalLayout_14)
        self.verticalLayout_11.addLayout(self.horizontalLayout_14)
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.verticalLayout_3.addWidget(self.tabWidget)
        self.verticalLayout_39 = QtGui.QVBoxLayout()
        self.verticalLayout_39.setObjectName(_fromUtf8("verticalLayout_39"))
        self.label_59 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_59.sizePolicy().hasHeightForWidth())
        self.label_59.setSizePolicy(sizePolicy)
        self.label_59.setObjectName(_fromUtf8("label_59"))
        self.verticalLayout_39.addWidget(self.label_59)
        self.horizontalLayout_24 = QtGui.QHBoxLayout()
        self.horizontalLayout_24.setObjectName(_fromUtf8("horizontalLayout_24"))
        self.verticalLayout_40 = QtGui.QVBoxLayout()
        self.verticalLayout_40.setObjectName(_fromUtf8("verticalLayout_40"))
        self.label_60 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_60.sizePolicy().hasHeightForWidth())
        self.label_60.setSizePolicy(sizePolicy)
        self.label_60.setObjectName(_fromUtf8("label_60"))
        self.verticalLayout_40.addWidget(self.label_60)
        self.label_61 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_61.sizePolicy().hasHeightForWidth())
        self.label_61.setSizePolicy(sizePolicy)
        self.label_61.setObjectName(_fromUtf8("label_61"))
        self.verticalLayout_40.addWidget(self.label_61)
        self.label_62 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_62.sizePolicy().hasHeightForWidth())
        self.label_62.setSizePolicy(sizePolicy)
        self.label_62.setObjectName(_fromUtf8("label_62"))
        self.verticalLayout_40.addWidget(self.label_62)
        self.label_63 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_63.sizePolicy().hasHeightForWidth())
        self.label_63.setSizePolicy(sizePolicy)
        self.label_63.setObjectName(_fromUtf8("label_63"))
        self.verticalLayout_40.addWidget(self.label_63)
        self.horizontalLayout_24.addLayout(self.verticalLayout_40)
        self.verticalLayout_41 = QtGui.QVBoxLayout()
        self.verticalLayout_41.setObjectName(_fromUtf8("verticalLayout_41"))
        self.leStartAngle = QtGui.QLineEdit(self.layoutWidget)
        self.leStartAngle.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartAngle.sizePolicy().hasHeightForWidth())
        self.leStartAngle.setSizePolicy(sizePolicy)
        self.leStartAngle.setObjectName(_fromUtf8("leStartAngle"))
        self.verticalLayout_41.addWidget(self.leStartAngle)
        self.leStopAngle = QtGui.QLineEdit(self.layoutWidget)
        self.leStopAngle.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopAngle.sizePolicy().hasHeightForWidth())
        self.leStopAngle.setSizePolicy(sizePolicy)
        self.leStopAngle.setObjectName(_fromUtf8("leStopAngle"))
        self.verticalLayout_41.addWidget(self.leStopAngle)
        self.leStartSaveAngle = QtGui.QLineEdit(self.layoutWidget)
        self.leStartSaveAngle.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStartSaveAngle.sizePolicy().hasHeightForWidth())
        self.leStartSaveAngle.setSizePolicy(sizePolicy)
        self.leStartSaveAngle.setReadOnly(True)
        self.leStartSaveAngle.setObjectName(_fromUtf8("leStartSaveAngle"))
        self.verticalLayout_41.addWidget(self.leStartSaveAngle)
        self.leStopSaveAngle = QtGui.QLineEdit(self.layoutWidget)
        self.leStopSaveAngle.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leStopSaveAngle.sizePolicy().hasHeightForWidth())
        self.leStopSaveAngle.setSizePolicy(sizePolicy)
        self.leStopSaveAngle.setReadOnly(True)
        self.leStopSaveAngle.setObjectName(_fromUtf8("leStopSaveAngle"))
        self.verticalLayout_41.addWidget(self.leStopSaveAngle)
        self.horizontalLayout_24.addLayout(self.verticalLayout_41)
        self.verticalLayout_42 = QtGui.QVBoxLayout()
        self.verticalLayout_42.setObjectName(_fromUtf8("verticalLayout_42"))
        self.cbAngle = QtGui.QComboBox(self.layoutWidget)
        self.cbAngle.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cbAngle.sizePolicy().hasHeightForWidth())
        self.cbAngle.setSizePolicy(sizePolicy)
        self.cbAngle.setObjectName(_fromUtf8("cbAngle"))
        self.verticalLayout_42.addWidget(self.cbAngle)
        self.pbSaveAngle = QtGui.QPushButton(self.layoutWidget)
        self.pbSaveAngle.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbSaveAngle.sizePolicy().hasHeightForWidth())
        self.pbSaveAngle.setSizePolicy(sizePolicy)
        self.pbSaveAngle.setObjectName(_fromUtf8("pbSaveAngle"))
        self.verticalLayout_42.addWidget(self.pbSaveAngle)
        self.chAngleCriterion = QtGui.QCheckBox(self.layoutWidget)
        self.chAngleCriterion.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chAngleCriterion.sizePolicy().hasHeightForWidth())
        self.chAngleCriterion.setSizePolicy(sizePolicy)
        self.chAngleCriterion.setObjectName(_fromUtf8("chAngleCriterion"))
        self.verticalLayout_42.addWidget(self.chAngleCriterion)
        self.pbUncheckCriterion = QtGui.QPushButton(self.layoutWidget)
        self.pbUncheckCriterion.setEnabled(False)
        self.pbUncheckCriterion.setObjectName(_fromUtf8("pbUncheckCriterion"))
        self.verticalLayout_42.addWidget(self.pbUncheckCriterion)
        self.horizontalLayout_24.addLayout(self.verticalLayout_42)
        self.verticalLayout_39.addLayout(self.horizontalLayout_24)
        self.verticalLayout_3.addLayout(self.verticalLayout_39)
        self.pbApply = QtGui.QPushButton(self.layoutWidget)
        self.pbApply.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbApply.sizePolicy().hasHeightForWidth())
        self.pbApply.setSizePolicy(sizePolicy)
        self.pbApply.setMinimumSize(QtCore.QSize(450, 0))
        self.pbApply.setObjectName(_fromUtf8("pbApply"))
        self.verticalLayout_3.addWidget(self.pbApply)
        self.layoutCentralWidget.addLayout(self.verticalLayout_3)
        self.layoutCentralWidget.setStretch(2, 1)
        self.verticalLayout_2.addLayout(self.layoutCentralWidget)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        spacerItem6 = QtGui.QSpacerItem(20, 1, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem6)
        self.label_27 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_27.sizePolicy().hasHeightForWidth())
        self.label_27.setSizePolicy(sizePolicy)
        self.label_27.setObjectName(_fromUtf8("label_27"))
        self.horizontalLayout_6.addWidget(self.label_27)
        self.leOutputName = QtGui.QLineEdit(self.layoutWidget)
        self.leOutputName.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leOutputName.sizePolicy().hasHeightForWidth())
        self.leOutputName.setSizePolicy(sizePolicy)
        self.leOutputName.setMinimumSize(QtCore.QSize(500, 0))
        self.leOutputName.setObjectName(_fromUtf8("leOutputName"))
        self.horizontalLayout_6.addWidget(self.leOutputName)
        self.pbSaveSelected = QtGui.QPushButton(self.layoutWidget)
        self.pbSaveSelected.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbSaveSelected.sizePolicy().hasHeightForWidth())
        self.pbSaveSelected.setSizePolicy(sizePolicy)
        self.pbSaveSelected.setObjectName(_fromUtf8("pbSaveSelected"))
        self.horizontalLayout_6.addWidget(self.pbSaveSelected)
        spacerItem7 = QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem7)
        self.pbSaveSettings = QtGui.QPushButton(self.layoutWidget)
        self.pbSaveSettings.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbSaveSettings.sizePolicy().hasHeightForWidth())
        self.pbSaveSettings.setSizePolicy(sizePolicy)
        self.pbSaveSettings.setObjectName(_fromUtf8("pbSaveSettings"))
        self.horizontalLayout_6.addWidget(self.pbSaveSettings)
        self.pbLoadSettings = QtGui.QPushButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbLoadSettings.sizePolicy().hasHeightForWidth())
        self.pbLoadSettings.setSizePolicy(sizePolicy)
        self.pbLoadSettings.setObjectName(_fromUtf8("pbLoadSettings"))
        self.horizontalLayout_6.addWidget(self.pbLoadSettings)
        spacerItem8 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem8)
        self.pbAbout = QtGui.QPushButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pbAbout.sizePolicy().hasHeightForWidth())
        self.pbAbout.setSizePolicy(sizePolicy)
        self.pbAbout.setObjectName(_fromUtf8("pbAbout"))
        self.horizontalLayout_6.addWidget(self.pbAbout)
        spacerItem9 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem9)
        self.verticalLayout_2.addLayout(self.horizontalLayout_6)

        self.retranslateUi(MSMainWidget)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MSMainWidget)

    def retranslateUi(self, MSMainWidget):
        MSMainWidget.setWindowTitle(_translate("MSMainWidget", "sxgui_drift", None))
        self.pbFindDir.setText(_translate("MSMainWidget", "Shift Directory", None))
        self.label_7.setText(_translate("MSMainWidget", "Name Pattern:", None))
        self.leSuffix.setText(_translate("MSMainWidget", "*_shift.txt", None))
        self.pbSelectAll.setText(_translate("MSMainWidget", "Select All", None))
        self.pbInvert.setText(_translate("MSMainWidget", "Invert Selection", None))
        self.label_6.setText(_translate("MSMainWidget", "Info of Current Entry:", None))
        self.label_4.setText(_translate("MSMainWidget", "   Micrograph Name", None))
        self.label_3.setText(_translate("MSMainWidget", "   Overall Drift", None))
        self.label.setText(_translate("MSMainWidget", "   Drift per Frame", None))
        self.label_2.setText(_translate("MSMainWidget", "   End to End Length", None))
        self.label_30.setText(_translate("MSMainWidget", "   Maximum Distance between Frames", None))
        self.label_32.setText(_translate("MSMainWidget", "   Maximum Distance from Frame Start Frame", None))
        self.label_8.setText(_translate("MSMainWidget", "Drift Info for selected Micrographs:", None))
        self.label_13.setText(_translate("MSMainWidget", "   Nr. of Micrographs", None))
        self.label_29.setText(_translate("MSMainWidget", "   Checked Micrographs", None))
        self.label_9.setText(_translate("MSMainWidget", "   Average Overall Drift", None))
        self.label_12.setText(_translate("MSMainWidget", "   Average Drift per Frame", None))
        self.label_10.setText(_translate("MSMainWidget", "   Average End to End Length", None))
        self.label_28.setText(_translate("MSMainWidget", "   Average Maximum Distance between Frames", None))
        self.label_31.setText(_translate("MSMainWidget", "   Average Maximum Distance from Start Frame", None))
        self.label_16.setText(_translate("MSMainWidget", "Sort Entrys:", None))
        self.cbSort.setItemText(0, _translate("MSMainWidget", "File Name", None))
        self.cbSort.setItemText(1, _translate("MSMainWidget", "Overall Drift", None))
        self.cbSort.setItemText(2, _translate("MSMainWidget", "Drift per Frame", None))
        self.cbSort.setItemText(3, _translate("MSMainWidget", "End to End Length", None))
        self.cbSort.setItemText(4, _translate("MSMainWidget", "Maximum Distance between Frames", None))
        self.cbSort.setItemText(5, _translate("MSMainWidget", "Maximum Distance from Start Frame", None))
        self.chDescending.setText(_translate("MSMainWidget", "Descending", None))
        self.chSortSelected.setText(_translate("MSMainWidget", "Sort Selected", None))
        self.label_14.setText(_translate("MSMainWidget", "  Show Plots of Current Entry:   ", None))
        self.chPlotDriftMic.setText(_translate("MSMainWidget", "Drift", None))
        self.chPlotFrameMic.setText(_translate("MSMainWidget", "Drift per Frame", None))
        self.chPlotAngleMic.setText(_translate("MSMainWidget", "Angle per Frame", None))
        self.label_15.setText(_translate("MSMainWidget", "  Show Plots of all Micrographs:", None))
        self.chAverageDriftPerFrame.setText(_translate("MSMainWidget", "Average Drift per Frame", None))
        self.chPlotDrift.setText(_translate("MSMainWidget", "Overall Drift Histogram", None))
        self.chPlotFrame.setText(_translate("MSMainWidget", "Drift per Frame Histogram", None))
        self.chPlotAngle.setText(_translate("MSMainWidget", "Angle per Frame Histogram", None))
        self.label_5.setText(_translate("MSMainWidget", "Start Frame", None))
        self.label_11.setText(_translate("MSMainWidget", "End Frame", None))
        self.label_21.setText(_translate("MSMainWidget", "Threshold Overall Drift", None))
        self.label_22.setText(_translate("MSMainWidget", "Start (A):", None))
        self.label_23.setText(_translate("MSMainWidget", "Stop (A):", None))
        self.label_24.setText(_translate("MSMainWidget", "Registered Start (A):", None))
        self.label_25.setText(_translate("MSMainWidget", "Registered Stop (A):", None))
        self.pbSaveOverall.setText(_translate("MSMainWidget", "Register", None))
        self.chOverallCriterion.setText(_translate("MSMainWidget", "Use as criterion", None))
        self.label_37.setText(_translate("MSMainWidget", "Threshold if the drift of one Frame is", None))
        self.label_38.setText(_translate("MSMainWidget", "Larger then (A):", None))
        self.label_39.setText(_translate("MSMainWidget", "Smaller then (A):", None))
        self.pbSaveGeneral.setText(_translate("MSMainWidget", "Register", None))
        self.chGeneralCriterion.setText(_translate("MSMainWidget", "Use as criterion", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MSMainWidget", "General", None))
        self.label_26.setText(_translate("MSMainWidget", "Threshold Drift per Frame", None))
        self.label_33.setText(_translate("MSMainWidget", "Start (A):", None))
        self.label_34.setText(_translate("MSMainWidget", "Stop (A):", None))
        self.label_35.setText(_translate("MSMainWidget", "Registered Start (A):", None))
        self.label_36.setText(_translate("MSMainWidget", "Registered Stop (A):", None))
        self.pbSaveFrame.setText(_translate("MSMainWidget", "Register", None))
        self.chFrameCriterion.setText(_translate("MSMainWidget", "Use as criterion", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MSMainWidget", "Per Frame", None))
        self.label_59.setText(_translate("MSMainWidget", "Threshold Angle per Frame", None))
        self.label_60.setText(_translate("MSMainWidget", "Start (A):", None))
        self.label_61.setText(_translate("MSMainWidget", "Stop (A):", None))
        self.label_62.setText(_translate("MSMainWidget", "Registered Start (A):", None))
        self.label_63.setText(_translate("MSMainWidget", "Registered Stop (A):", None))
        self.pbSaveAngle.setText(_translate("MSMainWidget", "Register", None))
        self.chAngleCriterion.setText(_translate("MSMainWidget", "Use as criterion", None))
        self.pbUncheckCriterion.setText(_translate("MSMainWidget", "Uncheck Criteria", None))
        self.pbApply.setText(_translate("MSMainWidget", "Apply registered settings marked as criterion", None))
        self.label_27.setText(_translate("MSMainWidget", "  Output Name:", None))
        self.pbSaveSelected.setText(_translate("MSMainWidget", "Save Selection", None))
        self.pbSaveSettings.setText(_translate("MSMainWidget", "Save Settings", None))
        self.pbLoadSettings.setText(_translate("MSMainWidget", "Load Settings", None))
        self.pbAbout.setText(_translate("MSMainWidget", "About", None))



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
        self.toolbar = NavigationToolbar2QT(self.canvas, self.widgetPlot)
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

    def __init__(self, inputfile=None, parent=None):
        super(SXDriftUnblur, self).__init__(parent)

        # Setup Gui Elements
        self.setupUi(self)

        # Set variables
        self._set_variables()

        # Connect events
        self._connect_events()

        if inputfile is not None:
            self._fill_gui(inputfile)

    def _set_variables(self):
        """Set instance variables"""

        # Arrays
        self.arrData = None

        # Variables
        self.strInputDir = None
        self.intNrFrames = None
        self.idxFirstFrame = None
        self.idxLastFrame = None

        # Indices Threshold and Widget
        self.idxStart = 0
        self.idxStop = 1
        self.idxStartSave = 2
        self.idxStopSave = 3
        self.idxCriterion = 4
        self.idxSaved = 5
        self.idxName = 6

        # output directory
        self.outputDir = 'unblur_GUI_output'

        # Modes
        self.modeOverall = 'Overall'
        self.modeFrame = 'Frame'
        self.modeAngle = 'Angle'
        self.modeGeneral = 'General'
        self.modeAverage = 'Average'
        self.modeDrift = 'Drift'
        self.modeDriftPerFrame = 'Drift per Frame'
        self.modeAnglePerFrame = 'Angle per Frame'

        # DType
        self.dFile = 'fileName'
        self.dMic = 'micName'
        self.dOverall = 'driftOverall'
        self.dFrame = 'driftFrame'
        self.dEnd = 'endToEndDrift'
        self.dMax = 'maxDistance'
        self.dMaxFirst = 'maxDistanceFirst'

        # Sort
        self.sortFile = 'File Name'
        self.sortOverall = 'Overall Drift'
        self.sortFrame = 'Drift per Frame'
        self.sortEnd = 'End to End Length'
        self.sortMax = 'Maximum Distance between Frames'
        self.sortMaxFirst = 'Maximum Distance from Start Frame'

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
            self.chPlotAngle: self.msAllPlotAngle
            }
        self.dictHide = {
            self.msPlotDrift: self.chPlotDriftMic,
            self.msPlotFrame: self.chPlotFrameMic,
            self.msPlotAngle: self.chPlotAngleMic,
            self.msAllPlotFrameAvg: self.chAverageDriftPerFrame,
            self.msAllPlotDrift: self.chPlotDrift,
            self.msAllPlotFrame: self.chPlotFrame,
            self.msAllPlotAngle: self.chPlotAngle
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

        # Lists
        self.listChecked = []
        self.listUnchecked = []
        self.listFile = []
        self.listCoordX = []
        self.listCoordY = []
        self.listFrames = []
        self.listAngles = []
        self.listDType = []

        # List widget flags
        self.newItemFlags = \
            Qt.ItemFlags(Qt.ItemIsSelectable) | \
            Qt.ItemFlags(Qt.ItemIsEnabled) | \
            Qt.ItemFlags(Qt.ItemIsUserCheckable)

    def _connect_events(self):
        """Connect the widgets to the events"""

        # Connect buttons
        self.pbFindDir.clicked.connect(self._find_dir)
        self.pbSelectAll.clicked.connect(self._select_all)
        self.pbInvert.clicked.connect(self._invert_selection)
        self.pbAbout.clicked.connect(self._show_about)
        self.pbSaveSelected.clicked.connect(self._write_selection)
        self.pbSaveOverall.clicked.connect(
            lambda: self._save_settings(mode=self.modeOverall)
            )
        self.pbSaveGeneral.clicked.connect(
            lambda: self._save_settings(mode=self.modeGeneral)
            )
        self.pbSaveFrame.clicked.connect(
            lambda: self._save_settings(mode=self.modeFrame)
            )
        self.pbSaveAngle.clicked.connect(
            lambda: self._save_settings(mode=self.modeAngle)
            )
        self.pbApply.clicked.connect(self._apply_setting)
        self.pbUncheckCriterion.clicked.connect(self._uncheck_angle_criterion)
        self.pbSaveSettings.clicked.connect(self._write_settings)
        self.pbLoadSettings.clicked.connect(self._load_settings)

        # Connect list widget
        self.lsFiles.itemClicked.connect(self._current_info)

        # Connect entry widgets to calculation
        self.leFrameStart.returnPressed.connect(self._refresh_calculations)
        self.leFrameStop.returnPressed.connect(self._refresh_calculations)
        self.leStartOverall.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeOverall
                )
            )
        self.leStopOverall.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeOverall
                )
            )
        self.leStartGeneral.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeGeneral
                )
            )
        self.leStopGeneral.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeGeneral
                )
            )
        self.leStartFrame.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeFrame
                )
            )
        self.leStopFrame.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeFrame
                )
            )
        self.leStartAngle.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='start', mode=self.modeAngle
                )
            )
        self.leStopAngle.returnPressed.connect(
            lambda: self._refresh_histogram_widget(
                event='stop', mode=self.modeAngle
                )
            )

        # Connect entry widgets to change color
        self.leFrameStart.textChanged.connect(
            lambda: self._change_color(
                widget=self.leFrameStart
                )
            )
        self.leFrameStop.textChanged.connect(
            lambda: self._change_color(
                widget=self.leFrameStop
                )
            )
        self.leStartOverall.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartOverall
                )
            )
        self.leStopOverall.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopOverall
                )
            )
        self.leStartGeneral.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartGeneral
                )
            )
        self.leStopGeneral.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopGeneral
                )
            )
        self.leStartFrame.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartFrame
                )
            )
        self.leStopFrame.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopFrame
                )
            )
        self.leStartAngle.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStartAngle
                )
            )
        self.leStopAngle.textChanged.connect(
            lambda: self._change_color(
                widget=self.leStopAngle
                )
            )

        # Connect check boxes
        self.chDescending.stateChanged.connect(self._apply_sorting)
        self.chSortSelected.stateChanged.connect(self._apply_sorting)
        self.chOverallCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeOverall)
            )
        self.chGeneralCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeGeneral)
            )
        self.chFrameCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeFrame)
            )
        self.chAngleCriterion.stateChanged.connect(
            lambda: self._mark_as_criterion(mode=self.modeAngle)
            )
        self.chPlotDriftMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotDriftMic)
            )
        self.chPlotFrameMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotFrameMic)
            )
        self.chPlotAngleMic.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotAngleMic)
            )
        self.chAverageDriftPerFrame.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chAverageDriftPerFrame)
            )
        self.chPlotDrift.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotDrift)
            )
        self.chPlotFrame.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotFrame)
            )
        self.chPlotAngle.stateChanged.connect(
            lambda: self._show_plot(checkbox=self.chPlotAngle)
            )

        # Connect combo boxes
        self.cbSort.activated.connect(self._apply_sorting)
        self.cbFrame.activated.connect(
            lambda: self._plot_threshold(mode=self.modeFrame, fill=True)
            )
        self.cbAngle.activated.connect(
            lambda: self._plot_threshold(mode=self.modeAngle, fill=True)
            )

        # Connect canvas close events
        self.msPlotDrift.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msPlotDrift)
            )
        self.msPlotFrame.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msPlotFrame)
            )
        self.msPlotAngle.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msPlotAngle)
            )
        self.msAllPlotFrameAvg.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotFrameAvg)
            )
        self.msAllPlotDrift.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotDrift)
            )
        self.msAllPlotFrame.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotFrame)
            )
        self.msAllPlotAngle.sigClose.connect(
            lambda: self._hide_plot(msplot=self.msAllPlotAngle)
            )

        # Connect canvas refresh frame event
        self.msAllPlotFrame.sigFrame.connect(self._refresh_frame)
        self.msAllPlotAngle.sigFrame.connect(self._refresh_frame)

        # Connect canvas refresh plot event
        self.msAllPlotDrift.sigRefresh.connect(
            self._refresh_histogram_mouse
            )
        self.msAllPlotFrame.sigRefresh.connect(
            self._refresh_histogram_mouse
            )
        self.msAllPlotAngle.sigRefresh.connect(
            self._refresh_histogram_mouse
            )

    def _enable_all(self):
        """Enable all widgets"""

        self.pbFindDir.setEnabled(True)
        self.pbSelectAll.setEnabled(True)
        self.pbInvert.setEnabled(True)
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

    def _find_dir(self):
        """Open the drift files and do the drift calculations"""

        # Find Directory
        self.strInputDir = QtGui.QFileDialog.getExistingDirectory(
            directory=os.getcwd(),
            options=QtGui.QFileDialog.DontUseNativeDialog
            )

        # If the return value is not empty, fill the gui
        if self.strInputDir != '':
            self._fill_gui()

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

    def _fill_gui(self, inputfile=None):
        """Do the Calculations and fill the GUI"""

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
        self.lsFiles.clear()
        if inputfile is None:
            strDirectory = self.strInputDir
            strSuffix = str(self.leSuffix.text())
        else:
            strDirectory = os.getcwd()
            strSuffix = inputfile
        self.listFile = glob.glob('{:s}/{:s}'.format(strDirectory, strSuffix))

        if not self.listFile:
            messageBox2 = QtGui.QMessageBox()
            messageBox2.setText('No drift files found ({:s})'.format(strSuffix))
            messageBox2.exec_()
            messageBox.hide()
            return None

        # Fill list widget
        for file in self.listFile:
            file = file.split('/')[-1]
            newItem = QtGui.QListWidgetItem(file)
            newItem.setFlags(self.newItemFlags)
            newItem.setCheckState(Qt.Checked)
            self.lsFiles.addItem(newItem)

        # Load and calculate data for the first time
        self._first_time_calculations()

        # Fill the frame widgets
        self.leFrameStart.clear()
        self.leFrameStop.clear()
        self.leFrameStart.insert('1')
        self.leFrameStop.insert('{:d}'.format(self.intFrames))

        # Do Calculations
        self._calculations(oldfirst=1, oldlast=self.intFrames)

        # Fill the GUI
        self._refresh_gui()

        # Enable all widgets
        self._enable_all()

        # Go use the GUI
        messageBox.hide()
        print('Done')

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

    def _first_time_calculations(self):
        """Calculate the Drift"""

        # Check how many frames are there
        self.intFrames = 0
        with open(self.listFile[0], 'r') as f:
            for linenumber, line in enumerate(f):
                if linenumber == 3:
                    self.intFrames = int(line.split()[-1])
                    break

        # Fill the listDType with the right number of frames and angles.
        # Add entrys to the lists and fill the translation Dictionary.
        self.listDType = [
            (self.dFile, '|S200'),
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

        # Create an empty array
        self.arrData = numpy.empty(len(self.listFile), dtype=self.listDType)

        # Load the data
        idxX = 0
        idxY = 1
        for number, file in enumerate(self.listFile):

            # Get the micrograph name
            with open(file, 'r') as f:
                self.arrData[self.dMic][number] = f.readline().split()[-1]

            # Get the file name
            self.arrData[self.dFile][number] = file.split('/')[-1]

            # Import the data
            arrCoord = numpy.genfromtxt(file, unpack=True)

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
                intBins = self.lsFiles.count() / 3
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
        intBins = self.lsFiles.count() / 3

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
        for index in xrange(self.lsFiles.count()):
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
        self.dictHide[msplot].setCheckState(Qt.Unchecked)
        msplot.hide()

    def _show_plot(self, checkbox):
        """Show the Plot Widget"""

        # If the checkbox is checked, show the related widget
        if checkbox.checkState() == Qt.Checked:
            self.dictShow[checkbox].show()
        # Otherwise hide it
        elif checkbox.checkState() == Qt.Unchecked:
            self.dictShow[checkbox].hide()

    def _select_all(self):
        """Select all entrys"""

        # Set all items to checked and save the current selection state
        for index in xrange(self.lsFiles.count()):
            self.lsFiles.item(index).setCheckState(Qt.Checked)

        # Save new selection in cache
        self._save_selection()

    def _invert_selection(self):
        """Invert Selection"""

        # Invert the selection and save the current selection state
        for index in xrange(self.lsFiles.count()):
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
            '{:s}'.format(arrCurrentEntry[self.dMic][0])
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
                arrY[number - 1] = \
                    numpy.sum(self.arrData['frame{:d}'.format(int(number))]) / \
                    self.lsFiles.count()
            strTitle = r'Average drift per Frame'
            strXLabel = r'Frame'
            strYLabel = r'Average Drift / Angstrom'
            varOutput = self.msAllPlotFrameAvg
            pylab.xlim([arrX[0] - arrX[-1] * 0.1, arrX[-1] * 1.1])

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
        for index in xrange(self.lsFiles.count()):
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

    def _write_selection(self):
        """Write the selected micrographs to a file"""

        if not os.path.exists('unblur_GUI_output'):
            os.mkdir('unblur_GUI_output')

        # Create output directory
        strOutput = '{:s}/{:s}/{:s}'.format(
            os.getcwd(),
            self.outputDir,
            str(self.leOutputName.text())
            )

        # Write output
        with open(strOutput, 'w') as f:
            for name in sorted(self.listChecked):
                arrCurrentEntry = self.arrData[self.arrData[self.dFile] == name]
                f.write('{:s}\n'.format(arrCurrentEntry[self.dMic][0]))

        # Show message with the save path
        messageBox = QtGui.QMessageBox(self)
        messageBox.setText(
            'File saved!\n{:s}'.format(strOutput)
            )
        messageBox.exec_()

    def _write_settings(self):
        """Save settings"""

        # Create directory if it does not exists
        if not os.path.exists('unblur_GUI_output'):
            os.mkdir('unblur_GUI_output')

        # Output file name
        strSaveName = str(
            QtGui.QInputDialog.getText(
                self,
                'Save settings',
                'File name:'
                )[0]
            )

        # Save data to file
        if strSaveName != '':
            with open(self.outputDir + '/' + strSaveName, 'w') as f:
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
                f.write('# End\n')

    def _load_settings(self):

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
        if strLoadName != '':

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

            # Load checked list
            end = separate_list[idxEnd] - separate_list[idxUnchecked] + 1
            self.listChecked = list(numpy.genfromtxt(
                strLoadName,
                skiprows=separate_list[idxChecked] + 1,
                skip_footer=end, comments='$',
                dtype=None
                ))
            if self.listChecked[0] == 'No':
                self.listChecked = []

            # Load unchecked list
            end = separate_list[idxEnd] - separate_list[idxCoordX] + 1
            self.listUnchecked = list(numpy.genfromtxt(
                strLoadName,
                skiprows=separate_list[idxUnchecked] + 1,
                skip_footer=end, comments='$',
                dtype=None
                ))
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

            # Refresh calculations
            self._refresh_calculations(goon=True)

            for row in arrThresh:
                listElement = []
                for index in xrange(2, len(row)):
                    listElement.append(row[index])
                if row[0] != self.modeOverall:
                    self.dictThresh.update({'{0} {1}'.format(row[0], row[1]): listElement})
                else:
                    self.dictThresh.update({'{0}'.format(row[0]): listElement})

            self.leStartGeneral.setText(arrGeneral[0])
            self.leStopGeneral.setText(arrGeneral[1])
            self._fill_widgets(mode=self.modeFrame)
            self._fill_widgets(mode=self.modeAngle)
            self._fill_widgets(mode=self.modeOverall)
            self._save_selection()
            self._invert_selection()
            self._invert_selection()

            # Enable everything and color all black
            self._default_color()
            self._enable_all()

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
                self.dictVisible = {}
                self.dictVisible.update({self.msPlotDrift: self.msPlotDrift.isVisible()})
                self.dictVisible.update({self.msPlotFrame: self.msPlotFrame.isVisible()})
                self.dictVisible.update({self.msPlotAngle: self.msPlotAngle.isVisible()})
                self.dictVisible.update({self.msAllPlotFrameAvg: self.msAllPlotFrameAvg.isVisible()})
                self.dictVisible.update({self.msAllPlotDrift: self.msAllPlotDrift.isVisible()})
                self.dictVisible.update({self.msAllPlotFrame: self.msAllPlotFrame.isVisible()})
                self.dictVisible.update({self.msAllPlotAngle: self.msAllPlotAngle.isVisible()})

                for key in self.dictVisible:
                    if self.dictVisible[key]:
                        key.hide()

                self.minimized = True
            elif self.minimized:
                self.minimized = False
                for key in self.dictVisible:
                    if self.dictVisible[key]:
                        key.show()
                        key.activateWindow()

    def closeEvent(self, event):
        """Change the closeEvent to close the application cleanly"""

        # Ignore the incomming event when pressing the "X" and
        # quit the application instead
        event.ignore()
        print('Bye Bye!')
        QtCore.QCoreApplication.instance().quit()


def _main():
    app = QtGui.QApplication(sys.argv)

    if len(sys.argv) > 2:
        print(
            'Too many arguments!\n'
            'Please start the GUI without argument or use \n'
            'quote marks around the wildcard "example-*_shift.txt"'
            )
        return None
    elif len(sys.argv) == 2:
        arg = sys.argv[-1]
    else:
        arg = None

    msApp = SXDriftUnblur(arg)
    msApp.show()

    app.exec_()


if __name__ == '__main__':
    _main()
