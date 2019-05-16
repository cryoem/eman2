#!/usr/bin/env python
#
# Author: Orenbas, Turan 2019
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
#

import os
import EMAN2
import eman2_gui.emimage2d as emimage2d
import eman2_gui.emapplication as emapplication
import sp_fundamentals
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QVBoxLayout, QHBoxLayout
from PyQt5.QtWidgets import QPushButton, QListWidget, QListWidgetItem, QFileDialog
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtCore import Qt, QUrl, QFileInfo, pyqtSignal

input_file = EMAN2.EMData(700, 700)  # placeholder for EMImage2DWidget


class MainWindow(QMainWindow):
    def __init__(self, *args):
        QMainWindow.__init__(self, *args)
        self.current_image = None
        self.createComponents()
        self.loadPowerSpectrum(input_file)
        self.loadMicrographs(input_file)
        self.createLayout()
        # self.createMenu()
        self.createConnects()
        self.setWindowTitle('Micrograph Selection Tool')

    def createComponents(self):
        self.buttonLoad = QPushButton('Load')
        self.buttonSave = QPushButton('Save')
        self.buttonKeep = QPushButton('Keep (Right arrow key)')
        self.buttonKeep.setStyleSheet('QPushButton {color: green;}')
        self.buttonDiscard = QPushButton('Discard (Left arrow key)')
        self.buttonDiscard.setStyleSheet('QPushButton {color: red;}')
        self.fileList = QListWidget()
        self.micrograph = emimage2d.EMImage2DWidget(image=EMAN2.EMData(700, 700), application=app,
                                                    parent=self)
        self.micrograph.keyPressEvent = self.myKeyPressEvent
        self.powerSpectrum = emimage2d.EMImage2DWidget(image=EMAN2.EMData(700, 700),
                                                       application=app, parent=self)
        self.powerSpectrum.keyPressEvent = self.myKeyPressEvent
        self.powerSpectrum.set_FFT(2)

    def loadMicrographs(self, image_path):
        """creates an EMImage-Widget"""
        if self.current_image is None:
            self.current_image = EMAN2.EMData(image_path)
        self.micrograph.set_data(
            self.current_image)  # in project descriptopn as 'aaa' instead of 'micrograph'

    def loadPowerSpectrum(self, image_path):
        """creates power spectrum of a micrograph"""

        self.current_image = EMAN2.EMData(image_path)
        # self.current_image = sp_fundamentals.resample(self.current_image, 0.5)
        # image.process_inplace("math.realtofft")
        self.powerSpectrum.set_data(self.current_image)

        #
        # self.powerSpectrum.force_fft_redo()

    def createLayout(self):
        """layout: left, middle, right layouts, all three are vertically structured"""
        layoutMain = QHBoxLayout()

        layoutLeft = QVBoxLayout()
        layoutLeft.addWidget(self.buttonLoad)
        layoutLeft.addWidget(self.fileList)
        layoutLeft.addWidget(self.buttonSave)

        layoutMiddle = QVBoxLayout()
        layoutMiddle.addWidget(self.micrograph)
        layoutMiddle.addWidget(self.buttonKeep)

        layoutRight = QVBoxLayout()
        layoutRight.addWidget(self.powerSpectrum)
        layoutRight.addWidget(self.buttonDiscard)

        layoutMain.addLayout(layoutLeft, stretch=1)
        layoutMain.addLayout(layoutMiddle, stretch=2)
        layoutMain.addLayout(layoutRight, stretch=2)

        widgetMain = QWidget()
        widgetMain.setLayout(layoutMain)
        self.setCentralWidget(widgetMain)

    def createMenu(self):
        self.menuBar().addMenu("File")
        self.menuBar().addMenu("Help")

    def createConnects(self):
        self.buttonLoad.clicked.connect(buttonLoad_clicked)  # function call without ()-sign
        self.buttonSave.clicked.connect(buttonSave_clicked)
        self.buttonKeep.clicked.connect(buttonKeep_clicked)
        self.buttonDiscard.clicked.connect(buttonDiscard_clicked)
        self.fileList.clicked.connect(moveToImageOfCurrentRow)
        self.fileList.keyPressEvent = self.myKeyPressEvent

    # buttonSignal = pyqtSignal()

    def myKeyPressEvent(self, buttonSignal):
        if GUI.fileList.count() > 0:
            if buttonSignal.key() == Qt.Key_Right:
                print "Right"
                arrowKeyRight_clicked()
            elif buttonSignal.key() == Qt.Key_Left:
                print "Left"
                arrowKeyLeft_clicked()
            elif buttonSignal.key() == Qt.Key_Up:
                print "Up"
                arrowKeyUp_clicked()
            elif buttonSignal.key() == Qt.Key_Down:
                print "Down!"
                arrowKeyDown_clicked()


def buttonLoad_clicked():
    """"opens dialog for folder selection, """

    # dataLocation = QDesktopServices.storageLocation(QDesktopServices.HomeLocation)
    micrographs, _ = QFileDialog.getOpenFileNames(parent=None, caption=u"Open files",
                                                  filter=u"All files(*. *.mrc *.tiff *.tif *.hdf *.png)")  # name may be cofusing
    print "mics", micrographs
    global itemName  # maybe problematic?
    itemName = []
    for filePath in micrographs:
        """lists all selected micrographs (or other data)
        some changes may be good (either work work with filename or filePath??) (use consisting 
        variable names with capital letters for each word)"""
        # item = QListWidgetItem(filePath)
        # item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        # url = QUrl.fromLocalFile(filePath)
        print "Path", filePath
        filename = QFileInfo(filePath).fileName()
        item = QListWidgetItem(filename)
        item.setCheckState(Qt.Checked)
        GUI.fileList.addItem(item)

        filePath = str(filePath)
        itemName.append(filePath)
    if len(micrographs) > 0:
        GUI.fileList.setCurrentRow(0)
        GUI.loadPowerSpectrum(itemName[0])
        GUI.loadMicrographs(itemName[0])


def buttonSave_clicked():
    pathKeep, _ = QFileDialog.getSaveFileName(parent=None, caption='Save All Accepted Files as txt',
                                              filter="Text files (*.txt)", directory="keep.txt")

    try:
        fileKeep = open(pathKeep, 'w')
    except IOError:
        return

    pathDiscard, _ = QFileDialog.getSaveFileName(parent=None,
                                                 caption="Save All Rejected Files as txt",
                                                 filter="Text files (*.txt)",
                                                 directory="discard.txt")

    try:
        fileDiscard = open(pathDiscard, 'w')
    except IOError:
        return

    NumberOfItems = GUI.fileList.count()
    # selectedMicrographs = GUI.fileList.item
    for index in range(NumberOfItems):  # as often as number of files
        if (GUI.fileList.item(index).checkState()) == Qt.Checked:
            fileKeep.write(os.path.basename(itemName[index]) + "\n")
        else:
            fileDiscard.write(os.path.basename(itemName[index]) + "\n")
    fileKeep.close
    fileDiscard.close


def buttonKeep_clicked():
    if GUI.fileList.currentItem() is not None:
        setCheck()
        moveToNextItem()
        moveToImageOfCurrentRow()


def buttonDiscard_clicked():
    if GUI.fileList.currentItem() is not None:
        setUncheck()
        moveToNextItem()
        moveToImageOfCurrentRow()


def arrowKeyRight_clicked():
    setCheck()
    moveToNextItem()
    moveToImageOfCurrentRow()


def arrowKeyLeft_clicked():
    setUncheck()
    moveToNextItem()
    moveToImageOfCurrentRow()


def arrowKeyUp_clicked():
    moveToPrevItem()
    moveToImageOfCurrentRow()


def arrowKeyDown_clicked():
    moveToNextItem()
    moveToImageOfCurrentRow()


def setCheck():
    GUI.fileList.currentItem().setCheckState(Qt.Checked)


def setUncheck():
    GUI.fileList.currentItem().setCheckState(Qt.Unchecked)


def moveToNextItem():
    rowIndex = GUI.fileList.row(GUI.fileList.currentItem())
    GUI.fileList.setCurrentRow(rowIndex + 1)
    if rowIndex >= (GUI.fileList.count()) - 1:
        GUI.fileList.setCurrentRow(0)


def moveToPrevItem():
    rowIndex = GUI.fileList.row(GUI.fileList.currentItem())
    GUI.fileList.setCurrentRow(rowIndex - 1)
    if rowIndex >= (GUI.fileList.count()) - 1:
        GUI.fileList.setCurrentRow(0)


def moveToImageOfCurrentRow():
    rowIndex = GUI.fileList.row(GUI.fileList.currentItem())
    GUI.loadPowerSpectrum(itemName[rowIndex])
    GUI.loadMicrographs(itemName[rowIndex])

    print itemName[rowIndex]


if __name__ == '__main__':
    app = emapplication.EMApp()
    GUI = MainWindow()
    GUI.showMaximized()
    # GUI.show()
    app.exec_()
