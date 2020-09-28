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
####
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

import os
import EMAN2
import eman2_gui.emimage2d as emimage2d
import eman2_gui.emapplication as emapplication
from collections import deque
import time
import thread
from ..libpy import sp_fundamentals
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QVBoxLayout, QHBoxLayout
from PyQt5.QtWidgets import QPushButton, QListWidget, QListWidgetItem, QFileDialog, QMessageBox
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtCore import Qt, QUrl, QFileInfo, pyqtSignal
import eman2_gui.emapplication as emapplication



input_file = EMAN2.EMData(700, 700)  # placeholder for EMImage2DWidget
IMG_BUFFER_SIZE = 4
img_buffer = deque(maxlen=IMG_BUFFER_SIZE)
last_loaded_image = -1


class MainWindow(QMainWindow):

    def __init__(self,  *args):
        QMainWindow.__init__(self)
        self.current_image = None
        self.createComponents()
        self.loadFromBufferMemory()
        #print("Input file type", type(input_file))
        self.createMenu()
        self.createLayout()
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
        self.micrograph = emimage2d.EMImage2DWidget(image=EMAN2.EMData(700, 700), application=app, parent=self)
        self.powerSpectrum = emimage2d.EMImage2DWidget(image=EMAN2.EMData(700, 700), application=app, parent=self)
        self.all_items = []



    def loadFromBufferMemory(self):
        """loads the content from the buffer database"""
        image_list = []
        buffer_lines = []
        self.buffer_path = "/home/turi/Schreibtisch/buffermemory.txt"

        with open(self.buffer_path) as buffer:
            buffer = buffer.read()
            if len(buffer) > 0:
                buffer_lines = buffer.split(',')
                print("buffer", buffer_lines)
            else:
                print("Buffer is empty")


        for image_index in range(len(buffer_lines)):
            # cleans the image name and status to add them in the GUI
            image_path = str(buffer_lines[image_index].split(':')[0])
            print("imagepath", image_path)
            image_status = str(buffer_lines[image_index].split(':')[1])
            print("imagestatus", image_status)

            # remove symbols from the image_name string
            if "(u'" in image_path:
                image_path = image_path.replace("(u'", "")
            if "',)" in image_path:
                image_path = image_path.replace("',)", "")
            if "\n" in image_path:
                image_path = image_path.replace("\n", "")
            if "'" in image_path:
                image_path = image_path.replace("'", "")
            if "\"" in image_path:
                image_path = image_path.replace("\"", "")
            if "[" in image_path:
                image_path = image_path.replace("[", "")
            if "]" in image_path:
                image_path = image_path.replace("]", "")
            if "\\" in image_path:
                image_path = image_path.replace("\\", "")

            image_name = QFileInfo(image_path).fileName()
            image_list.append(image_path)     # adds the name
            item = QListWidgetItem(image_name)

            # remove symbols from the status string
            if "(u'" in image_status:
                image_status = image_status.replace("(u'", "")
            if "',)" in image_status:
                image_status = image_status.replace("',)", "")
            if "\n" in image_status:
                image_status = image_status.replace("\n", "")
            if "'" in image_status:
                image_status = image_status.replace("'", "")
            if "[" in image_status:
                image_status = image_status.replace("[", "")
            if "]" in image_status:
                image_status = image_status.replace("]", "")

            if image_status == "checked":
                item.setCheckState(Qt.Checked)
            if image_status == "unchecked":
                item.setCheckState(Qt.Unchecked)

            self.fileList.addItem(item)

        if len(image_list) > 0:
            self.uploadMicrographs(image_list=image_list, upload_from_buffer=True)
        else:
            pass



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
        self.actionOpenData = QAction(("Open file.."), self)
        self.actionSaveData = QAction(("Save"), self)
        self.actionSelectAll = QAction(("Select all"), self)
        self.actionDeselectAll = QAction(("Deselect all"), self)
        self.actionSetCheck_fromKeepFile = QAction("Upload Keep/Discard Files", self)
        self.actionQuit = QAction(("Quit"), self)
        self.actionQuit.setMenuRole(QAction.QuitRole)

        menu_file = self.menuBar().addMenu("File")
        menu_file.addAction(self.actionOpenData)
        menu_file.addAction(self.actionSaveData)
        menu_file.addAction(self.actionSelectAll)
        menu_file.addAction(self.actionDeselectAll)
        menu_file.addAction(self.actionSetCheck_fromKeepFile)
        menu_file.addSeparator()
        menu_file.addAction(self.actionQuit)

        self.menuBar().addMenu("Help")



    def createConnects(self):
        self.buttonLoad.clicked.connect(self.buttonLoad_clicked)  # function call without ()-sign
        self.buttonSave.clicked.connect(self.buttonSave_clicked)
        self.buttonKeep.clicked.connect(self.buttonKeep_clicked)
        self.buttonDiscard.clicked.connect(self.buttonDiscard_clicked)
        self.fileList.clicked.connect(self.changeImageByMouse)
        self.fileList.keyPressEvent = self.myKeyPressEvent
        self.micrograph.keyPressEvent = self.myKeyPressEvent
        self.powerSpectrum.keyPressEvent = self.myKeyPressEvent
        self.actionSelectAll.triggered.connect(self.selectAll)
        self.actionDeselectAll.triggered.connect(self.deselectAll)
        self.actionSetCheck_fromKeepFile.triggered.connect(self.setCheck_fromKeepFile)



    def myKeyPressEvent(self, buttonSignal):
        if self.fileList.count() > 0:
            if buttonSignal.key() == Qt.Key_Right:
                print
                "Right"
                self.arrowKeyRight_clicked()
            elif buttonSignal.key() == Qt.Key_Left:
                print
                "Left"
                self.arrowKeyLeft_clicked()
            elif buttonSignal.key() == Qt.Key_Up:
                print
                "Up"
                self.arrowKeyUp_clicked()
            elif buttonSignal.key() == Qt.Key_Down:
                print
                "Down!"
                self.arrowKeyDown_clicked()



    def loadMicrographsFromItemList(self, image_index):
        """creates an EMImage-Widget"""
        if self.current_image is None:
            image_path = self.all_items[image_index]  #
            self.current_image = EMAN2.EMData(image_path)
        self.micrograph.set_data(self.current_image)  # in project description as 'aaa' instead of 'micrograph'



    def loadPowerSpectrumFromItemList(self, image_index):
        """creates power spectrum of a micrograph"""
        # Try to load image from buffer
        # If buffer is empty, wait the image appears
        while True:
            load_successfull = False
            try:
                preload_image_path, preload_image_index, img, fft_img = img_buffer.popleft()
                if preload_image_index == image_index:
                    load_successfull = True
            except IndexError as e:
                print("Index error", e)
                time.sleep(1)

            if load_successfull:
                break

        self.current_image = img
        self.powerSpectrum.set_data(fft_img)



    def preload_images(self, micrograph_list, img_buffer):
        """
        Preloads IMG_BUFFER_SIZE number of images into the memory.

        :param micrograph_list: list of all micrographs ######## mistake --> here it gets just a string not a list
        :param img_buffer: deque as buffer for images
        :return: None
        """
        print("preload_images micrograph_list", micrograph_list)
        offset = 0
        last_index = -1
        while True:
            index = self.fileList.row(self.fileList.currentItem())
            if last_index != -1:
                if index - last_index == 1:
                    offset = offset - 1
                elif index - last_index != 0:
                    offset = 0

            if len(img_buffer) < IMG_BUFFER_SIZE and (index + offset) < len(micrograph_list):
                start = time.time()
                print("in", index + offset)
                print("micrograph_list", micrograph_list)
                filename = str(micrograph_list[index + offset])

                print("filename", filename)
                image = EMAN2.EMData(filename)

                fft_img = image.do_fft()

                fft_img.set_value_at(0, 0, 0, 0)
                fft_img.set_value_at(1, 0, 0, 0)
                fft_img.process_inplace("xform.phaseorigin.tocorner")
                fft_img.process_inplace("xform.fourierorigin.tocenter")
                fft_img = fft_img.get_fft_amplitude()

                img_buffer.append((filename, index + offset, image, fft_img))
                end = time.time()
                print("Put new image:", filename, "Current index", index, "Image+offset:",
                      index + offset, "offset", offset, "time", end - start)
                offset = offset + 1
            else:
                time.sleep(1)

            last_index = index



    def buttonLoad_clicked(self):
        """"opens dialog for folder selection, """

        image_dir = QFileDialog.getExistingDirectory(parent=None, caption=u"Open directory")
        image_dir = str(image_dir)

        image_list = []
        for root, dirs, files in os.walk(image_dir):
            for file in files:
                if (".mrc" or ".tiff" or ".tif" or ".hdf" or ".png") in file:
                    image_list.append(image_dir + "/" + file)
        #if no micrograph loaded
        if len(image_list) == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setInformativeText("No images in the selected folder")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

        if len(image_list) > 0:
            self.uploadMicrographs(image_list=image_list, upload_from_buffer=False)
        else:
            pass



    def buttonSave_clicked(self):
        """the keep and discard image names are added to two independent text files one for keep one for discard"""

        save_dir = QFileDialog.getExistingDirectory(parent=None, caption="Directory of selected and discarded file with image names")
        path_selected = str(save_dir) + "/" + "selected_images.txt"
        path_discarded = str(save_dir) + "/" + "discarded_images.txt"

        try:
            file_selected = open(path_selected, "w+")
            file_discarded = open(path_discarded, "w+")
        except IOError:
            return

        number_of_items = self.fileList.count()
        for index in range(number_of_items):  # as often as number of files
            if (self.fileList.item(index).checkState()) == Qt.Checked:
                file_selected.write(os.path.basename(self.all_items[index]) + "\n")
            else:
                file_discarded.write(os.path.basename(self.all_items[index]) + "\n")
        file_selected.close()
        file_discarded.close()



    def uploadMicrographs(self, image_list, upload_from_buffer):
        """loads the micrograph into the GUI and saves their path and status in buffer memory"""
        #global all_items  # maybe problematic?
        #all_items = []
        self.all_items = []
        image_path_and_status = []

        if upload_from_buffer == True:
            for image_path in image_list:
                # loads the micrograph names and status into the GUI
                self.all_items.append(image_path)

        if upload_from_buffer == False:
            self.fileList.clear()
            # add a try except tree
            for image_path in image_list:
                print
                "Path", image_path
                image_name = QFileInfo(image_path).fileName()
                item = QListWidgetItem(image_name)
                image_path = str(image_path)  # contains the path and the name of the image
                self.all_items.append(image_path)

                image_path_and_status.append(image_path + ":" + "checked")

                item.setCheckState(Qt.Checked)
                self.fileList.addItem(item)

            # remove symbols from the list
            image_path_and_status = str(image_path_and_status)
            if "'" in image_path_and_status:
                image_path_and_status = image_path_and_status.replace("'", "")
            if "[" in image_path_and_status:
                image_path_and_status = image_path_and_status.replace("[", "")
            if "]" in image_path_and_status:
                image_path_and_status = image_path_and_status.replace("]", "")
            if " " in image_path_and_status:
                image_path_and_status = image_path_and_status.replace(" ", "")

            with open(self.buffer_path, "w") as buffer:
                buffer.write(image_path_and_status)

        if len(image_list) > 0:
            # loads the micrographs into the GUI
            self.fileList.setCurrentRow(0)
            thread.start_new_thread(self.preload_images, (image_list, img_buffer,))
            self.loadPowerSpectrumFromItemList(0)
            self.loadMicrographsFromItemList(0)
            last_loaded_image = 0



    def writeToBuffer_check(self):
        """the status of the image is written to the buffer"""
        row_index = self.fileList.row(self.fileList.currentItem())

        with open(self.buffer_path, "r") as buffer:
            buffer_lines = buffer.read().split(",")
            buffer_line = buffer_lines[row_index].split(":")
            print(buffer_line)
            image_path = buffer_line[0]
            image_status = buffer_line[1]
            # order of checking importance: checked is included in unchecked !!!
            if "unchecked" in image_status:
                image_status = image_status.replace("unchecked", "checked")
            elif "checked" in image_status:
                pass
            buffer_lines[row_index] = str(image_path + ":" + image_status)

        buffer_lines = str(buffer_lines)
        if "'" in buffer_lines:
            buffer_lines = buffer_lines.replace("'", "")
        if "\"" in buffer_lines:
            buffer_lines = buffer_lines.replace("\"", "")
        if "[" in buffer_lines:
            buffer_lines = buffer_lines.replace("[", "")
        if "]" in buffer_lines:
            buffer_lines = buffer_lines.replace("]", "")
        if " " in buffer_lines:
            buffer_lines = buffer_lines.replace(" ", "")
        if "\\" in buffer_lines:
            buffer_lines = buffer_lines.replace("\\", "")

        print("Bufferline", str(buffer_lines))
        with open(self.buffer_path, "w") as buffer:
            buffer.write(str(buffer_lines))



    def writeToBuffer_uncheck(self):
        """the status of the image is written to the buffer"""
        row_index = self.fileList.row(self.fileList.currentItem())

        with open(self.buffer_path, "r") as buffer:
            buffer_lines = buffer.read().split(",")
            buffer_line = buffer_lines[row_index].split(":")
            image_path = buffer_line[0]
            image_status = buffer_line[1]
            # order of checking importance: checked is included in unchecked !!!
            if "unchecked" in image_status:
                pass
            elif "checked" in image_status:
                image_status = image_status.replace("checked", "unchecked")
            buffer_lines[row_index] = str(image_path + ":" + image_status)
            print("Puffer", buffer_lines)

        buffer_lines = str(buffer_lines)
        if "'" in buffer_lines:
            buffer_lines = buffer_lines.replace("'", "")
        if "\"" in buffer_lines:
            buffer_lines = buffer_lines.replace("\"", "")
        if "[" in buffer_lines:
            buffer_lines = buffer_lines.replace("[", "")
        if "]" in buffer_lines:
            buffer_lines = buffer_lines.replace("]", "")
        if " " in buffer_lines:
            buffer_lines = buffer_lines.replace(" ", "")
        if "\\" in buffer_lines:
            buffer_lines = buffer_lines.replace("\\", "")

        print("Bufferline", str(buffer_lines))
        with open(self.buffer_path, "w") as buffer:
            buffer.write(str(buffer_lines))



    def setCheck_fromKeepFile(self):
        """images are automatically selected or discarded according to saved keep files"""
        file_dirs = QFileDialog.getExistingDirectory(parent=None, caption=u"Open files")

        file_list = []
        files_selected = []
        for root, dirs, files in os.walk(file_dirs):
            for file in files:
                if ("selected_images.txt") in file:
                    file_list.append(file_dirs + "/" + file)
                    files_selected.append(file_dirs + "/" + file)
                if ("discarded_images.txt") in file:
                    file_list.append(file_dirs + "/" + file)

        # if no micrograph loaded
        if len(file_list) == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setInformativeText("No textfile in the selected folder")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

        else:
            for index in range(len(files_selected)):

                with open(files_selected[index], "r") as file:
                    selected_images = file.readlines()

                if len(selected_images) > 0:

                    if self.fileList.currentItem() is not None:
                        self.fileList.setCurrentRow(0)

                        for uploaded_image in self.all_items:
                            # image contains name and path
                            # all_items has the same order like self.fileList
                            uploaded_image = str(uploaded_image)
                            uploaded_image = uploaded_image.split("/")[-1]    # removes the path, extracts the name

                            for index in range(len(selected_images)):
                                if "\n" in selected_images[index]:
                                    selected_images[index] = selected_images[index].replace("\n", "")

                            if uploaded_image in selected_images:
                                self.setCheck()
                                self.moveToNextItem()
                            else:
                                self.setUncheck()
                                self.moveToNextItem()

                    else:
                        pass



    def buttonKeep_clicked(self):
        if self.fileList.currentItem() is not None:
            self.setCheck()
            self.moveToNextItem()
            self.showImageOfCurrentRow()



    def buttonDiscard_clicked(self):
        if self.fileList.currentItem() is not None:
            self.setUncheck()
            self.moveToNextItem()
            self.showImageOfCurrentRow()



    def arrowKeyRight_clicked(self):
        self.setCheck()
        self.moveToNextItem()
        self.showImageOfCurrentRow()



    def arrowKeyLeft_clicked(self):
        self.setUncheck()
        self.moveToNextItem()
        self.showImageOfCurrentRow()



    def arrowKeyUp_clicked(self):
        self.moveToPrevItem()
        self.showImageOfCurrentRow()



    def arrowKeyDown_clicked(self):
        self.moveToNextItem()
        self.showImageOfCurrentRow()



    def setCheck(self):
        self.fileList.currentItem().setCheckState(Qt.Checked)
        self.writeToBuffer_check()



    def setUncheck(self):
        self.fileList.currentItem().setCheckState(Qt.Unchecked)
        self.writeToBuffer_uncheck()



    def selectAll(self):
        """all images get the checked status"""
        print("bbb")
        self.fileList.setCurrentRow(0)
        self.loadPowerSpectrumFromItemList(0)
        self.loadMicrographsFromItemList(0)
        print(len(self.all_items))

        for row_index in range(len(self.all_items)):
            print("bbbbb")
            self.setCheck()
            self.moveToNextItem()



    def deselectAll(self):
        """all images get the unchecked status"""
        print("ccc")
        self.fileList.setCurrentRow(0)
        self.loadPowerSpectrumFromItemList(0)
        self.loadMicrographsFromItemList(0)
        print(len(self.all_items))

        for row_index in range(len(self.all_items)):
            print("ccccc")
            self.setUncheck()
            self.moveToNextItem()



    def moveToNextItem(self):
        row_index = self.fileList.row(self.fileList.currentItem())
        self.fileList.setCurrentRow(row_index + 1)
        if row_index >= (self.fileList.count()) - 1:
            self.fileList.setCurrentRow(0)



    def moveToPrevItem(self):
        row_index = self.fileList.row(self.fileList.currentItem())
        self.fileList.setCurrentRow(row_index - 1)
        if row_index >= (self.fileList.count()) - 1:
            self.fileList.setCurrentRow(0)



    def changeImageByMouse(self):
        self.showImageOfCurrentRow()



    def showImageOfCurrentRow(self):
        """the image of current row is shown"""
        global last_loaded_image    # maybe problematic?
        row_index = self.fileList.row(self.fileList.currentItem())
        if last_loaded_image - row_index != 0:
            self.loadPowerSpectrumFromItemList(row_index)
            self.loadMicrographsFromItemList(row_index)
            last_loaded_image = row_index


def main():
    app = emapplication.EMApp()
    GUI = MainWindow()
    GUI.showMaximized()
    # GUI.show()
    app.exec_()


if __name__ == '__main__':
    main()
