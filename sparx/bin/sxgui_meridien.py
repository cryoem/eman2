#!/usr/bin/env python
from __future__ import print_function
# sxgui_meridien for analysing meridien outputs.
# Author 2017  Thorsten Wagner (thorsten.wagner@mpi-dortmund.mpg.de)
# Copyright (C) 2017 Max planck institute for molecular physiology, Dortmund
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

from builtins import range
from builtins import object
from PyQt4.QtCore import QObject, pyqtSignal, pyqtSlot, QThread, QString, QThreadPool, QTimer
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QDialog, QGridLayout, QTreeWidget, QMessageBox, QFontMetrics
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import os
import os.path
import json
import sys
import numpy as np


class DriverFileReader(QObject):

    sig_readfsc = pyqtSignal(object)
    sig_readfolders = pyqtSignal(str)
    sig_sendfolders = pyqtSignal(object)
    sig_sendfsc = pyqtSignal(object,object,object,object)


    def __init__(self):
        super(DriverFileReader, self).__init__()
        self.sig_readfsc.connect(self.handle_read_fsc_triggered)
        self.sig_readfolders.connect(self.handle_read_refinement_folders)

    @pyqtSlot(str)
    def handle_read_fsc_triggered(self, paths):
        """
        Slot for reading fsc triggered
        :return: 1D array with fsc values
        """
        runnames = [os.path.join(*word.split("/")[-2:]) for word in paths]
        runnumbers = [word[-3:] for word in paths]
        pixelsizes = []
        nnxos = []
        for i,path in enumerate(paths):
            tracker_reader = TrackerFileReader()
            trackerPath = "{0}Tracker_{1}.json".format(path+"/", runnumbers[i])
            v = tracker_reader.read_pixelsize_and_nnxo(trackerPath)
            pixelsizes.append(v[0])
            nnxos.append(v[1])

        self.sig_sendfsc.emit(self.read_fsc_values(paths),runnames,pixelsizes,nnxos)



    def read_fsc_values(self,paths):
        """
        Reads the fsc value from the driver.txt file


        :param paths: Paths to run-folders (eg. /my/path/main001)
        :return: nd array with fsc values
        """
        dataOfRuns = []

        for path in paths:
            driver_path = "{0}/driver_{1}.txt".format(path,path[-3:])
            data = np.genfromtxt(driver_path, usecols=0)
            dataOfRuns.append(data)

        return dataOfRuns

    @pyqtSlot(str)
    def handle_read_refinement_folders(self,path):

        self.sig_sendfolders.emit(path)

    @staticmethod
    def read_refinement_folders(path):
        """
        Reads the folder structure of a refinement folder


        :param path: Path to the refinement folder which contains the main folders
        :return: list of string with mainXXX folder names
        """
        main_dict = []
        for dict in sorted(os.listdir(path)):
            combined_path = os.path.join(str(path),dict)


            if os.path.isdir(combined_path) and ("main" in dict) and ("000" not in dict):
                valid = False
                for files in os.listdir(combined_path):
                    if "Tracker" in files:
                        valid = True
                if valid:
                    main_dict.append(dict)
        return main_dict











class FSCPlot(QDialog):


    def __init__(self,parent=None):
        super(FSCPlot,self).__init__(parent)
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot(self, fscvalues, names, pixelsizes, boxsizes):
        """
        Plots the FSC curves for several runs
        :param fscvalues: nd-array with fsc values. Row corresponds to an run.
        :param names: list of names of each run. These names are snown as legend.
        :param pixelsizes: list of used pixelsize for each run
        :param boxsizes: list of the used boxsize for each run
        :return: none, just plotting..
        """
        if not fscvalues:
            return
        # create an axis
        ax = self.figure.add_subplot(111)

        for i,entry in enumerate(fscvalues):
            fsc_x = np.empty(len(entry), dtype=float)
            for fsc_value in range(len(entry)):
                fsc_x[fsc_value] = float(fsc_value) / float(pixelsizes[i] * boxsizes[i])

            ax.plot(fsc_x, entry,label=names[i])

        plt.grid()
        plt.title("Fourier shell correlation curves")
        plt.ylabel("FSC")
        plt.xlabel("Spatial frequency / 1/A")
        plt.legend()
        plt.gca().set_ylim(bottom=-0.05)
        self.canvas.draw()
        self.show()



class MonitorRefinementFolder(QObject):

    def __init__(self, path, sig_update_tree, parent=None):
        super(MonitorRefinementFolder, self).__init__(parent)
        self.refinement_folder = path
        self.reader = DriverFileReader()
        self.mainDicts = set(self.reader.read_refinement_folders(path))
        self.sig_update_tree = sig_update_tree

    @pyqtSlot()
    def update(self):
        '''
        Checks if there is a new main folder in a specified refinement folder and updates the QTree. If the
        refinement folder was deleted, the corresponding element in the QTree is deleted.
        :return:
        '''
        if os.path.exists(self.refinement_folder):
            current_dictionaries = set(self.reader.read_refinement_folders(self.refinement_folder))

            new_dictionaries = current_dictionaries.difference(self.mainDicts)
            if len(new_dictionaries) >= 1:
                list_new_dictionaries = list(new_dictionaries)

                for dictionary in list_new_dictionaries:
                    if os.path.isfile("{0}/{1}/driver_{2}.txt".format(self.refinement_folder, dictionary, dictionary[-3:])) and \
                            os.path.isfile("{0}/{1}/Tracker_{2}.json".format(self.refinement_folder, dictionary, dictionary[-3:])):
                        next_item = QtGui.QTreeWidgetItem([dictionary])
                        next_item.setCheckState(0, QtCore.Qt.Unchecked)
                        delete=False
                        self.sig_update_tree.emit(next_item, self.refinement_folder, delete)
                        self.mainDicts = current_dictionaries
        else:
            delete = True
            self.sig_update_tree.emit(None, self.refinement_folder, delete)
            self.setParent(None)


class ResolutionOverviewPlot(QDialog):
    def __init__(self, parent=None):
        super(ResolutionOverviewPlot, self).__init__(parent)
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot(self, resolutions_05, resolution_0143):
        """
        Plots the resolution over all runs for FSC 0.5 and FSC 0.143

        :param resolutions_05: list of spatial resolution with FSC 0.5
        :param resolution_0143: list of spatial resolution with FSC 0.143
        :param pixelsizes: list of used pixelsize for each run
        :param boxsizes: list of used boxsizes
        :return: none
        """
        if not resolutions_05 and not resolution_0143:
            return
        # create an axis
        ax = self.figure.add_subplot(111)
        xruns = list(range(1,len(resolution_0143)+1))

        ax.plot(xruns, resolution_0143, "x-",label="Resolution FSC 0.143")
        xruns = list(range(1, len(resolutions_05)+1))
        ax.plot(xruns, resolutions_05, "x-", label="Resolution FSC 0.5")

        plt.grid()
        plt.title("Resolution curves")
        plt.ylabel("Resolution [A]")
        plt.xlabel("Runs")
        plt.legend()
        self.canvas.draw()
        self.show()




class TrackerFileReader(object):

    def __init__(self):
        pass

    def read_pixelsize_and_nnxo(self, path):
        """

        :param path:  Path to run folder
        :return: list of pixel and boxsize
        """
        with open(path, 'r') as read:
            tracker = json.load(read)
        pixel_size = tracker['constants']['pixel_size']
        nnxo = tracker['constants']['nnxo']

        return [pixel_size, nnxo]

    def _read_res0143_and_res05(self, path):
        """
        :param path: Paths to run folder
        """
        with open(path, 'r') as read:
            tracker = json.load(read)
        res_05 = tracker['currentres']
        res_0143 = tracker['fsc143']


        return [res_0143, res_05]

    def read_res0143_and_res05(self, path):

        """
        :param path: Path to refinement folder
        """
        res_0143 = []
        res_05 = []

        main_paths = DriverFileReader.read_refinement_folders(path)

        for main_path in main_paths:
            runnumber = main_path[-3:]
            tracker_path = "{0}/{1}/Tracker_{2}.json".format(path, main_path, runnumber)
            if os.path.isfile(tracker_path):
                pixelsize,boxsize = self.read_pixelsize_and_nnxo(tracker_path)
                res0143, res05 = self._read_res0143_and_res05(tracker_path)
                res0143 = float(pixelsize * boxsize) / float(res0143)
                res05 = float(pixelsize * boxsize) / float(res05)
                res_0143.append(res0143)
                res_05.append(res05)

        return [res_0143, res_05]







class MainWindow(QtGui.QMainWindow):

    sig_update_tree = pyqtSignal(object, object, object)
    sig_show_overview_plot = pyqtSignal()

    def __init__(self, font, parent=None):
        super(MainWindow, self).__init__(parent)
        self.font = font
        self.setWindowTitle("Meridien")
        central_widget = QtGui.QWidget(self)
        self.setCentralWidget(central_widget)

        #Center on screen
        resolution = QtGui.QDesktopWidget().screenGeometry()
        self.move((resolution.width() / 2) - (self.frameSize().width() / 2),
                  (resolution.height() / 2) - (self.frameSize().height() / 2))

        """
        Setting up menu bar
        """

        close_action = QtGui.QAction('Close', self)
        close_action.setShortcut("Ctrl+Q")
        close_action.setStatusTip('Leave the app')
        close_action.triggered.connect(lambda: self.close())

        open_refinement_folder = QtGui.QAction('Open Refinement Folder', self)
        open_refinement_folder.triggered.connect(self.open_refinement_folder)

        self.mainMenu = self.menuBar()
        self.fileMenu = self.mainMenu.addMenu('&File')
        self.fileMenu.addAction(open_refinement_folder)
        self.fileMenu.addAction(close_action)
        self.refinement_folder = ""

        create_new_fsc_plot = QtGui.QAction('&New FSC plot', self)
        create_new_fsc_plot.triggered.connect(self.event_ontriggered_show_fsc_plot)

        create_new_overview_plot = QtGui.QAction('&New resolution overview plot', self)
        create_new_overview_plot.triggered.connect(self.event_show_resolution_overview_plot)
        self.plotMenu = self.mainMenu.addMenu('&Plot')
        self.plotMenu.addAction(create_new_fsc_plot)
        self.plotMenu.addAction(create_new_overview_plot)

        """
        Setup other components
        """
        self.layout = QGridLayout(central_widget)
        self.setMenuBar(self.mainMenu)

        self.tree = QTreeWidget(self)
        self.tree.setHeaderHidden(True)
        self.layout.addWidget(self.tree, 1, 0)

        self.root_items_path_dictionary = {}


        # Threads
        self.threadpool = QThreadPool()
        self.thread_list = []
        thr = QThread(self)
        thr.start()

        self.reader = DriverFileReader()
        self.reader.moveToThread(thr)
        self.thread_list.append(thr)
        self.timer = QTimer(self)

        # Connect signals
        self.reader.sig_sendfolders.connect(self.fill_tree)
        self.reader.sig_sendfsc.connect(self.show_dialog_fsc)
        self.tree.itemChanged.connect(self._event_select_deselect_all)
        self.sig_update_tree.connect(self.update_tree)
        self.sig_show_overview_plot.connect(self.event_show_resolution_overview_plot)
        self.show()
        self.open_refinement_folder()

        self.monitor = None

    @pyqtSlot(object,object,object)
    def update_tree(self, item, monitored_folder, delete):
        '''

        :param item: main folder to add to tree item belonging to folder montiored_folder
        :param monitored_folder: Refinement folder
        :param delete: if this is true, it will delete the monitored folder
        :return:
        '''
        if delete:
            root = self.root_items_path_dictionary[monitored_folder]
            index = self.tree.indexOfTopLevelItem(root)
            self.tree.takeTopLevelItem(index)
        else:
            root = self.root_items_path_dictionary[monitored_folder]
            root.addChild(item)
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Meridien - There is data of a new run: " + str(item.text(0)))
            msg.setWindowTitle("Update")
            msg.exec_()

    @pyqtSlot(object, object, object, object)
    def show_dialog_fsc(self, *args):
        """
        Shows fsc plot
        :param args: list of arguments
        :return: None
        """
        fsc_plot = FSCPlot(self)
        fsc_plot.plot(*args)

    @pyqtSlot()
    def event_show_resolution_overview_plot(self):
        """
        Trigger for fsc plot
        :return: None
        """
        tracker_reader = TrackerFileReader()
        res_0143, res_05 = tracker_reader.read_res0143_and_res05(self.refinement_folder)

        res_plot = ResolutionOverviewPlot(parent=self)
        res_plot.setWindowTitle(self.refinement_folder)
        res_plot.plot(res_05, res_0143)


    def event_ontriggered_show_fsc_plot(self):
        """
        Trigger for fsc plot
        :return: none
        """
        checked_runs = self._get_checked_items()
        paths = []
        for checked_run in checked_runs:

            refinement_folder = checked_run.parent().text(0)

            path = os.path.join(str(refinement_folder), str(checked_run.text(0)))

            paths.append(path)
        self.reader.sig_readfsc.emit(paths)

    def _event_select_deselect_all(self, root_tree_item):
        '''
        Checks all childs of root_tree_item
        :param root_tree_item: QTreeWidgetItem which represents a refinement folder
        :return:
        '''

        if root_tree_item in list(self.root_items_path_dictionary.values()):
            number_of_runs = root_tree_item.childCount()
            for i in range(number_of_runs):
                child_run = root_tree_item.child(i)
                child_run.setCheckState(0, root_tree_item.checkState(0))

    def _get_checked_items(self):
        """
        Finds the checked elements in tree
        :return: List of checked QItemWidget
        """
        checked_runs = []
        for root in list(self.root_items_path_dictionary.values()):
            number_of_runs = root.childCount()

            for i in range(number_of_runs):
                main_run = root.child(i)
                if main_run.checkState(0) == QtCore.Qt.Checked:
                    checked_runs.append(main_run)
        return checked_runs

    def open_refinement_folder(self):
        """
        Let the user choose the refinement folder and adds it to the RefinementFolder-Tree
        :return: none
        """
        self.refinement_folder = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Refinement Directory"))
        if self.refinement_folder == "":
            return
        if self.refinement_folder in self.root_items_path_dictionary:
            return
        self._open_refinement_folder(self.refinement_folder)
        self.sig_show_overview_plot.emit()

    def _open_refinement_folder(self, path):
        """
        Reads the refinement folder, setup the folder daemon and signals
        :param path: Path to refinement folder

        """

        if path != '':

            #for i in reversed(range(self.root.childCount())):
            #    self.root.removeChild(self.root.child(i))
            name = os.path.basename(path)
            qname = QString(name)
            root = QtGui.QTreeWidgetItem([str(path)])
            self.root_items_path_dictionary[str(path)] = root
            self.tree.addTopLevelItem(root)
            fm = QFontMetrics(self.font)
            w = fm.width(path)
            self.tree.setMinimumWidth(w+150)
            #self.root.setText(0, qname)
            self.reader.sig_readfolders.emit(path)
            self.monitor = MonitorRefinementFolder(path, self.sig_update_tree,self)
            self.timer = QTimer(self)
            self.timer.timeout.connect(self.monitor.update)
            self.timer.start(2000)

    def closeEvent(self, close_event):
        '''
        Closes all threads.
        '''
        for thr in self.thread_list:
            thr.quit()
            thr.wait()

    @pyqtSlot(object)
    def fill_tree(self, path_to_refinement_folder):
        '''
        Reads all runs in path_to_refinement_folder and add them as child to the corresponding root element in the tree
        :param path_to_refinement_folder: Path to refinement folder
        :return: none
        '''
        root = self.root_items_path_dictionary[str(path_to_refinement_folder)]
        root.setCheckState(0,QtCore.Qt.Unchecked)
        main_dicts = DriverFileReader.read_refinement_folders(path_to_refinement_folder)
        for dictionary in main_dicts:
            next_item = QtGui.QTreeWidgetItem([dictionary])
            next_item.setCheckState(0, QtCore.Qt.Unchecked)
            root.addChild(next_item)

    def close_application(self):
        '''
        Close the application
        :return: none
        '''
        sys.exit()


def run(args=None):
    app = QtGui.QApplication(sys.argv)

    gui = MainWindow(app.font())
    sys.exit(app.exec_())

if __name__ == '__main__':
    run()
