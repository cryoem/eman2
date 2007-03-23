import sys
from topology_hunter import TopologyHunter
from PyQt4 import QtCore, QtGui


class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)

        self.modules = []

        self.createActions()
        self.createMenus()
        
        self.statusBar().showMessage(self.tr("Protein Visualization Suite"))
        self.setWindowTitle(self.tr("Gorgon Prototype"))

    def addModule(self, module):
        self.modules.append(module)
        
    def createActions(self):
        self.exitAct = QtGui.QAction(self.tr("E&xit"), self)
        self.exitAct.setShortcut(self.tr("Ctrl+Q"))
        self.exitAct.setStatusTip(self.tr("Exit the application"))
        self.connect(self.exitAct, QtCore.SIGNAL("triggered()"), 
                     QtGui.qApp.closeAllWindows)
        
    def createMenus(self):
        self.fileMenu = self.menuBar().addMenu(self.tr("&File"))
        self.openMenu = self.fileMenu.addMenu(self.tr("&Open"))
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        self.optionsMenu = self.menuBar().addMenu(self.tr("&Options"))

    def menuOpen(self):
        return self.openMenu

    def menuOptions(self):
        return self.optionsMenu


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.addModule(TopologyHunter(window))
    window.show()
    sys.exit(app.exec_())
    
