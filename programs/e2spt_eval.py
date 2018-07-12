#!/usr/bin/env python
# Muyuan Chen 2018-04
from __future__ import print_function
from EMAN2 import *
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
#from eman2_gui.emimage2d import EMImage2DWidget
#from eman2_gui.emplot2d import EMPlot2DWidget
#from eman2_gui.emimagemx import EMImageMXWidget
#from eman2_gui.valslider import ValSlider,CheckBox,ValBox
from eman2_gui import embrowser
from eman2_gui.emapplication import EMApp
#import subprocess

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--path", type=str,help="path", default="")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	app=EMApp()
	gui=SptEvalGUI(options)
	gui.show()
	gui.raise_()
	app.exec_()
	E2end(logid)
	
	
class SptEvalGUI(QtGui.QWidget):

	
	def __init__(self, options):
		
		self.path="tomograms/"
		QtGui.QWidget.__init__(self,None)


		self.win_size=[1000,680]
		self.setMinimumSize(self.win_size[0], self.win_size[1])

		# This object is itself a widget we need to set up
		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setMargin(8)
		self.gbl.setSpacing(6)
		self.gbl.setColumnStretch(0,4)
		self.gbl.setColumnStretch(1,1)
		#self.gbl.setColumnStretch(2,1)
		self.gbl.setRowStretch(0,1)

		# Micrograph list
		self.imglst=QtGui.QTableWidget(1, 3, self)
		self.imglst.verticalHeader().hide()
		
		self.gbl.addWidget(self.imglst,0,0,10,1)
		self.imglst.setColumnWidth(0,50)
		self.imglst.setColumnWidth(1,200)
		self.imglst_srtby=0
		hdr=self.imglst.horizontalHeader()
		QtCore.QObject.connect(self.imglst,QtCore.SIGNAL("cellClicked(int,int)"),self.select_folder)
		QtCore.QObject.connect(hdr,QtCore.SIGNAL("sectionPressed(int)"),self.sortlst)
		
		self.dp_folder=QtGui.QComboBox()
		self.dp_folder.setToolTip("Folder suffix")
		self.gbl.addWidget(self.dp_folder, 0,1,1,1)
		sfxlst=["spt", "sptsgd"]
		for i in sfxlst:
			self.dp_folder.addItem(i)
		QtCore.QObject.connect(self.dp_folder,QtCore.SIGNAL("currentIndexChanged(int)"),self.set_sfx)
		
		
		self.bt_showbs=QtGui.QPushButton("ShowBrowser")
		self.bt_showbs.setToolTip("Show Browser")
		self.gbl.addWidget(self.bt_showbs, 1,1,1,1)
		self.bt_showbs.clicked[bool].connect(self.show_browser)
		
		self.setspanel=QtGui.QListWidget(self)
		self.gbl.addWidget(self.setspanel, 2,1,8,1)
		self.itemflags=	Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)
		
		#self.wg_notes=QtGui.QTextEdit(self)
		#self.gbl.addWidget(self.wg_notes, 10,1,1,2)
				
		QtCore.QObject.connect(self.setspanel,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.click_set)
		#QtCore.QObject.connect(self.wg_notes,QtCore.SIGNAL("textChanged()"),self.noteupdate)
		
		#self.wg_plot2d=EMPlot2DWidget()
		
		self.update_files()
		self.browser=embrowser.EMBrowserWidget(withmodal=False,multiselect=False)
		self.browser.show()
		
	def update_files(self):
		self.fldlst=[]
		self.paramlst={}
		self.initialized=False
		#foldersfx=["spt"]
		
		#for sfx in foldersfx:
		sfx=self.dp_folder.currentText()
		#print(sfx)
			
		flds=[f for f in os.listdir('.') if (os.path.isdir(f) and f.startswith("{}_".format(sfx)))]
		flds=sorted(flds)
		for f in flds:
			dic={"ID":int(f[len(sfx)+1:])}
			jsfile="{}/0_spt_params.json".format(f)
			if not os.path.isfile(jsfile):
				continue
			
			js=js_open_dict(jsfile)
			dic.update(js.data)
			for k in js.keys():
				if str(k) not in self.paramlst:
					self.paramlst[str(k)]=True
			
			self.fldlst.append(dic)
		
		self.setspanel.clear()
		for k in sorted(self.paramlst.keys()):
			item=QtGui.QListWidgetItem(k)
			item.setFlags(self.itemflags)
			self.setspanel.addItem(item)
			item.setCheckState(Qt.Checked)
		
		
		self.update_list()
		self.initialized=True
		return
	
	def update_list(self):
		
		#### update file list
		self.imglst.clear()
		self.imglst.setRowCount(len(self.fldlst))
		
		hdrs=["ID"]+sorted([i for i in self.paramlst if self.paramlst[i]])
		
		self.imglst.setColumnCount(len(hdrs))
		self.imglst.setHorizontalHeaderLabels(hdrs)
		
		for i,info in enumerate(self.fldlst):
			#### use Qt.EditRole so we can sort them as numbers instead of strings
			for j,pm in enumerate(hdrs):
				v=info[pm]
				try:
					v=float(v)
				except:
					v=str(v)
				it=QtGui.QTableWidgetItem()
				it.setData(Qt.EditRole, v)
				self.imglst.setItem(i,j,it)
				
	
	def click_set(self, item):
		if not self.initialized:
			return
		name=str(item.text())
		chk=item.checkState() == Qt.Checked
		if self.paramlst[name]!=chk:
			self.paramlst[name]=chk
			self.update_list()
		
	def set_sfx(self, i):
		if not self.initialized:
			return
		self.update_files()
		
		
	def show_browser(self):
		self.browser.show()
		
		
	def select_folder(self, row, col):
		idx=int(self.imglst.item(row, 0).text())
		sfx=self.dp_folder.currentText()
		path="{}_{:02d}".format(sfx, idx)
		print("Showing {} in browser.".format(path))
		self.browser.setPath(path)
		
		
	def sortlst(self,col):
		print("Sort by",self.imglst.horizontalHeaderItem(col).text())
		self.imglst_srtby=1-self.imglst_srtby
		self.imglst.sortItems(col, self.imglst_srtby)
		
	def closeEvent(self,event):
		print("Exiting")
		self.browser.close()
		
		#self.wg_2dimage.close()
		#self.wg_plot2d.close()
		#self.wg_thumbnail.close()
		#self.wg_2dimage.close()	
		
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	