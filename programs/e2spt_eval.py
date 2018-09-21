#!/usr/bin/env python
# Muyuan Chen 2018-04
from __future__ import print_function
from __future__ import division
from EMAN2 import *
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from eman2_gui import embrowser
from eman2_gui.emapplication import EMApp
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
from eman2_gui.emplot2d import EMPlot2DWidget


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--path", type=str,help="path", default="")
	parser.add_header(name="orblock1", help='', title="Click launch to evaluate subtomogram averages", row=1, col=0, rowspan=1, colspan=2)
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
		#self.imglst.setColumnWidth(0,50)
		self.imglst.setColumnWidth(1,200)
		self.imglst_srtby=0
		hdr=self.imglst.horizontalHeader()
		self.imglst.cellClicked[int, int].connect(self.select_folder)
		hdr.sectionPressed[int].connect(self.sortlst)
		
		self.dp_folder=QtGui.QComboBox()
		self.dp_folder.setToolTip("Folder suffix")
		self.gbl.addWidget(self.dp_folder, 0,1,1,1)
		sfxlst=["spt", "sptsgd"]
		for i in sfxlst:
			self.dp_folder.addItem(i)
		self.dp_folder.currentIndexChanged[int].connect(self.set_sfx)

		self.wg_thumbnail=EMScene3D()#parent=self)
		#self.wg_thumbnail.set_scale(1)
		self.wg_thumbnail_width=old_div(self.size().width(),4)*.9
		self.wg_thumbnail.resize(self.wg_thumbnail_width,self.wg_thumbnail_width)
		#print self.wg_thumbnail_width
		self.wg_thumbnail.setMinimumHeight(330)
		self.wg_thumbnail.setMinimumWidth(330)
		self.gbl.addWidget(self.wg_thumbnail, 1,1,2,2)

		self.wg_thumbnail.newnode=None

		self.bt_showbs=QtGui.QPushButton("ShowBrowser")
		self.bt_showbs.setToolTip("Show Browser")
		self.gbl.addWidget(self.bt_showbs, 2,1,1,2)
		self.bt_showbs.clicked[bool].connect(self.show_browser)

		self.bt_plotParms=QtGui.QPushButton("PlotParams")
		self.bt_plotParms.setToolTip("Examine particle orientations")
		self.gbl.addWidget(self.bt_plotParms, 3,1,1,2)
		self.bt_plotParms.clicked[bool].connect(self.plot_params)

		self.paramplot = EMPlot2DWidget()
		
		self.bt_plotFSC=QtGui.QPushButton("PlotFSCs")
		self.bt_plotFSC.setToolTip("Examine tightly masked FSCs from this SPT refinement")
		self.gbl.addWidget(self.bt_plotFSC, 4,1,1,2)
		self.bt_plotFSC.clicked[bool].connect(self.plot_fscs)

		self.fscplot = EMPlot2DWidget()

		self.setspanel=QtGui.QListWidget(self)
		self.gbl.addWidget(self.setspanel, 5,1,5,2)
		self.itemflags=	Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)
		
		#self.wg_notes=QtGui.QTextEdit(self)
		#self.gbl.addWidget(self.wg_notes, 10,1,1,2)
				
		self.setspanel.itemChanged[QtGui.QListWidgetItem].connect(self.click_set)
		#QtCore.QObject.connect(self.wg_notes,QtCore.SIGNAL("textChanged()"),self.noteupdate)
		
		#self.wg_plot2d=EMPlot2DWidget()
		
		self.update_files()
		self.browser=SPTBrowserWidget(withmodal=False,multiselect=False)

	def json_to_array(self,json):
		""" 
		Note: This is only meant to work with particle params json files. 
		It does NOT generalize in current state. 
		"""
		jd = js_open_dict(json)
		rows = []
		for k in jd.keys():
			dct = jd[k]
			tf = dct[u'xform.align3d']
			t = tf.get_trans()
			r = tf.get_rotation()
			row = [r["az"],r["alt"],r["phi"],t[0],t[1],t[2],dct["score"],dct["coverage"]]

			rows.append([float(x) for x in row])
		jd.close()
		return np.array(rows).T.tolist()

	def plot_params(self):
		for f in sorted(os.listdir(self.path)):
			if "particle_parms_" in f:
				parms = self.json_to_array("{}/{}".format(self.path,f))
				self.paramplot.set_data(parms,f)
		self.paramplot.show()
		return

	def plot_fscs(self):
		for f in sorted(os.listdir(self.path)):
			if "fsc_maskedtight_" in f and "ali.txt" in f:
				self.fscplot.set_data_from_file("{}/{}".format(self.path,f))
		self.fscplot.show()
		return

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
			for k in list(js.keys()):
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
				if not info.has_key(pm): continue
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
		self.browser.isup = True

	def select_folder(self, row, col):
		idx=int(self.imglst.item(row, 0).text())
		sfx=self.dp_folder.currentText()
		self.path="{}_{:02d}".format(sfx, idx)
		#if self.browser.isup == False: 
		#	self.show_browser()
		#	print("Showing {} in browser.".format(path))
		self.browser.setPath(self.path)

		vols = []
		for f in os.listdir(self.path):
			if "threed" in f:
				if "even" not in f and "odd" not in f:
					vols.append(f)#"{}".format(self.path,f))
		lastvol = sorted(vols)[-1]
		print("Displaying {} from {}".format(lastvol,self.path))
		volpath = "{}/{}".format(self.path,lastvol)
		name = os.path.basename(str(volpath))

		if self.wg_thumbnail.newnode:
			self.wg_thumbnail.newnode.setData(volpath)
		else:
			self.wg_thumbnail.newnode = EMDataItem3D(volpath)
			self.cur_3d_node=self.wg_thumbnail.newnode
			self.wg_thumbnail.insertNewNode(name, self.wg_thumbnail.newnode)
			self.wg_thumbnail.newnode.setTransform(self.wg_thumbnail.newnode.getParentMatrixProduct().inverse()*self.wg_thumbnail.newnode.getTransform())
			self.wg_thumbnail.isonode = EMIsosurface(self.wg_thumbnail.newnode, transform=Transform())
			self.wg_thumbnail.insertNewNode("Isosurface", self.wg_thumbnail.isonode, parentnode=self.wg_thumbnail.newnode)
		self.wg_thumbnail.updateSG()

	def sortlst(self,col):
		print("Sort by",self.imglst.horizontalHeaderItem(col).text())
		self.imglst_srtby=1-self.imglst_srtby
		self.imglst.sortItems(col, self.imglst_srtby)
		
	def closeEvent(self,event):
		print("Exiting")
		self.browser.close()
		if self.wg_thumbnail.main_3d_inspector != None:
			self.wg_thumbnail.main_3d_inspector.close()#closeEvent()
		self.paramplot.close()
		self.fscplot.close()

class SPTBrowserWidget(embrowser.EMBrowserWidget):

	isup = False

	def closeEvent(self, event) :
		E2saveappwin("e2display", "main", self)
		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d+self.viewhist :
			w.close()

		if self.infowin != None :
			self.infowin.close()

		self.isup = False

		event.accept()
		self.updtimer.stop()
		# self.app().close_specific(self)
		self.module_closed.emit()

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	
