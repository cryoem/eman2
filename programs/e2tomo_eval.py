#!/usr/bin/env python
# Muyuan Chen 2018-04
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from EMAN2 import *
import os
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emplot2d import EMPlot2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
from eman2_gui.valslider import ValSlider,CheckBox,ValBox
from eman2_gui.emshape import EMShape
from eman2_gui.emapplication import EMApp
import subprocess


def main():
	
	usage="This is a GUI program that allows users inspect tomograms easily. Simply run without argument in a tomogram project directory."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ppid", type=int,help="", default=None)

	parser.add_header(name="orblock1", help='', title="Click launch to evaluate reconstructed tomograms", row=1, col=0, rowspan=1, colspan=2, mode="")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	if not os.path.isdir("tomograms"): os.mkdir("tomograms")
		#print("No tomograms found. You must perform at least one reconstruction or manually populate the 'tomograms' directory with at least one reconstruction.")
		#sys.exit()

	app=EMApp()
	gui=TomoEvalGUI(options)
	gui.show()
	gui.raise_()
	app.exec_()
	E2end(logid)
	
	
class TomoEvalGUI(QtGui.QWidget):

	
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
		self.gbl.setColumnStretch(2,1)
		self.gbl.setRowStretch(0,1)

		# Micrograph list
		self.imglst=QtGui.QTableWidget(1, 3, self)
		self.imglst.verticalHeader().hide()
		
		self.gbl.addWidget(self.imglst,0,0,11,1)
		self.imglst.setColumnWidth(0,50)
		self.imglst.setColumnWidth(1,200)
		self.imglst_srtby=0
		hdr=self.imglst.horizontalHeader()
		self.imglst.cellClicked[int, int].connect(self.selimg)
		hdr.sectionPressed[int].connect(self.sortlst)
		
		self.wg_thumbnail=EMImage2DWidget(parent=self)
		self.wg_thumbnail.set_scale(1)
		self.wg_thumbnail_width=self.size().width()/3*.9
		self.wg_thumbnail.resize(self.wg_thumbnail_width,self.wg_thumbnail_width)
		#print self.wg_thumbnail_width
		self.wg_thumbnail.setMinimumHeight(330)
		self.gbl.addWidget(self.wg_thumbnail, 0,1,3,2)
		
		self.bt_show2d=QtGui.QPushButton("Show2D")
		self.bt_show2d.setToolTip("Show 2D images")
		self.gbl.addWidget(self.bt_show2d, 4,1,1,2)
		
		self.bt_runboxer=QtGui.QPushButton("Boxer")
		self.bt_runboxer.setToolTip("Run spt_boxer")
		self.gbl.addWidget(self.bt_runboxer, 5,1)
		
		self.bt_refresh=QtGui.QPushButton("Refresh")
		self.bt_refresh.setToolTip("Refresh")
		self.gbl.addWidget(self.bt_refresh, 5,2)
		
		self.bt_showtlts=QtGui.QPushButton("ShowTilts")
		self.bt_showtlts.setToolTip("Show raw tilt series")
		self.gbl.addWidget(self.bt_showtlts, 6,1)
		
		self.bt_plottpm=QtGui.QPushButton("TiltParams")
		self.bt_plottpm.setToolTip("Plot tilt parameters")
		self.gbl.addWidget(self.bt_plottpm, 6,2)
		
		self.bt_plotloss=QtGui.QPushButton("PlotLoss")
		self.bt_plotloss.setToolTip("Plot alignment loss")
		self.gbl.addWidget(self.bt_plotloss, 7,1)
		
		
		self.bt_plotctf=QtGui.QPushButton("PlotCtf")
		self.bt_plotctf.setToolTip("Plot CTF estimation")
		self.gbl.addWidget(self.bt_plotctf, 7,2)
		

		self.bt_show2d.clicked[bool].connect(self.show2d)
		self.bt_runboxer.clicked[bool].connect(self.runboxer)
		self.bt_plotloss.clicked[bool].connect(self.plot_loss)
		self.bt_plottpm.clicked[bool].connect(self.plot_tltparams)
		self.bt_showtlts.clicked[bool].connect(self.show_tlts)
		self.bt_refresh.clicked[bool].connect(self.update_files)
		self.bt_plotctf.clicked[bool].connect(self.plot_ctf)
		
		self.wg_2dimage=EMImage2DWidget()
		self.wg_2dimage.setWindowTitle("Tomo2D")
		self.cur_data=None
		
		self.wg_tltimage=EMImage2DWidget()
		self.wg_tltimage.setWindowTitle("Tiltseries")
		self.wg_tltimage.set_scale(.2)
		self.cur_tlt=None
		
		self.setspanel=QtGui.QListWidget(self)
		self.gbl.addWidget(self.setspanel, 8,1,2,2)
		self.itemflags=	Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)
		
		self.wg_notes=QtGui.QTextEdit(self)
		self.wg_notes.setText("Comments:")
		#self.wg_notes.setStyleSheet("color: rgb(150, 150, 150);")
		self.gbl.addWidget(self.wg_notes, 10,1,1,2)
		
		self.setspanel.itemChanged[QtGui.QListWidgetItem].connect(self.clickset)
		self.wg_notes.textChanged.connect(self.noteupdate)
		
		self.wg_plot2d=EMPlot2DWidget()
		
		self.update_files()

	def update_files(self):
		self.imginfo=[]
		
		#### getting information from json files
		files=sorted([os.path.join(self.path,f) for f in os.listdir(self.path)])
		ptclcls={}
		for i,name in enumerate(files):
			info=info_name(base_name(name))
			if os.path.isfile(info):
				dic={}
				nbox=0
				bxcls={}
				js=js_open_dict(info)
				
				if "boxes_3d" in js:
					boxes=js["boxes_3d"]
					nbox=len(boxes)
					
				if "ali_loss" in js:
					dic["loss"]=np.array(js["ali_loss"])
				else:
					dic["loss"]=[]
				
				if "tlt_params" in js:
					dic["tlt_params"]=np.array(js["tlt_params"])
				else:
					dic["tlt_params"]=[]
				
				if "tlt_file" in js:
					dic["tltfile"]=str(js["tlt_file"])
				else:
					dic["tltfile"]=""
					
				if "notes" in js:
					dic["notes"]=str(js["notes"])
				else:
					dic["notes"]=""
					
				if "defocus" in js:
					dic["defocus"]=np.array(js["defocus"])
				else:
					dic["defocus"]=[]
					
				if "phase" in js:
					dic["phase"]=np.array(js["phase"])
				else:
					dic["phase"]=[]
					
				if nbox>0 and "class_list" in js:
					cls=js["class_list"]
					for k in list(cls.keys()):
						vname=str(cls[k]["name"])
						n=np.sum([b[-1]==int(k) for b in boxes])
						
						if n>0:
							bxcls[vname]=n
							if vname in ptclcls:
								ptclcls[vname][1]+=n
							else:
								ptclcls[vname]=[1,n]
				
				dic["basename"]= os.path.basename(name).split(".")[0] #base_name(name)
				dic["e2basename"] = base_name(name)
				dic["filename"]=name
				dic["nbox"]=nbox
				dic["boxcls"]=bxcls
				dic["id"]=len(self.imginfo)
				self.imginfo.append(dic)
				js=None
		
		self.ptclcls=ptclcls
		
		self.update_list()
		
		#### update particle type list
		self.setspanel.clear()
		for k in list(self.ptclcls.keys()):
			v=self.ptclcls[k]
			kname="    {}\t:  {}".format(k, v[1])
			item=QtGui.QListWidgetItem(kname)
			item.setFlags(self.itemflags)
			self.setspanel.addItem(item)
			item.setCheckState(Qt.Checked)
	
	def update_list(self):
		#### update file list
		self.imglst.clear()
		self.imglst.setRowCount(len(self.imginfo))
		self.imglst.setColumnCount(4)
		self.imglst.setHorizontalHeaderLabels(["ID", "file name", "#box", "loss"])
		for i,info in enumerate(self.imginfo):
			#### use Qt.EditRole so we can sort them as numbers instead of strings
			it=QtGui.QTableWidgetItem()
			it.setData(Qt.EditRole, int(info["id"]))
			self.imglst.setItem(i,0,it)
			self.imglst.setItem(i,1,QtGui.QTableWidgetItem(str(info["basename"])))
			nbox=0
			for kname in list(info["boxcls"].keys()):
				if self.ptclcls[kname][0]==1:
					nbox+=info["boxcls"][kname]
			it=QtGui.QTableWidgetItem()
			it.setData(Qt.EditRole, int(nbox))
			self.imglst.setItem(i,2, it)
			if len(info["loss"])==0:
				loss=-1
			else: 
				loss=np.round(np.mean(info["loss"]), 2)
			
			it=QtGui.QTableWidgetItem()
			it.setData(Qt.EditRole, float(loss))
			self.imglst.setItem(i,3, it)
		
	def get_id_info(self):
		#### utility function to get the info of current selected row.
		crow=self.imglst.currentRow()
		idx=self.imglst.item(crow, 0).text()
		info=self.imginfo[int(idx)]
		return idx, info
	
	def plot_loss(self):
		idx, info=self.get_id_info()
		self.wg_plot2d.set_data(info["loss"], info["e2basename"], replace=True)
		self.wg_plot2d.show()
	
	def plot_tltparams(self):
		idx, info=self.get_id_info()
		if len(info["tlt_params"])==0: return 
		tpm=info["tlt_params"].T
		tpm=np.vstack([np.arange(len(tpm[0])), tpm])
		self.wg_plot2d.set_data(tpm, info["e2basename"], replace=True, linetype=0,symtype=0)
		self.wg_plot2d.setAxes(info["e2basename"], 4, 2)
		self.wg_plot2d.show()
		
	def show2d(self):
		idx, info=self.get_id_info()
		print("Showing 2D for image {} : {}".format(int(idx), info["filename"]))
		
		self.cur_data=EMData(info["filename"])
		self.wg_2dimage.list_idx=int(old_div(self.cur_data["nz"],2))		
		self.wg_2dimage.set_data(self.cur_data)
		self.wg_2dimage.show()
	
	def show_tlts(self):
		idx, info=self.get_id_info()
		print("Showing tilt series for image {} : {}".format(int(idx), info["filename"]))
		
		self.cur_tlt=EMData(info["tltfile"])
		self.wg_tltimage.list_idx=int(old_div(self.cur_tlt["nz"],2))		
		self.wg_tltimage.set_data(self.cur_tlt)
		self.wg_tltimage.show()
	
	def plot_ctf(self):
		idx, info=self.get_id_info()
		if len(info["defocus"])>0:
			data=info["defocus"]
			ln=len(data)
			if len(info["tlt_params"])==ln:
				data=np.vstack([info["tlt_params"][:,3], data]) 
			data=np.vstack([np.arange(ln), data]) 
			if len(info["phase"])>0:
				data=np.vstack([data, info["phase"]])
			
			self.wg_plot2d.set_data(data, info["e2basename"], replace=True, linetype=0,symtype=0)
			self.wg_plot2d.setAxes(info["basename"], 1, 2)
			self.wg_plot2d.show()
			
	
	def runboxer(self):
		idx, info=self.get_id_info()
		#### not doing this via launch_childprocess so the boxer wont be killed when one kill the browser...
		subprocess.Popen(["e2spt_boxer22.py",info["filename"]] )

	def clickset(self, item):
		name=str(item.text())
		k=name.split(':')[0].strip()
		chk=int(item.checkState() == Qt.Checked)
		if self.ptclcls[k][0]!=chk:
			self.ptclcls[k][0]=chk
			self.update_list()
		
		
	def noteupdate(self):
		notes=self.wg_notes.toPlainText()
		idx, info=self.get_id_info()
		if notes==info["notes"]:
			return
		try:	
			info["notes"]=notes
			infoname=info_name(info["e2basename"])
			js=js_open_dict(infoname)
			js["notes"]=notes
			js=None
		except:
			print("cannot write notes...")
		
	def selimg(self, row, col):
		idx=self.imglst.item(row, 0).text()
		info=self.imginfo[int(idx)]
		hdr=EMData(info["filename"], 0,True)
		iz=old_div(hdr["nz"],2)
		e=EMData(info["filename"], 0, False, Region(0,0,iz, hdr["nx"], hdr["ny"],1))
		
		fac=float(hdr["nx"])/self.bt_show2d.width()*1.01
		e.process_inplace('math.fft.resample',{"n":fac})
		self.wg_thumbnail.set_data(e)
		self.wg_notes.setText(str(info["notes"]))
		
		self.bt_plotctf.setEnabled(len(info["defocus"])>0)
		self.bt_plotloss.setEnabled(len(info["loss"])>0)
		self.bt_plottpm.setEnabled(len(info["tlt_params"])>0)
		
		
	def sortlst(self,col):
		print("Sort by",self.imglst.horizontalHeaderItem(col).text())
		self.imglst_srtby=1-self.imglst_srtby
		self.imglst.sortItems(col, self.imglst_srtby)
		
	def closeEvent(self,event):
		print("Exiting")
		
		self.wg_2dimage.close()
		self.wg_plot2d.close()
		self.wg_thumbnail.close()
		self.wg_2dimage.close()	
		
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()

