#!/usr/bin/env python
# Muyuan Chen 2018-04
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from EMAN2 import *
from EMAN2_utils import natural_sort
import os
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
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
	
	
class TomoEvalGUI(QtWidgets.QWidget):

	
	def __init__(self, options):
		
		self.path="tomograms/"
		QtWidgets.QWidget.__init__(self,None)


		self.win_size=[1000,680]
		self.setMinimumSize(self.win_size[0], self.win_size[1])

		# This object is itself a widget we need to set up
		self.gbl = QtWidgets.QGridLayout(self)
		self.gbl.setContentsMargins(8, 8, 8, 8)
		self.gbl.setSpacing(6)
		self.gbl.setColumnStretch(0,4)
		self.gbl.setColumnStretch(1,1)
		self.gbl.setColumnStretch(2,1)
		self.gbl.setRowStretch(0,1)

		# Micrograph list
		self.imglst=QtWidgets.QTableWidget(1, 3, self)
		#self.imglst.verticalHeader().hide()
		
		self.gbl.addWidget(self.imglst,0,0,11,1)
		self.imglst.setColumnWidth(0,50)
		self.imglst.setColumnWidth(1,200)
		self.imglst_srtby=0
		hdr=self.imglst.horizontalHeader()
		self.imglst.cellClicked[int, int].connect(self.selimg)
		hdr.sectionPressed[int].connect(self.sortlst)
		
		self.wg_thumbnail=EMImage2DWidget(parent=self)
		self.wg_thumbnail.set_scale(1)
		self.wg_thumbnail_width=int(self.size().width()/3*.9)
		self.wg_thumbnail.resize(self.wg_thumbnail_width,self.wg_thumbnail_width)
		#print self.wg_thumbnail_width
		self.wg_thumbnail.setMinimumHeight(330)
		self.gbl.addWidget(self.wg_thumbnail, 0,1,3,2)
		
		self.bt_show2d=QtWidgets.QPushButton("Show2D")
		self.bt_show2d.setToolTip("Show 2D images")
		self.gbl.addWidget(self.bt_show2d, 4,1,1,2)
		
		self.bt_runboxer=QtWidgets.QPushButton("Boxer")
		self.bt_runboxer.setToolTip("Run spt_boxer")
		self.gbl.addWidget(self.bt_runboxer, 5,1)
		
		self.bt_refresh=QtWidgets.QPushButton("Refresh")
		self.bt_refresh.setToolTip("Refresh")
		self.gbl.addWidget(self.bt_refresh, 5,2)
		
		self.bt_showtlts=QtWidgets.QPushButton("ShowTilts")
		self.bt_showtlts.setToolTip("Show raw tilt series")
		self.gbl.addWidget(self.bt_showtlts, 6,1)
		
		self.bt_plottpm=QtWidgets.QPushButton("TiltParams")
		self.bt_plottpm.setToolTip("Plot tilt parameters")
		self.gbl.addWidget(self.bt_plottpm, 6,2)
		
		self.bt_plotloss=QtWidgets.QPushButton("PlotLoss")
		self.bt_plotloss.setToolTip("Plot alignment loss")
		self.gbl.addWidget(self.bt_plotloss, 7,1)
		
		
		self.bt_plotctf=QtWidgets.QPushButton("PlotCtf")
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
		
		self.setspanel=TomoListWidget(self)
		self.gbl.addWidget(self.setspanel, 8,1,2,2)
		
		
		self.wg_notes=QtWidgets.QLineEdit(self)
		self.wg_notes.setText("Comments:")
		#self.wg_notes.setStyleSheet("color: rgb(150, 150, 150);")
		self.gbl.addWidget(self.wg_notes, 10,1,1,2)
		
		#self.setspanel.itemClicked[QtWidgets.QListWidgetItem].connect(self.clickset)
		self.wg_notes.textChanged.connect(self.noteupdate)
		
		self.wg_plot2d=EMPlot2DWidget()
		self.wg_plot2d.show()
		self.wg_plot2d.hide()
		
		
		self.update_files()

	def update_files(self):
		self.imginfo=[]
		
		#### getting information from json files
		files=natural_sort([os.path.join(self.path,f) for f in os.listdir(self.path)])
		ptclcls={}
		infonames=[]
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
					
				if (nbox>0) and ("class_list" in js):
					cls=js["class_list"]
					for k in list(cls.keys()):
						vname=str(cls[k]["name"])
						n=np.sum([b[-1]==int(k) for b in boxes])
						if n>0:
							bxcls[vname]=n
							if (info in infonames):
								continue
							if vname in ptclcls:
								ptclcls[vname][1]+=n
							else:
								ptclcls[vname]=[1,n]
				if ("curves" in js) and len(js["curves"])>0:
					cv=dic["curves"]=np.array(js["curves"])
					cids=np.unique(cv[:,-1])
					
					for ci in cids:
						ctag="_curves_{:02d}".format(int(ci))
					
						if ctag not in ptclcls:
							ptclcls[ctag]=[1,0]
					
						ptclcls[ctag][1]+=np.sum(cv[:,-1]==ci)
						bxcls[ctag]=np.sum(cv[:,-1]==ci)
						
				else:
					dic["curves"]=[]
				dic["basename"]= os.path.basename(name).split(".")[0] #base_name(name)
				dic["e2basename"] = base_name(name)
				dic["filename"]=name
				dic["nbox"]=nbox
				dic["boxcls"]=bxcls
				dic["id"]=len(self.imginfo)
				self.imginfo.append(dic)
				infonames.append(info)
				js=None
				
		
		self.ptclcls=ptclcls
		
		self.update_list()
		
		#### update particle type list
		self.setspanel.update_list(self.ptclcls)
	
	def update_list(self):
		#### update file list
		
		self.imglst.clear()
		self.imglst.setRowCount(len(self.imginfo))
		self.imglst.setColumnCount(5)
		self.imglst.setHorizontalHeaderLabels(["ID", "file name", "#box", "loss", "defocus"])
		self.imglst.setColumnHidden(0, True)
		for i,info in enumerate(self.imginfo):
			#### use Qt.EditRole so we can sort them as numbers instead of strings
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, int(info["id"]))
			self.imglst.setItem(i,0,it)
			self.imglst.setItem(i,1,QtWidgets.QTableWidgetItem(str(info["basename"])))
			nbox=0
			for kname in list(info["boxcls"].keys()):
				if self.ptclcls[kname][0]==1:
					nbox+=info["boxcls"][kname]
			#nbox+=len(info["curves"])
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, int(nbox))
			self.imglst.setItem(i,2, it)
			if len(info["loss"])==0:
				loss=-1
			else: 
				loss=np.round(np.mean(info["loss"]), 2)
				
			
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, float(loss))
			self.imglst.setItem(i,3, it)
			
			if len(info["defocus"])==0:
				df=-1
			else: 
				df=np.round(np.mean(info["defocus"]), 1)
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, float(df))
			self.imglst.setItem(i,4, it)
			#print(str(info["basename"]), nbox, loss, df)
			
		self.imglst.setVerticalHeaderLabels([str(i) for i in range(len(self.imginfo))])
		
	def get_id_info(self):
		#### utility function to get the info of current selected row.
		crow=self.imglst.currentRow()
		if crow>=0:
			idx=self.imglst.item(crow, 0).text()
			info=self.imginfo[int(idx)]
			return idx, info
		else:
			return None, None
	
	def plot_loss(self):
		idx, info=self.get_id_info()
		if idx==None: return
		self.wg_plot2d.set_data(info["loss"], info["e2basename"], replace=True)
		self.wg_plot2d.show()
	
	def plot_tltparams(self):
		idx, info=self.get_id_info()
		if idx==None: return
		if len(info["tlt_params"])==0: return 
		tpm=info["tlt_params"].T
		tpm=np.vstack([np.arange(len(tpm[0])), tpm])
		self.wg_plot2d.set_data(tpm, info["e2basename"], replace=True, linetype=0,symtype=0)
		self.wg_plot2d.setAxes(info["e2basename"], 4, 2)
		self.wg_plot2d.show()
		
	def show2d(self):
		idx, info=self.get_id_info()
		if idx==None: return
		print("Showing 2D for image {} : {}".format(int(idx), info["filename"]))
		
		self.cur_data=EMData(info["filename"])
		self.wg_2dimage.list_idx=int(old_div(self.cur_data["nz"],2))
		self.wg_2dimage.setWindowTitle(info["filename"])
		self.wg_2dimage.set_data(self.cur_data)
		self.wg_2dimage.show()
	
	def show_tlts(self):
		idx, info=self.get_id_info()
		if idx==None: return
		print("Showing tilt series for image {} : {}".format(int(idx), info["filename"]))
		
		if EMUtil.get_image_count(info["tltfile"])==1:
			self.cur_tlt=EMData(info["tltfile"])
			self.wg_tltimage.list_idx=int(old_div(self.cur_tlt["nz"],2))
		else:
			self.cur_tlt=EMData.read_images(info["tltfile"])
			self.wg_tltimage.list_idx=int(len(self.cur_tlt)/2)

		self.wg_tltimage.set_data(self.cur_tlt)
		self.wg_tltimage.show()
	
	def plot_ctf(self):
		idx, info=self.get_id_info()
		if idx==None: return
		if len(info["defocus"])>0:
			data=info["defocus"]
			ln=len(data)
			if len(info["tlt_params"])==ln:
				data=np.vstack([info["tlt_params"][:,3], data]) 
			data=np.vstack([np.arange(ln), data]) 
			if len(info["phase"])>0:
				data=np.vstack([data, info["phase"]])
			
			self.wg_plot2d.set_data(data, info["e2basename"], replace=True, linetype=0,symtype=0)
			self.wg_plot2d.setAxes(info["e2basename"], 1, 2)
			self.wg_plot2d.show()
			
	
	def runboxer(self):
		idx, info=self.get_id_info()
		if idx==None: return
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		### do not use launch_childprocess so the gui wont be frozen when boxer is opened
		if modifiers == QtCore.Qt.ShiftModifier:
			subprocess.Popen("e2tomo_drawcurve.py {} --ppid {}".format(info["filename"], os.getpid()),shell=True)
		else:
			subprocess.Popen("e2spt_boxer.py {} --ppid {}".format(info["filename"], os.getpid()),shell=True)
		#launch_childprocess()

	#def clickset(self, item):
		#name=str(item.text())
		#print(name)
		#k=name.split(':')[0].strip()
		#chk=int(item.checkState() == Qt.Checked)
		#if self.ptclcls[k][0]!=chk:
			#self.ptclcls[k][0]=chk
			#self.update_list()
		
		
	def noteupdate(self):
		notes=self.wg_notes.text()
		idx, info=self.get_id_info()
		if idx==None: return
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
		self.setspanel.highlight(info["boxcls"])
		
		
	def sortlst(self,col):
		if col<0: return
		print("Sort by",self.imglst.horizontalHeaderItem(col).text())
		self.imglst_srtby=1-self.imglst_srtby
		self.imglst.sortItems(col, self.imglst_srtby)
		
	def closeEvent(self,event):
		print("Exiting")
		
		self.wg_2dimage.close()
		self.wg_plot2d.close()
		self.wg_thumbnail.close()
		self.wg_2dimage.close()	
		self.wg_tltimage.close()
		
class TomoListWidget(QtWidgets.QListWidget):
	def __init__(self, parent=None):
		super(TomoListWidget, self).__init__(parent)
		self.parent=parent
		self.itemflags=Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)
		self.itemlst=[]
		self.lblen=10
		#self.itemClicked.connect(self.on_item_clicked)
		
	def get_text(self, label, ni=0, curi=-1):
		if curi<0:
			txt="{i:<{l}}:{n:>6}".format(i=label, l=self.lblen, n=ni)
		else:
			txt="{i:<{l}}:{n:>6} : {c:>5}".format(i=label, l=self.lblen, n=ni, c=curi)
		return txt
		
	def update_list(self, lst):
		self.clear()
		if len(lst.keys())>0:
			self.lblen=max([len(l) for l in lst.keys()])+3
		else:
			self.lblen=0
		self.itemlst=[]
		for k in list(lst.keys()):
			ni=lst[k][1]
			self.itemlst.append([k, ni])
			#print(txt)
			txt=self.get_text(k, ni)
			item=QtWidgets.QListWidgetItem(txt)
			item.setFlags(self.itemflags)
			item.setCheckState(Qt.Checked)
			try: item.setFont(QtGui.QFont("Monospace"))
			except:pass
			self.addItem(item)
			#v=self.ptclcls[k]
			
			#self.setspanel.add_item_number(k, v[1])
		
	
	def highlight(self, dic):
		lst=dic.keys()
		for i,key in enumerate(self.itemlst):
			if key[0] in lst:
				self.item(i).setText(self.get_text(key[0], key[1], dic[key[0]]))
				self.item(i).setForeground(QtGui.QColor("blue"))
			else:
				self.item(i).setText(self.get_text(key[0], key[1]))
				self.item(i).setForeground(QtGui.QColor("black"))

		
	
	def mousePressEvent(self, event):
		self._mouse_button = event.buttons()
		item=self.itemAt(event.pos())
		self.setCurrentItem(item)
		
		#print(item.text(), event.button())
		
		if event.button()==1:
			if item==None:
				for i in range(self.count()):
					self.item(i).setCheckState(Qt.Checked)
			else:
				if item.checkState()==Qt.Checked:
					item.setCheckState(Qt.Unchecked)
				else:
					item.setCheckState(Qt.Checked)
			#return
		
		elif event.button()==2:
			for i in range(self.count()):
				self.item(i).setCheckState(Qt.Unchecked)
			if item!=None:
				item.setCheckState(Qt.Checked)
		
		
		for i,key in enumerate(self.itemlst):
			#txt=str(self.item(i).text().split(":")[0]).strip()
			txt=key[0]
			self.parent.ptclcls[txt][0]=(self.item(i).checkState()==Qt.Checked)
			
			
		self.parent.update_list()
		#print(self.currentItem(),event)
		#print(item.text(), self._mouse_button)
		#super(TomoListWidget, self).mousePressEvent(event)
		
		
	def mouseReleaseEvent(self, event):
		#print("aaa")
		return

	#def on_item_clicked(self, item):
		#print(item.text(), self._mouse_button)


	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()

