#!/usr/bin/env python
# Muyuan Chen 2020-05
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import range
from EMAN2 import *
import numpy as np
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from EMAN2_utils import *


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--jsout", type=str,help="json output", default=None)
	parser.add_argument("--inplace", action="store_true", default=False ,help="overwrite input.")
	parser.add_argument("--invert", action="store_true", default=False ,help="invert direction.")
	parser.add_argument("--mode", type=str,help="choose from rad/line, default is line", default="line")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	app = EMApp()
	win=EMSptEval(app,args[0], options)
	win.show()
	app.execute()
	
	E2end(logid)
	


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None):
        fig = Figure(figsize=(4, 4))
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        
        a=np.random.randn(1000,2)
        self.axes.plot(a[:,0], a[:,1],'.')


class EMSptEval(QtWidgets.QMainWindow):

	def __init__(self,application,jsname, options):
		QtWidgets.QWidget.__init__(self)
		self.jsname=jsname
		self.options=options
		self.load_json()
		
		self.setMinimumSize(700,200)
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		self.gbl.setColumnStretch( 2, 10 ) 
		
		self.imglst=QtWidgets.QTableWidget(1, 1, self)
		self.imglst.setMinimumSize(450, 100)
		self.gbl.addWidget(self.imglst, 0,1,10,1)
		self.imglst.cellClicked[int, int].connect(self.on_list_selected)
		
		self.bt_save=QtWidgets.QPushButton("Save")
		self.gbl.addWidget(self.bt_save, 0,0,1,1)
		self.bt_save.clicked[bool].connect(self.save_json)

		self.update_list()
		self.plotwiny=MplCanvas(self)
		self.plotwiny.setMinimumSize(500, 500)
		self.ploty=self.plotwiny.axes
		self.gbl.addWidget(self.plotwiny, 0,2,7,1)
		self.plotwiny.mpl_connect('button_press_event', self.onclick_ploty)
		
		self.plotwinz=MplCanvas(self)
		self.plotwiny.setMinimumSize(500, 200)
		self.plotz=self.plotwinz.axes
		self.gbl.addWidget(self.plotwinz, 7,2,3,1)
		self.plotwinz.mpl_connect('button_press_event', self.onclick_plotz)
		self.cursel=-1
		
		
	def load_json(self):
		
		self.jsdict=js=dict(js_open_dict(self.jsname)).copy()
		keys=js.keys()
		self.data=[]
		print("Reading particles from json...")
		for k in keys:
			src, ii = eval(k)
			e=EMData(src, ii, True)
			srcname=e["file_twod"]#e["data_source"]
			self.data.append({"k":k,
					"srcfull":srcname,
					"src":base_name(srcname),
					"srci":e["data_n"],
					"pos":e["ptcl_source_coord"],
					"xf":js[k]["xform.align3d"], 
					"score":js[k]["score"],
					"toinv":False,
					})
			sys.stdout.write("\r  {}/{}        ".format(len(self.data), len(keys)))
			sys.stdout.flush()
		print()
		
		fs=natural_sort([d["src"] for d in self.data])
		files, count=np.unique(fs, return_counts=True)
		self.filenames=files
		self.fileinfo={f:{"id":i, "point":None, "nptcl":count[i]} for i,f in enumerate(files)}
		for fname in files:
			js=js_open_dict(info_name(fname))
			if js.has_key("global_ptcl_point"):
				self.fileinfo[fname]["point"]=js["global_ptcl_point"]
			else:
				self.fileinfo[fname]["point"]=[0,0,0]
			js.close()
			
		print("load {} particles from {} tomograms".format(len(self.data), len(files)))
		
	def save_json(self):
			
		if self.options.inplace:
			outname=self.jsname
		elif self.options.jsout:
			outname=self.options.jsout
		else:
			outname=self.jsname.replace(".json", "_mod.json")
			
		print("Writting output to {}".format(outname))
		rot=Transform({"type":"eman","alt":180, "az":0})
		jsout={}
		suminv=0
		for d in self.data:
			k=d["k"]
			
			val=self.jsdict[k].copy()
			if d["toinv"]:
				suminv+=1
				xf=val["xform.align3d"]
				val["xform.align3d"]=rot*xf
				
			jsout[k]=val
			
		if os.path.isfile(outname):
			os.remove(outname)
		js=js_open_dict(outname)
		js.update(jsout)
		js.close()
		print("{} particles inverted".format(suminv))
		
		for fname in self.filenames:
			v=self.fileinfo[fname]["point"]
			if type(v)==list:
				js=js_open_dict(info_name(fname))
				js["global_ptcl_point"]=v
				js.close()
		
	def onclick_ploty(self,event):
		if self.cursel>=0:
			fname=self.filenames[self.cursel]
			self.fileinfo[fname]["point"][0]=event.xdata
			self.fileinfo[fname]["point"][1]=event.ydata
			self.check_orient(fname, True)
		
	def onclick_plotz(self,event):
		if self.cursel>=0:
			fname=self.filenames[self.cursel]
			self.fileinfo[fname]["point"][0]=event.xdata
			self.fileinfo[fname]["point"][2]=event.ydata
			self.check_orient(fname, True)

		
	def on_list_selected(self, row, col):
		self.cursel=idx=int(self.imglst.item(row, 0).text())
		self.check_orient(self.filenames[idx], True)
		
		
	def check_orient(self, fname, doplot=True):
		data=[d for d in self.data if d["src"]==fname]
		
		pos=np.array([d["pos"] for d in data])
		ks=[d["k"] for d in data]
		xfs=[Transform(d["xf"]) for d in data]
		for x in xfs:
			x.set_trans(0,0,0)
			x.invert()
		vecs=np.array([xf.transform([0,0,1]) for xf in xfs])
		#vecs=np.array([xf.transform([-.44,-.89,0]) for xf in xfs])
		pc=np.mean(pos, axis=0)
		
		glbpt=self.fileinfo[fname]["point"]
		
		if self.options.mode=="rad":
			v=pos-np.array(glbpt)
			#v[:,:2]=0
			v=v/np.linalg.norm(v, axis=1)[:, None]
			dt=np.sum(vecs*v, axis=1)
		else:
			v=np.array(glbpt)-pc
			dt=np.dot(vecs, v/np.linalg.norm(v))
			
		if self.options.invert:
			dt*=-1
		inv=dt<0
		print("point: ({:.1f}, {:.1f}, {:.1f}), {:.1f} flipped".format(glbpt[0], glbpt[1], glbpt[2], np.mean(inv)*100))
		for i,d in enumerate(data):
			d["toinv"]=inv[i]
		
		if doplot:
			for j,plot in enumerate([self.ploty, self.plotz]):
				plot.cla()
				ii=inv
				plot.quiver(pos[ii,0], pos[ii,j+1], vecs[ii,0], vecs[ii,j+1], scale=50, color='r',width=.01)
				ii=np.logical_not(inv)
				plot.quiver(pos[ii,0], pos[ii,j+1], vecs[ii,0], vecs[ii,j+1], scale=50, color='b',width=.01)
				if self.options.mode=="rad":
					plot.plot(glbpt[0], glbpt[j+1],'.c', ms=10)
				else:
					plot.plot(pc[0], pc[j+1],'.c', ms=10)
					plot.plot([pc[0], glbpt[0]], [pc[j+1], glbpt[j+1]],'c')
					
			self.plotwiny.draw()
			self.plotwinz.draw()
		
	def update_list(self):
		#### update file list
		files=self.filenames
		self.imglst.clear()
		self.imglst.setRowCount(len(files))
		self.imglst.setColumnCount(3)
		self.imglst.setHorizontalHeaderLabels(["ID", "FileName", "Ptcls"])
		self.imglst.setColumnHidden(0, True)
		for i,fname in enumerate(files):
			#### use Qt.EditRole so we can sort them as numbers instead of strings
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, i)
			self.imglst.setItem(i,0,it)
			
			self.imglst.setItem(i,1,QtWidgets.QTableWidgetItem(fname))
			
			c=self.fileinfo[fname]["nptcl"]
			self.imglst.setItem(i,2,QtWidgets.QTableWidgetItem(str(c)))
			
			
		for i,w in enumerate([50,300]):
			self.imglst.setColumnWidth(i,w)
		#self.imglst.setVerticalHeaderLabels([str(i) for i in range(len(self.imginfo))])
	
	
if __name__ == '__main__':
	main()
	