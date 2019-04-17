#!/usr/bin/env python
# Muyuan Chen 2018-04
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *
from EMAN2_utils import *
import numpy as np
import weakref
import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets, QtCore, QtOpenGL
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
from eman2_gui.emshape import EMShape
import scipy.spatial.distance as scipydist

def main():

	usage="""This shows the alignment result from e2tomogram.py uing the information from a tomorecon_xx folder.
	
	Example: e2tomo_showali.py --path <tomorecon_xx> 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path to tomorecon_xx directory to examine fiducial error.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=2, mode="fiderr")
	#parser.add_argument("--iteration", type=int,help="Refinement iteration number", default=2,guitype="intbox",row=1, col=0, rowspan=1, colspan=1,mode="fiderr")
	parser.add_argument("--ppid", type=int,help="ppid", default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	app = EMApp()

	drawer=EMDrawWindow(app,options)


	drawer.show()
	app.execute()
	E2end(logid)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

def get_circle(p,r):
	t=np.arange(0, np.pi*2, old_div(np.pi,10))
	pts=r*np.vstack([np.cos(t), np.sin(t)]).T
	pts+=p
	return pts

class Boxes(EMShape):
	def __init__(self, img, pks2d, dirs):
		self.isanimated = False
		self.shape=["boxes",0]
		self.image=img
		self.pks2d=pks2d
		#self.xfs=xfs
		self.dirs=dirs
		self.selid=-1
		#self.triangles=[]



	def draw(self,d2s=None,col=None):
		mi=self.image.list_idx
		#xf=self.xfs[mi]
		#print(mi,xf)
		dx=self.image.data["nx"]//2
		dy=self.image.data["ny"]//2
		ps=[]
		glColor3f( 1., 1., 1. );
		glPointSize(3)
		for pid, ptx in enumerate(self.pks2d[mi]):
			ps.append(ptx.tolist())
			a=np.array(ptx)-self.dirs[mi, pid]
			ps.append(a.tolist())

			lns=get_circle(ptx, 32)
			
			if pid==self.selid:
				glColor3f( .2, 1, 1 );
			else:
				glColor3f( .2, .2, 1 );
				
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(lns)
			glDrawArrays(GL_LINES, 0, len(lns))


		glColor3f( .2, .2, 1 );
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointerf(ps)
		glDrawArrays(GL_LINES, 0, len(ps))
		glPointSize(5)
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointerf([ps[i] for i in range(0,len(ps),2)])
		glDrawArrays(GL_POINTS, 0, old_div(len(ps),2))

		return


class EMDrawWindow(QtWidgets.QMainWindow):

	def __init__(self,application,options):
		
		self.options=options
		self.check_path(options.path)
		self.get_data(0)

		QtWidgets.QWidget.__init__(self)
		self.imgview = EMImage2DWidget()
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		
		self.lb_name=QtWidgets.QLabel(self.tomoname)
		self.lb_name.setWordWrap(True)
		self.gbl.addWidget(self.lb_name, 0,0,1,2)
		
		self.iterlst=QtWidgets.QListWidget()
		self.iterlst.itemflags=Qt.ItemFlags(Qt.ItemIsSelectable)
		
		for i in sorted(self.losses.keys()):
			txt="{:d}  :  loss = {:.1f}".format(i, self.losses[i])
			item=QtWidgets.QListWidgetItem(txt)
			self.iterlst.addItem(item)
			
		
		self.iterlst.currentRowChanged[int].connect(self.update_list)
		self.gbl.addWidget(self.iterlst,1,0,1,2)
		

		
		self.app=weakref.ref(application)


		self.imgview = EMImage2DWidget()
		self.boxes=Boxes(self.imgview, self.pks2d, self.dirs)
		self.shape_index = 0
		
		self.imgview.set_data(self.datafile)
		self.imgview.shapes = {0:self.boxes}
		self.imgview.show()
		self.imgview.mouseup.connect(self.on_mouseup)
		
		self.boxesviewer=EMImageMXWidget()
		self.boxesviewer.show()
		self.boxesviewer.set_mouse_mode("App")
		self.boxesviewer.setWindowTitle("Landmarks")
		self.boxesviewer.rzonce=True


		#glEnable(GL_POINT_SMOOTH)
		#glEnable( GL_LINE_SMOOTH );
		#glEnable( GL_POLYGON_SMOOTH );
		#glEnable(GL_BLEND);
		#glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
	def check_path(self, path):
		js=js_open_dict(os.path.join(path, "0_tomorecon_params.json"))
		self.tomoname=base_name(js["inputname"])
		js.close()
		
		self.losses={}
		for itr in range(5):
			fname=os.path.join(self.options.path, "loss_{:02d}.txt".format(itr))
			#print(fname)
			if os.path.isfile(fname):
				l=np.loadtxt(fname)
				self.losses[itr]=np.mean(l[:,1])
				
		#print(self.losses)
				
		
	def update_list(self, row):
		
		idx=sorted(self.losses.keys())[row]
		if idx==self.cur_iter:
			return
		print("Showing iteration {:d}...".format(sorted(self.losses.keys())[row]))
		
		self.get_data(idx)
		self.boxes.pks2d=self.pks2d
		self.boxes.dirs=self.dirs
		self.boxes.selid=-1
		self.shape_index = 0
		
		self.imgview.set_data(self.datafile)
		self.imgview.shapes = {0:self.boxes}
		#self.imgview.show()
		self.imgview.shapechange=1
		self.imgview.updateGL()
		self.boxesviewer.set_data([])
		self.boxesviewer.update()
		
	def select_landmark(self, x, y):
		nid=self.imgview.list_idx
		pks=self.pks2d[nid]
		d=np.sqrt(np.mean((pks-[x,y])**2, axis=1))
		if np.min(d)>64:
			return
		pid=np.argmin(d)
		print("Selecting landmark #{}...".format(pid))
		self.boxes.selid=pid
		self.imgview.shapechange=1
		self.imgview.updateGL()
		
		aname=os.path.join(self.options.path, "ptclali_{:02d}.hdf".format(self.cur_iter))
		#self.boxesviewer.set_data(aname)
		self.ptcls=EMData.read_images(aname)
		#self.ptcls=[p for p in self.ptcls if p["pid"]==pid]
		self.boxesviewer.set_data([p for p in self.ptcls if p["pid"]==pid])
		self.boxesviewer.update()
		
		
	def on_mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		self.select_landmark(x,y)
		
		#if event.button()&Qt.LeftButton:

		#if event.modifiers()&Qt.ControlModifier:
			##### interpolate curve from previous slices
			#print('new contour')
			#self.contour.add_point([x, y], True)

		#elif event.modifiers()&Qt.ShiftModifier:
			##### remove point
			#pts=np.array(self.contour.points)
		
		
	def get_data(self, itr):
		
		
		fname=os.path.join(self.options.path, "ali_{:02d}.hdf".format(itr))
		pname=os.path.join(self.options.path, "landmarks_{:02d}.txt".format(itr))
		pks=np.loadtxt(pname)
		imgs = EMData.read_images(fname)
		img=EMData(imgs[0]["nx"],imgs[0]["ny"], len(imgs))
		xfs=[]
		for i,m in enumerate(imgs):
			img.insert_clip(m, (0,0,i))
			xfs.append(m["xform.projection"])


		aname=os.path.join(self.options.path, "ptclali_{:02d}.hdf".format(itr))
		n=EMUtil.get_image_count(aname)
		dirs=np.zeros((len(imgs),len(pks),2))
		for i in range(n):
			e=EMData(aname, i, True)
			s=e["score"]
			dirs[e["nid"], e["pid"]]=[s[0], s[1]]
		#print(dirs)
		self.dirs=dirs=dirs*e["apix_x"]/imgs[0]["apix_x"]
		pks2d=[]
		
		dx=img["nx"]//2
		dy=img["ny"]//2
		for nid, xf in enumerate(xfs):
			p2d=[]
			for pid, p in enumerate(pks):
				pt=[p[0], p[1], p[2]]
				ptx=xf.transform(pt)
				ptx=[i/2 for i in ptx]
				ptx=[ptx[0]+dx, ptx[1]+dy]
				p2d.append(ptx)
				
			pks2d.append(p2d)
		
		self.pks2d=np.array(pks2d)
		#print(self.pks2d.shape)
				
		
		self.xfs=xfs
		self.pks=pks
		self.datafile=img
		self.cur_iter=itr
		

	def closeEvent(self, event):
		self.imgview.close()
		self.boxesviewer.close()


if __name__ == '__main__':
	main()
