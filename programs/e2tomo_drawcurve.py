#!/usr/bin/env python
# Muyuan Chen 2017-03
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
from PyQt5 import QtGui, QtWidgets, QtCore, QtOpenGL
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emshape import EMShape
import scipy.spatial.distance as scipydist

def main():

	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomogram",help="Specify the tomogram to be segmented.", default="", guitype='filebox', browser="EMTomoTable(withmodal=True,multiselect=False)",  row=0, col=0,rowspan=1, colspan=2, mode="tomoseg")
	#parser.add_argument("--load", type=str,help="Load previous contour segmentation.", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)",  row=1, col=0,rowspan=1, colspan=2, mode="tomoseg")
	#parser.add_argument("--noupdate", action="store_true",help="do not erase shapes", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	

	app = EMApp()

	drawer=EMDrawWindow(app,options,datafile=args[0])


	drawer.show()
	app.execute()
	E2end(logid)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

class Contour(EMShape):
	def __init__(self, img, points):
		self.points=points
		self.isanimated = False
		self.shape=["scr_contour",0]
		self.image=img
		self.triangles=[]
		self.select=0
		self.classid=0
		
		
	
		

	def add_point(self, newpt=[], newcontour=False):



		zpos=self.image.list_idx
		if len(self.points)==0 and len(newpt)>0:
			self.points.append([newpt[0], newpt[1], zpos, self.select, self.classid])

		if newcontour==False:

			ci=self.select
			#### separate the points on the current contour and everything else
			nppts=np.array(self.points)
			sel=np.logical_and(nppts[:,3]==ci, nppts[:,4]==self.classid)
			pts=nppts[sel,:3].copy()
			otherpts=nppts[sel==False].copy()

			if len(pts)<3:
				if len(newpt)>0:
					self.points.append([newpt[0], newpt[1], zpos, ci, self.classid])
					return
			else:
				thr=1000.
				newpt=np.array([newpt[0], newpt[1], zpos])
			

			if len(newpt)>0:
				#### add a point
				dst=np.array(scipydist.cdist([newpt], pts)[0])
				rg=np.arange(len(dst), dtype=int)+1
				#print rg
				rg[-1]=rg[-2]
				d0=dst[0]
				d1=dst[-1]
				dst+=dst[rg]
				dst=np.append(d0*3, dst)
				dst[-1]=d1*3

				mi=np.argmin(dst)
				#print(mi, len(pts))
				#ci=pts[mi, 3]
				if mi==0:
					pts=np.vstack([ newpt, pts])
				if mi==len(pts):
					pts=np.vstack([pts, newpt])
				pts=np.insert(pts, mi+0, newpt, axis=0)


			allpts=[]
			#### a simple tsp solver...
			#idx=np.where(pts[:,3]==ci)[0]
			pp=pts.copy()
			path=np.arange(len(pp), dtype=int)

			dmat=scipydist.squareform(scipydist.pdist(pp))
			dmat+=np.eye(len(dmat))*thr
			if len(pp)>=3:
				niter=2000
			else:
				niter=0

			calc_dist=lambda path: np.sum([dmat[path[i], path[i+1]] for i in range(len(path)-1)])
			dst=calc_dist(path)
			nochange=0
			for k in range(niter):
				p0=path.copy()
				i0,i1= np.sort(np.random.randint(0,len(pp)+1, 2))
				if abs(i0-i1)%len(path)<2: continue
				path= np.hstack([path[:i0], path[i0:i1][::-1], path[i1:]])
				d=calc_dist(path)
				#print(i0, i1, d, dst)
				#print(path)
				if d>=dst:
					path=p0
					nochange+=1
					if nochange>200: break

				else:
					dst=d
					nochange=0


			allpts=[[pp[i][0], pp[i][1], pp[i][2], ci, self.classid] for i in path]

			self.points=otherpts.tolist()
			self.points.extend(allpts)
			#self.points=allpts#.tolist()

		else:
			#### start a new contour
			pts=np.array(self.points)
			ci=np.max(pts[:,3])+1
			self.select=ci
			self.points.append([newpt[0], newpt[1], zpos, ci, self.classid])




	def draw(self,d2s=None,col=None):
		zpos=self.image.list_idx
		curpts=[p for p in self.points if p[4]==self.classid]
		
		#print np.array(self.points)
		#print "#########"

		cid=np.unique([p[3] for p in curpts])
		for ci in cid:
			#### draw lines
			pts=np.array([[p[0], p[1], p[2]] for p in curpts if ci==p[3]])
			dzs=pts[:,2]-zpos
			lns=[]
			for i in range(len(pts)-1):
				if dzs[i]*dzs[i+1]<=0: ### line cross current plane
					p0=pts[i]
					p1=pts[i+1]
					if (p1[2]==p0[2]):
						q0=p0
						q1=p1
					else:
						dz=abs(old_div(float(zpos-p0[2]),(p1[2]-p0[2])))
						dp=old_div(1.,abs((p1[2]-p0[2])))
						#print dp, dz
						d0=max(0, (dz-dp))
						d1=min(1, (dz+dp))
						q0=p1*d0+p0*(1-d0)
						q1=p1*d1+p0*(1-d1)

					lns.append([q0[0], q0[1]])
					lns.append([q1[0], q1[1]])

			if ci==self.select:
				clr0=.3
			else:
				clr0=.7

			glColor3f( 1., clr0, clr0 );
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(lns)
			glDrawArrays(GL_LINES, 0, len(lns))


		for p in curpts:
			#### draw nodes on the plane
			#pts=[[p[0], p[1]] for p in self.points if p[2]==zpos]
			s=14.-abs(p[2]-zpos)
			if s<=0: continue
			if p[3]==self.select:
				clr0=.3
			else:
				clr0=.7
			glColor3f( clr0, clr0, 1 );
			glPointSize(s)
			glBegin(GL_POINTS)
			glVertex(p[0], p[1], 0)

			glEnd()
			#glEnableClientState(GL_VERTEX_ARRAY)
			#glVertexPointerf(pts)
			#glDrawArrays(GL_POINTS, 0, len(pts))

		return


class EMDrawWindow(QtWidgets.QMainWindow):

	def __init__(self,application,options,datafile=None):
		QtWidgets.QWidget.__init__(self)
		self.imgview = EMImage2DWidget()
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		
		self.lb_txt0=QtWidgets.QLabel("ClassID")
		self.gbl.addWidget(self.lb_txt0, 0,0,1,1)
		self.classidbox=QtWidgets.QSpinBox()
		self.classidbox.setMinimum(0)
		self.classidbox.setMaximum(5)
		self.classidbox.setValue(0)
		self.gbl.addWidget(self.classidbox, 0,1,1,1)
		
		self.lb_lines=QtWidgets.QLabel("")
		self.lb_lines.setWordWrap(True)
		self.gbl.addWidget(self.lb_lines, 1,0,1,2)
		
		
		self.bt_showimg=QtWidgets.QPushButton("Show tomogram")
		self.bt_showimg.setToolTip("Show tomogram window")
		self.gbl.addWidget(self.bt_showimg, 2,0,1,2)
		
		self.bt_savepdb=QtWidgets.QPushButton("Save PDB")
		self.bt_savepdb.setToolTip("Save curves as PDB")
		self.gbl.addWidget(self.bt_savepdb, 3,0,1,2)
		
		self.bt_clear=QtWidgets.QPushButton("Clear")
		self.bt_clear.setToolTip("Clear all points")
		self.gbl.addWidget(self.bt_clear, 4,0,1,2)
		
		self.bt_interp=QtWidgets.QPushButton("Interpolate")
		self.bt_interp.setToolTip("Interpolate points")
		self.gbl.addWidget(self.bt_interp, 5,0,1,1)
		
		self.tx_interp=QtWidgets.QLineEdit(self)
		self.tx_interp.setText("20")
		self.gbl.addWidget(self.tx_interp, 5,1,1,1)
		
		self.classidbox.valueChanged[int].connect(self.classid_change)
		self.bt_showimg.clicked[bool].connect(self.show_tomo)
		self.bt_savepdb.clicked[bool].connect(self.save_pdb)
		self.bt_interp.clicked[bool].connect(self.interp_points)
		self.bt_clear.clicked[bool].connect(self.clear_points)

		#self.gbl.addWidget(self.imgview,0,0)
		self.options=options
		self.app=weakref.ref(application)
		
		
		
		self.datafile=datafile
		self.data=EMData(datafile)
		self.imgview.setWindowTitle(base_name(datafile))
		self.imgview.list_idx=self.data["nz"]//2
		self.imgview.set_data(self.data)
		self.imgview.show()
		self.infofile=info_name(datafile)


		
		pts=[]
		self.apix_scale=1
		self.tomocenter=np.zeros(3)

		js=js_open_dict(self.infofile)
		if "apix_unbin" in js:
			apix_cur=apix=self.data["apix_x"]
			apix_unbin=js["apix_unbin"]
			self.apix_scale=apix_cur/apix_unbin
			self.tomocenter= np.array([self.data["nx"],self.data["ny"],self.data["nz"]])/2
		
			
		if js.has_key("curves") and len(js["curves"])>0:
			pts=np.array(js["curves"]).copy()
			if len(pts[0])<5:
				pts=np.hstack([pts, np.zeros((len(pts),1))])
			pts[:,:3]=pts[:,:3]/self.apix_scale + self.tomocenter
			pts=pts.tolist()
				
		else:
			pts=[]
		js.close()
		
		self.contour=Contour(img=self.imgview, points=pts )
		self.shape_index = 0
		self.imgview.shapes = {0:self.contour}


		self.imgview.mouseup.connect(self.on_mouseup)
		self.imgview.keypress.connect(self.key_press)
		self.update_label()

		glEnable(GL_POINT_SMOOTH)
		glEnable( GL_LINE_SMOOTH );
		glEnable( GL_POLYGON_SMOOTH );
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	def update_label(self):
		pts=np.array([p for p in self.contour.points if p[4]==self.contour.classid])
		
		if len(pts)==0:
			return
		lb=np.unique(pts[:,3])
		txt="{:d} curves\n".format(len(lb))
		txt+="{:d} points\n".format(len(pts))
		txt+='   '+','.join([str(np.sum(pts[:,3]==i)) for i in lb])
		self.lb_lines.setText(txt)

	def classid_change(self):
		idx=int(self.classidbox.value())
		self.contour.classid=idx
		self.do_update()
		
		
	def show_tomo(self):
		self.imgview.show()
		
	def save_pdb(self):
		pts=np.array(self.contour.points)
		if len(pts)==0:
			return
		
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save PDB', os.getcwd(), "(*.pdb)")[0]
		numpy2pdb(data=pts[:,:3], fname=filename, chainid=pts[:,3])
	
	def interp_points(self):
		nppts=np.array(self.contour.points)
		sel=nppts[:,4]==self.contour.classid
		pts=nppts[sel].copy()
		otherpts=nppts[sel==False].copy()
		
		if len(pts)==0:
			return
		
		try:
			density=float(self.tx_interp.text())
			print("Interpolate to one points per {:.1f} pixel...".format(density))
		except:
			return
			
		pts_intp=[]
		kk=0
		for li in np.unique(pts[:,3]):
			pt=pts[pts[:,3]==li][:,:3].copy()
			pt=pt[np.append(True, np.linalg.norm(np.diff(pt, axis=0), axis=1)>0.1)]
			ln=np.linalg.norm(pt[-1]-pt[0])
			#     print np.round(ln)//2
			if len(pt)<2: continue
			#print len(pt), ln
			ipt=interp_points(pt, npt=np.round(ln)//density)
			
			pts_intp.append(np.hstack([ipt, kk+np.zeros((len(ipt), 1)), self.contour.classid+np.zeros((len(ipt), 1))]))
			kk+=1
		
		pts_intp=np.vstack(pts_intp)
		self.contour.points=otherpts.tolist()
		self.contour.points.extend(pts_intp.tolist())
		self.do_update()
		self.save_points()
		
	def save_points(self):
		js=js_open_dict(self.infofile)
		if self.apix_scale==0:
			js["curves"]=self.contour.points
		else:
			pts=np.array(self.contour.points)
			if len(pts)>0:
				pts[:,:3]=(pts[:,:3]- self.tomocenter) * self.apix_scale
			js["curves"]=pts.tolist()
			
		js.close()
		
		
	def key_press(self, event):
		#print(event.key())
		if event.key()==96:
			self.imgview.increment_list_data(1)
			self.do_update()
		elif event.key()==49:	
			self.imgview.increment_list_data(-1)
			self.do_update()
		return

	def clear_points(self):
		choice = QtWidgets.QMessageBox.question(self, 'Clear points', 
			'Clear all points in the tomogram?', QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
		if choice == QtWidgets.QMessageBox.Yes:
			self.contour.points=[]
		
		self.do_update()
		self.save_points()
		
		return

	def on_mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		if not event.button()&Qt.LeftButton:
			return

		if event.modifiers()&Qt.ControlModifier or event.modifiers()&Qt.ShiftModifier:

			#### check if clicking on an existing point
			pts=np.array(self.contour.points)
			if len(pts)<1:
				return

			ofpln=abs(pts[:,2]-self.imgview.list_idx)>3

			res=np.sqrt(np.sum((pts[:,:2]-np.array([x,y]))**2, axis=1))
			res[ofpln]+=1e5
			ii=np.argmin(res)

			if event.modifiers()&Qt.ShiftModifier:
				if res[ii]<20:
					#### remove point
					self.contour.points=[p for i,p in enumerate(self.contour.points) if i!=ii]

			elif event.modifiers()&Qt.ControlModifier:
				if res[ii]<20:
					#### select contour
					ci=self.contour.points[ii][3]
					print('select contour {:d}'.format(int(ci)))
					self.contour.select=ci
				else:
					print('new contour')
					self.contour.add_point([x, y], True)


		else:
			#### add point
			self.contour.add_point([x, y]) #, self.imgview.list_idx
		
		self.do_update()
		self.save_points()
	
	def do_update(self):
		self.imgview.shapechange=1
		self.update_label()
		self.imgview.updateGL()
		
	def closeEvent(self, event):
		self.imgview.close()

if __name__ == '__main__':
	main()
