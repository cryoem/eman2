#!/usr/bin/env python
# Muyuan Chen 2017-03
from builtins import range
from EMAN2 import *
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
	parser.add_argument("--load", type=str,help="Load previous contour segmentation.", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)",  row=1, col=0,rowspan=1, colspan=2, mode="tomoseg")
	#parser.add_argument("--noupdate", action="store_true",help="do not erase shapes", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()

	if len(args) == 0:
		print("You must specify a tomogram")
		sys.exit(1)

	logid=E2init(sys.argv)
	img = EMData(args[0])

	app = EMApp()

	drawer=EMDrawWindow(app,options,datafile=img)


	drawer.show()
	app.execute()
	E2end(logid)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

class Contour(EMShape):
	def __init__(self, img=None, points=[]):
		self.points=points
		self.isanimated = False
		self.shape=["scr_contour",0]
		self.image=img
		self.triangles=[]
		#self.lines=[]

#		if len(points)>3:

#		self.make_triangle()


	def add_point(self, newpt=[], newcontour=False):



		zpos=self.image.list_idx

		pts=np.array(self.points)
		if len(pts)>=3:
			pts=pts[pts[:,2]==zpos,:]

		if len(pts)<3:
			if len(newpt)>0:
				self.points.append([newpt[0], newpt[1], zpos, 0])

		else:

			thr=1000.



			if newcontour==False:

				cid=np.unique(pts[:,3])

				if len(newpt)>0:
					#### add a point
					dst=np.array(scipydist.cdist([newpt], pts[:,:2])[0])
					rg=np.arange(len(dst), dtype=int)+1
					for ci in cid:
						idx=np.where(pts[:,3]==ci)[0]
						rg[idx[-1]]=idx[0]
					#print rg
					dst+=dst[rg]

					mi=np.argmin(dst)
					ci=pts[mi, 3]
					pts=np.insert(pts, mi+1, np.append(newpt, [zpos, ci]), axis=0)


				allpts=[]
				for ci in cid:
					#### a simple tsp solver...
					idx=np.where(pts[:,3]==ci)[0]
					pp=pts[idx].copy()
					path=np.arange(len(pp), dtype=int)

					dmat=scipydist.squareform(scipydist.pdist(pp[:,:2]))
					dmat+=np.eye(len(dmat))*thr
					if len(pp)>3:
						niter=2000
					else:
						niter=0

					calc_dist=lambda path: np.sum([dmat[path[i], path[(i+1)%len(path)]] for i in range(len(path))])
					dst=calc_dist(path)
					nochange=0
					for k in range(niter):
						p0=path.copy()
						i0,i1= np.sort(np.random.randint(0,len(pp), 2))
						if abs(i0-i1)%len(path)<2: continue
						path= np.hstack([path[:i0+1], path[i0+1:i1+1][::-1], path[i1+1:]])
						d=calc_dist(path)
						if d>=dst:
							path=p0
							nochange+=1
							if nochange>200: break

						else:
							dst=d
							nochange=0


					allpts.extend([[pp[i][0], pp[i][1], zpos, ci] for i in path])

				self.points=[p for p in self.points if p[2]!=zpos]
				self.points.extend(allpts)

			else:
					#### start a new contour
					ci=np.max(pts[:,3])+1
					self.points.append([newpt[0], newpt[1], zpos, ci])

		np.savetxt("pts.save", self.points, fmt='%d')

	def next_slice(self):
		pts=np.array(self.points)
		mi=self.image.list_idx
		ii=np.argmin(abs(pts[:,2]-mi))
		if pts[ii,2]==mi: return
		last=pts[ii,2]
		pts=pts[pts[:,2]==last]
		print(mi, last, pts.shape)
		img=self.image.data.numpy()

		vec=[]
		cid=np.unique(pts[:,3])
		rg0=np.arange(len(pts), dtype=int)-1
		rg1=np.arange(len(pts), dtype=int)+1
		for ci in cid:
			idx=np.where(pts[:,3]==ci)[0]
			rg1[idx[-1]]=idx[0]
			rg0[idx[0]]=idx[-1]

		#for i in range(len(pts)):
		vec= pts[rg0]-pts[rg1]
		vec=np.vstack([vec[:,1], -vec[:,0]]).T
		vec/=np.linalg.norm(vec, axis=1)[:, None]
		#vec.append(v)

		#vec=np.array(vec)

		pval=[]
		rg=np.arange(-2,2.1)
		for i in rg:
			p=(pts[:,:2]-vec*i).astype(int)
			out=(p[:,0]<0)+(p[:,1]<0)+(p[:,0]>=img.shape[0])+(p[:,1]>=img.shape[1])
			out=out>0
			p[out,:]=0
			pval.append(img[p[:,1], p[:,0]])
			pval[-1][out]=1000

		pval=np.array(pval)

		vp=rg[np.argmin(pval, 0)]
		p1=pts.copy()
		p1[:,:2]-=vp[:,None]*vec
		p1[:,:2]=np.round(p1[:,:2])
#		p1[-1]=p1[0]
		p1[:,2]=mi
		self.points.extend(p1.tolist())



	def draw(self,d2s=None,col=None):

		#print "aaaaaaaa"
		zpos=self.image.list_idx
		allpts=[[p[0], p[1], p[3]] for p in self.points if p[2]==zpos]
		print(np.array(self.points))
		print("#########")
		cid=np.unique([p[2] for p in allpts])
		for ci in cid:
			pts=[[p[0], p[1]] for p in allpts if p[2]==ci]
			if len(pts)>2:
				pts.append(pts[0])
				
			area=0.
			for i in range(len(pts)):
				p0=pts[i]
				p1=pts[(i+1)%len(pts)]
				area+=p0[0]*p1[1]-p0[1]*p1[0]
			area=abs(area/2.)
			print("Contour {:d}, area {:.1f}".format(int(ci), area))


			glColor3f( 1, .3, .3 );
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts)
			glDrawArrays(GL_LINE_STRIP, 0, len(pts))


			glColor3f( .3, .3, 1 );
			glPointSize(7.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts)
			glDrawArrays(GL_POINTS, 0, len(pts))
		#print "bbbbbbbbb"

class EMDrawWindow(QtWidgets.QMainWindow):

	def __init__(self,application,options,datafile=None):
		QtWidgets.QWidget.__init__(self)
		self.imgview = EMImage2DWidget()
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())

		self.gbl.addWidget(self.imgview,0,0)
		self.options=options
		self.app=weakref.ref(application)

		self.datafile=datafile
		self.imgview.set_data(datafile)

		#self.all_shapes=[]
		#self.all_points=[]
		pts=[]
		if options.load:
			pts=np.loadtxt(options.load).tolist()

		self.contour=Contour(img=self.imgview, points=pts)
		self.shape_index = 0
		self.imgview.shapes = {0:self.contour}


		self.imgview.mouseup.connect(self.on_mouseup)
		self.imgview.keypress.connect(self.key_press)

		glEnable(GL_POINT_SMOOTH)
		glEnable( GL_LINE_SMOOTH );
		glEnable( GL_POLYGON_SMOOTH );
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	#def update_view(self):
		#self.contour.cur_z=self.imgview.list_idx
		#print self.imgview.list_idx
		#pts=Contour(self.all_points)
		#shps={0:pts}


		#self.imgview.shapes=shps
		#self.imgview.shapechange=1
		#self.imgview.updateGL()

	def key_press(self, event):
		if event.key()==Qt.Key_Shift:
			self.contour.next_slice()
			self.imgview.shapechange=1
			self.imgview.updateGL()

	def on_mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#if event.button()&Qt.LeftButton:

		if event.modifiers()&Qt.ControlModifier:
			#### interpolate curve from previous slices
			print('new contour')
			self.contour.add_point([x, y], True)

		elif event.modifiers()&Qt.ShiftModifier:
			#### remove point
			pts=np.array(self.contour.points)
			pts=pts[pts[:,2]==self.imgview.list_idx,:]
			if len(pts)<1:
				return
			res=np.sqrt(np.sum((pts[:,:2]-np.array([x,y]))**2, axis=1))
			ii=np.argmin(res)
			if res[ii]<20:
				pts=[[p[0], p[1], self.imgview.list_idx, p[3]] for i,p in enumerate(pts) if i!=ii]
				self.contour.points=[p for p in self.contour.points if p[2]!=self.imgview.list_idx]
				self.contour.points.extend(pts)
		else:
			#### add point
			self.contour.add_point([x, y]) #, self.imgview.list_idx
		self.imgview.shapechange=1
		self.imgview.updateGL()

if __name__ == '__main__':
	main()
