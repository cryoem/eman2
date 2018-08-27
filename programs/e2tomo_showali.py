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
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtGui, QtCore, QtOpenGL
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emshape import EMShape
import scipy.spatial.distance as scipydist

def main():

	usage="""This shows the alignment result from e2tomogram.py uing the information from a tomorecon_xx folder.
	
	Example: e2tomo_showali.py --path <tomorecon_xx> --iteration <iteration number>
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path to tomorecon_xx directory to examine fiducial error.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=2, mode="fiderr")
	parser.add_argument("--iteration", type=int,help="Refinement iteration number", default=2,guitype="intbox",row=1, col=0, rowspan=1, colspan=1,mode="fiderr")
	parser.add_argument("--ppid", type=int,help="ppid", default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	fname=os.path.join(options.path, "ali_{:02d}.hdf".format(options.iteration))
	pname=os.path.join(options.path, "landmarks_{:02d}.txt".format(options.iteration))
	pks=np.loadtxt(pname)
	imgs = EMData.read_images(fname)
	img=EMData(imgs[0]["nx"],imgs[0]["ny"], len(imgs))
	xfs=[]
	for i,m in enumerate(imgs):
		img.insert_clip(m, (0,0,i))
		xfs.append(m["xform.projection"])


	aname=os.path.join(options.path, "ptclali_{:02d}.hdf".format(options.iteration))
	n=EMUtil.get_image_count(aname)
	dirs=np.zeros((len(imgs),len(pks),2))
	for i in range(n):
		e=EMData(aname, i, True)
		s=e["score"]
		dirs[e["nid"], e["pid"]]=[s[1], -s[0]]
	print(dirs)

	app = EMApp()

	drawer=EMDrawWindow(app,options,img, pks, xfs, dirs)


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
	def __init__(self, img, pks, xfs, dirs):
		self.isanimated = False
		self.shape=["boxes",0]
		self.image=img
		self.pks=pks
		self.xfs=xfs
		self.dirs=dirs
		#self.triangles=[]



	def draw(self,d2s=None,col=None):
		mi=self.image.list_idx
		xf=self.xfs[mi]
		#print(mi,xf)
		dx=old_div(self.image.data["nx"],2)
		dy=old_div(self.image.data["ny"],2)
		ps=[]
		glColor3f( 1., 1., 1. );
		glPointSize(3)
		for pid, p in enumerate(self.pks):
			pt=[p[0], p[1], p[2]]
			ptx=xf.transform(pt)
			ptx=[old_div(i,2) for i in ptx]
			ptx=[ptx[0]+dx, ptx[1]+dy]
			#print(ptx[0]+dx, ptx[1]+dy)
			ps.append(ptx)
			a=np.array(ptx)+self.dirs[mi, pid]
			ps.append(a.tolist())

			lns=get_circle(ptx, 32)
			#print(lns)
			glColor3f( 1., 1., 1. );
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(lns)
			glDrawArrays(GL_LINES, 0, len(lns))

			#arr=[ptx]
			#a=np.array(pt)+self.dirs[mi, pid]
			#arr.append(a.tolist())
			#print arr
			#glEnableClientState(GL_VERTEX_ARRAY)
			#glVertexPointerf(arr)
			#glDrawArrays(GL_LINES, 0, len(arr))



		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointerf(ps)
		glDrawArrays(GL_LINES, 0, len(ps))
		glPointSize(5)
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointerf([ps[i] for i in range(0,len(ps),2)])
		glDrawArrays(GL_POINTS, 0, old_div(len(ps),2))

		return


class EMDrawWindow(QtGui.QMainWindow):

	def __init__(self,application,options,datafile, pks, xfs, dirs):
		QtGui.QWidget.__init__(self)
		self.imgview = EMImage2DWidget()
		self.setCentralWidget(QtGui.QWidget())
		self.gbl = QtGui.QGridLayout(self.centralWidget())

		self.gbl.addWidget(self.imgview,0,0)
		self.options=options
		self.app=weakref.ref(application)

		self.datafile=datafile
		self.imgview.set_data(datafile)

		#self.all_shapes=[]
		#self.all_points=[]

		self.boxes=Boxes(self.imgview, pks, xfs, dirs)
		self.shape_index = 0
		self.imgview.shapes = {0:self.boxes}



		glEnable(GL_POINT_SMOOTH)
		glEnable( GL_LINE_SMOOTH );
		glEnable( GL_POLYGON_SMOOTH );
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);



if __name__ == '__main__':
	main()
