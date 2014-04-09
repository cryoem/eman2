#!/usr/bin/env python

#
# Author: David Woolford 11/25/08 (woolford@bcm.edu
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
from EMAN2db import db_open_dict, db_check_dict
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtGui,QtCore
from PyQt4.QtCore import Qt
from emanimationutil import OrientationListAnimation,Animator
from emapplication import EMApp, get_application, error
from emglobjects import EM3DModel
from emimage2d import EMImage2DWidget
from emimage3dsym import EM3DSymModel, EMSymInspector, EMSymViewerWidget
from emimagemx import EMImageMXWidget, EMLightWeightParticleCache
import os
import sys
import weakref

def get_eulers_from(filename):
	eulers = []
	n = EMUtil.get_image_count(filename)
	for i in range(n):
		h = get_header(filename,i)
		try: p = h["xform.projection"]
		except:
			continue
			#print "image",i,"doesn't have the xform.projection attribute"
			#return None

		eulers.append(p)

	return eulers

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog

	Presents a graphical representation of the orientation distribution of the particles in a single particle
	reconstruction. This is displayed as a single asymmetric unit on a sphere, with cylinders of varying height
	representing the number of particles found in each orientation. Middle-click will produce a control-panel
	as usual, and clicking on a single peak will permit viewing the class-average and related projection
	(and particles).

	Normally launched without arguments from a e2workflow project directory.

"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_header(name="e2eulerxplorheader", help="Click Launch to Run e2eulerxplor.py", title="### Click Launch to Run e2eulerxplor.py ###", row=0, col=0, rowspan=1, colspan=3)
	parser.add_argument("--eulerdata", "-e", type=str,help="File for Eulerdata, Ryan style, if none is given, data is read from the DB.",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	global options
	(options, args) = parser.parse_args()


	logid=E2init(sys.argv,options.ppid)

	em_app = EMApp()
	window = EMEulerWidget(file_name = options.eulerdata)
	em_app.show_specific(window)
	em_app.execute()

	E2end(logid)

class EMEulerWidget(EMSymViewerWidget):
	def __init__(self, auto=True,sparse_mode=False, file_name = ""):
		EMSymViewerWidget.__init__(self, filename=file_name)

		#Replacing the EM3DSymModel that was created in the base class
		euler_explorer = EMEulerExplorer(self, auto, sparse_mode, file_name = file_name)
		self.model = euler_explorer
		euler_explorer.regen_dl()

def sadd(d,a,b):
	if a==None : return None
	if b==None : return d+"/"+a
	return d+"/"+a+b

class EMEulerExplorer(EM3DSymModel,Animator):

	def mousePressEvent(self,event):
		if self.events_mode == "inspect":
			self.current_hit = self.get_hit(event)
			if self.current_hit == None:
				EM3DSymModel.mousePressEvent(self,event)
		else:
			EM3DSymModel.mousePressEvent(self,event)

	def mouseReleaseEvent(self,event):
		if self.events_mode == "inspect":
			if self.current_hit != None:
				self.updateGL() # there needs to be a clear or something  in order for the picking to work. This is  bit of hack but our rendering function doesn't take long anyhow
				hit = self.get_hit(event)
				if hit == self.current_hit:
					self.emit(QtCore.SIGNAL("point_selected"),self.current_hit,event)
			else:
				#EM3DSymModel.mouseReleaseEvent(self,event)
				EM3DModel.mouseReleaseEvent(self, event) #behavior in EM3DSymModel is not what we want (needed in sibling classes?)

			self.current_hit = None
		else:
				#EM3DSymModel.mouseReleaseEvent(self,event)
				EM3DModel.mouseReleaseEvent(self, event) #behavior in EM3DSymModel is not what we want (needed in sibling classes?)


	def mouseMoveEvent(self,event):
		if self.events_mode == "inspect" and self.current_hit:
			pass
		else:
			EM3DSymModel.mouseMoveEvent(self,event)

	def get_hit(self,event):
		v = self.vdtools.wview.tolist()
		self.get_gl_widget().makeCurrent() # prevents a stack underflow
#		x = event.x()
#		y = v[-1]-event.y()
#		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
#		vals = self.render(color_picking=True)
#		glFlush()
#		vv = glReadPixels(x,y,1,1,GL_RGB,GL_FLOAT)
#		reslt = Vec3f(float(vv[0][0][0]),float(vv[0][0][1]),float(vv[0][0][2]))
#		for i,val in enumerate(vals):
##			print val,reslt,(reslt-val).length(),vv[0][0]
#			if (reslt-val).length() < 0.01:
#				print i
##				print (reslt-val).length()
#				return i
#		print vv
#
		# the problem with this approach is that depth testing is not part of picking
		sb = [0 for i in xrange(0,512)]
		glSelectBuffer(512)
		glRenderMode(GL_SELECT)
		glInitNames()
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		gluPickMatrix(event.x(),v[-1]-event.y(),5,5,v)
		self.get_gl_widget().load_perspective()
		glMatrixMode(GL_MODELVIEW)
		glInitNames()
		self.render()
		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		glFlush()

		intersection = None
		hits = list(glRenderMode(GL_RENDER))
		for hit in hits:
			a,b,c=hit
			if len(c) > 0:
				intersection = c[0]-1
				break

		return intersection


	def keyPressEvent(self,event):

		if event.key() == Qt.Key_F1:
			self.display_web_help("http://blake.bcm.edu/emanwiki/EMAN2/Programs/e2eulerxplor")
		elif event.key() == Qt.Key_F :
			if self.flatten>0 : self.flatten=0.0
			else: self.flatten=1.0
			self.generate_current_display_list(True)
			self.updateGL()
		else:
			EM3DSymModel.keyPressEvent(self,event)

	def __init__(self, gl_widget=None, auto=True,sparse_mode=False, file_name = ""):
		self.current_hit = None
		self.events_mode_list = ["navigate", "inspect"]
		self.events_mode = self.events_mode_list[1]


		self.init_lock = True # a lock indicated that we are still in the __init__ function
		self.au_data = None # This will be a dictionary, keys will be refinement directories, values will be something like available iterations for visual study
		if auto: # this a flag that tells the eulerxplorer to search for refinement data and automatically add elements to the inspector, if so
			self.gen_refinement_data()

		EM3DSymModel.__init__(self, gl_widget, eulerfilename=file_name)
		Animator.__init__(self)
		self.height_scale = 8.0 # This is a value used in EM3DSymModel which scales the height of the displayed cylinders - I made it 8 because it seemed fine. The user can change it anyhow
		self.projection_file = None  # This is a string - the name of the projection images file
		self.average_file = None # This is a string - the name of the class averages file
		self.proj_class_viewer = None # This will be an EMImageMXWidget that shows the class and/or projection
		self.particle_viewer = None  # This will be an EMImageMXWidget that shows the particles in a class
		self.clsdb = None # I think this will become redundant - it used to be the old database that stores which particles are in a class, but now that's stored in the header
		self.particle_file = None # This will be a string - the name of the file that has the particle files in it. This might be made redundant with the new approach
		self.alignment_file = None # This will be a string - the name of the file containing the alignment parameters - this is essential if you we want to show the aligned particles
		self.refine_dir = None # This will be a string - the name of the current refinement directory that is being studied
		self.dx = None # This is an EMData object storing the x shifts of the alignments for all particles. Generated by e2classaverage
		self.dy = None # This is an EMData object storing the y shifts of the alignments for all particles. Generated by e2classaverage
		self.da = None# This is an EMData object storing the angle of the alignments for all particles. Generated by e2classaverage
		self.dflip = None # This is an EMData object storing whether or not tthe alignment involved a flip, for all particles. Generated by e2classaverage
		self.classes = None # This is an EMData object storing which class(es) a particle belongs to. Generated by e2classaverage
		self.inclusions = None # This is and EMDAta storing a boolean that indicates the particle was actually included in the final average. Generated by e2classaverage

		self.average = None # This the class average itself, an EMData object
		self.projection = None # This is the projection itelse, an EMData object
		self.class_idx = None # This is the idx of the current class being studied in the interface

		self.previous_len = -1 # To keep track of the number of class averages that were previously viewable. This helps to make sure we can switch to the same class average in the context of a different refinement iteration
		self.mirror_eulers = False
		if sparse_mode:
			self.mirror_eulers = True # If True the drawn Eulers are are also rendered on the opposite side of the sphere - see EM3DSymModel.make_sym_dl_lis

		# Grab the symmetry from the workflow database if possible
		sym = "c1"
		if js_check_dict("refine_01/0_refine_parms.json"):
			try: sym = str(js_open_dict("refine_01/0_refine_parms.json")["sym"])
			except: pass

		# Have to tell the EM3DSymModel that there is a new sym
		self.set_symmetry(sym)

		# this object will have
		if self.au_data != None:
			combo_entries = self.au_data.keys()
			combo_entries.sort()
			combo_entries.reverse()

			if len(combo_entries) > 0:
				au = combo_entries[0]
				cls = self.au_data[au][0][0]
				self.au_selected(au,cls)
				self.mirror_eulers = True

		self.init_lock = False
		self.force_update=True # Force a display udpdate in EMImage3DSymModule

		QtCore.QObject.connect(self, QtCore.SIGNAL("point_selected"), self.au_point_selected)

	def __del__(self):
		EM3DSymModel.__del__(self) # this is here for documentation purposes - beware that the del function is important

	def initializeGL(self):
		glEnable(GL_NORMALIZE)

	def generate_current_display_list(self,force=False):
		'''
		Redefinition of EMImage3DSymModule.generate_current_display_list

		'''
		if self.init_lock: return 0
		if self.au_data == None or len(self.au_data) == 0:
			EM3DSymModel.generate_current_display_list(self,force)

		self.init_basic_shapes()
		if self.nomirror == True : val = 0
		else: val = 1
		self.trace_great_arcs(self.sym_object.get_asym_unit_points(val))
		self.trace_great_triangles(val)

		self.eulers = self.specified_eulers
		if self.eulers == None:	return 0

#		if not self.colors_specified: self.point_colors = []
#		else: self.point_colors = self.specified_colors
#		self.points = []
#		for i in self.eulers:
#			p = i.transpose()*Vec3f(0,0,self.radius)
#			self.points.append(p)
#			if not self.colors_specified: self.point_colors.append((0.34615, 0.3143, 0.0903,1))

		self.make_sym_dl_list(self.eulers)
		return 1
	def get_data_dims(self):
		return (2*self.radius,2*self.radius,2*self.radius)

	def width(self): return 2*self.radius
	def height(self): return 2*self.radius
	def depth(self): return 2*self.radius

	def gen_refinement_data(self):
		dirs,files = get_files_and_directories()

		dirs.sort()
		for i in range(len(dirs)-1,-1,-1):
			if len(dirs[i]) != 9:
				dirs.pop(i)
			elif dirs[i][:7] != "refine_":
				dirs.pop(i)
			else:
				try: int(dirs[i][7:])
				except: dirs.pop(i)

		self.dirs = dirs
		print dirs

		self.au_data = {}
		for dir in self.dirs:
			d = self.check_refine_db_dir(dir)
			if len(d) != 0 and len(d[dir]) != 0: self.au_data.update(d)

	def check_refine_db_dir(self,dir,s1="classes",s2=None,s3="cls_result",s4="threed",s5="projections"):
		# s2 used to be class_indices
		names = [s1,s2,s3,s4,s5]
		data = {}
		data[dir] = []
		register_js_name = "{}/0_refine_parms.json".format(dir)

		files=os.listdir(dir)
		try:
			nums=[int(i[7:9]) for i in files if "threed" in i and "even" not in i and "odd" not in i]
			maxnum=max(nums)
		except :
			print "Nothing in ",dir
			return {}

		for i in xrange(1,maxnum+1):
			exte="_{:02d}_even.hdf".format(i)
			exto="_{:02d}_odd.hdf".format(i)
			data[dir].append([sadd(dir,s1,exte),sadd(dir,s2,exte),sadd(dir,s3,exte),sadd(dir,s4,exte),sadd(dir,s5,exte)])
			data[dir].append([sadd(dir,s1,exto),sadd(dir,s2,exto),sadd(dir,s3,exto),sadd(dir,s4,exto),sadd(dir,s5,exto)])

		return data

	def set_projection_file(self,projection_file): self.projection_file = projection_file
	def get_inspector(self):
		if not self.inspector :
			if (self.au_data == None or len(self.au_data) == 0) and self.mirror_eulers == False: #self.mirror_eulers thing is a little bit of a hack, it's tied to the sparse_mode flag in the init function, which is used by euler_display in EMAN2.py
				self.inspector=EMAsymmetricUnitInspector(self,True,True)
			else:
				self.inspector=EMAsymmetricUnitInspector(self)
			QtCore.QObject.connect(self.inspector,QtCore.SIGNAL("au_selected"),self.au_selected)
		return self.inspector


	def au_selected(self,refine_dir,cls):
		self.refine_dir = refine_dir
		get_application().setOverrideCursor(Qt.BusyCursor)
		data = []
		for d in self.au_data[refine_dir]:
			if d[0] == cls:
				data = d;
				break

		if len(data) == 0:
			error("error, no data for %s %s, returning" %(refine_dir,cls))
#			print "error, no data for",au,cls,"returning"
			self.events_handlers["inspect"].reset()
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return

		try :
			self.particle_file=js_open_dict(refine_dir+"/0_refine_parms.json")["input"]
		except:
			error("No data in "+refine_dir )
			self.events_handlers["inspect"].reset()
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return

		self.average_file = cls
		self.projection_file = data[4]
		self.alignment_file = data[2]
		self.clsdb = data[1]

		self.dx = None
		self.dy = None
		self.da = None
		self.dflip = None
		self.classes = None

		eulers = get_eulers_from(self.average_file)
		#s = Symmetries.get("d7")
		#eulers = s.gen_orientations("rand",{"n":EMUtil.get_image_count(self.average_file)})

		self.specify_eulers(eulers)
		#from emimagemx import EMDataListCache
		#a = EMData.read_images(self.average_file)
		#a = [test_image() for i in range(EMUtil.get_image_count(self.average_file))]
		#print len(a),len(eulers)
		#b = [a[i].set_attr("xform.projection",eulers[i]) for i in range(len(eulers))]
		#b = [a[i].set_attr("ptcl_repr",1) for i in range(len(eulers))]

		self.set_emdata_list_as_data(EMLightWeightParticleCache.from_file(self.average_file),"ptcl_repr")
		#self.set_emdata_list_as_data(EMDataListCache(self.average_file),"ptcl_repr")
#		self.set_emdata_list_as_data(a,"ptcl_repr")
		self.force_update = True
		self.au_point_selected(self.class_idx,None)
		# if we have the same number of Eulers we can update everything
#		if self.previous_len == len(eulers) : self.events_handlers["inspect"].repeat_event()
#		else:self.events_handlers["inspect"].reset()
		self.previous_len = len(eulers)
		if not self.init_lock:self.updateGL()
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def __get_file_headers(self,filename):
		headers = []
		n = EMUtil.get_image_count(filename)
		for i in range(n):
			e = EMData()
			e.read_image(filename,i,True)
			headers.append(e)
		return headers

	def au_point_selected(self,i,event=None):
		if i == None:
			if event != None and event.modifiers()&Qt.ShiftModifier:
				if self.special_euler != None:
					self.special_euler = None
					if not self.init_lock:self.regen_dl()
			return
#		self.arc_anim_points = None
		self.projection = None
		if self.euler_data:
#			db = db_open_dict(self.average_file)
#			a = db.get(i)
#			print a["nx"]
#			print self.average_file,i
#			self.average = EMData(self.average_file,i)
#			self.average["nx"]
			self.average = self.euler_data[i]#
			self.projection = EMData(self.projection_file,self.average.get_attr("projection_image_idx"))
			self.average.process_inplace("normalize.toimage",{"to":self.projection})
			try:
				self.class_idx = self.average.get_attr("projection_image_idx")
				print "%d (%d)"%(self.class_idx,self.average["ptcl_repr"])
			except:
				self.class_idx = -1
		else: return

		#if self.projection  == None and self.average == None: return
		first = False
		if self.proj_class_viewer == None:
			first = True
			self.proj_class_viewer = EMImageMXWidget(data=None,application=get_application())
#			self.proj_class_viewer = EMImage2DWidget(image=None,application=get_application())
			QtCore.QObject.connect(self.proj_class_viewer,QtCore.SIGNAL("module_closed"),self.on_mx_view_closed)
#			self.proj_class_viewer.set_mouse_mode("App" )
			QtCore.QObject.connect(self.proj_class_viewer,QtCore.SIGNAL("mx_image_selected"), self.mx_image_selected)
			get_application().show_specific(self.proj_class_viewer)

			self.proj_class_single = EMImage2DWidget(image=None,application=get_application())
			QtCore.QObject.connect(self.proj_class_single,QtCore.SIGNAL("module_closed"),self.on_mx_view_closed)
#			QtCore.QObject.connect(self.proj_class_single,QtCore.SIGNAL("mx_image_selected"), self.mx_image_selected)
			get_application().show_specific(self.proj_class_single)

		disp = []
		if self.projection != None: disp.append(self.projection)
		if self.average != None and self.projection!=None:
			# ok, this really should be put into its own processor
			#dataf = self.projection.do_fft()
			#apix=self.projection["apix_x"]
			#curve = dataf.calc_radial_dist(dataf["ny"], 0, 0.5,True)
			#curve=[i/(dataf["nx"]*dataf["ny"])**2 for i in curve]
			#xcurve=[i/(apix*2.0*dataf["ny"]) for i in range(len(curve))]
			#xyd=XYData()
			#xyd.set_xy_list(xcurve,curve)
			#filt=self.average.process("filter.setstrucfac",{"apix":apix,"strucfac":xyd})
			#filt.process_inplace("normalize.toimage",{"to":self.average})
			self.projection["apix_x"]=self.average["apix_x"]
			self.projection["apix_y"]=self.average["apix_y"]
			self.projection["apix_z"]=self.average["apix_z"]
			filt=self.projection.process("threshold.notzero")
			filt.mult(self.average)
			filt.process_inplace("filter.matchto",{"to":self.projection})

			disp.append(filt)

		if self.average!=None:
			disp.append(self.average)

		self.proj_class_viewer.set_data(disp)
		self.proj_class_single.set_data(disp)

		self.proj_class_viewer.updateGL()
		self.proj_class_single.updateGL()
		if self.particle_viewer != None:
			self.mx_image_selected(None,None)
		if first: self.proj_class_viewer.optimally_resize()

		if i != self.special_euler:
			self.special_euler = i
			self.force_update = True

		if not self.init_lock: self.updateGL()


	def on_mx_view_closed(self):
		self.proj_class_viewer = None
		self.proj_class_single = None

	def on_particle_mx_view_closed(self):
		self.particle_viewer = None

	def animation_done_event(self,animation):
		pass

	def alignment_time_animation(self,transforms):
		if len(transforms) < 2: return
		animation = OrientationListAnimation(self,transforms,self.radius)
		self.register_animatable(animation)

	def particle_selected(self,event,lc):
		if lc != None:
			d = lc[3]
			ptcl_idx = d["Img #"]
			data = self.au_data[self.refine_dir]
			prj = []
			cls_result = []
			for l in data:
				for s in l:
					stag = base_name(s)

					if len(stag) > 11 and stag[:11] == "projections":
						prj.append(s)
					elif len(stag) > 10 and stag[:10] == "cls_result":
						cls_result.append(s)

			transforms = []
			if len(prj) != len(cls_result): RunTimeError("The number of cls_result files does not match the number of projection files?")

			e = EMData()
			for i,cr in enumerate(cls_result):
				r = Region(0,ptcl_idx,1,1)
				e.read_image(cr,0,False,r)
				p = int(e.get(0))
				e.read_image(prj[i],p,True)
				transforms.append(e["xform.projection"])

			self.alignment_time_animation(transforms)

	def mx_image_selected(self,event,lc):
#		self.arc_anim_points = None
		get_application().setOverrideCursor(Qt.BusyCursor)
		if lc != None: self.sel = lc[0]

		if self.average != None:
			included = []
			if self.average.has_attr("class_ptcl_idxs"):
				included = self.average["class_ptcl_idxs"]
			excluded = []
			if self.average.has_attr("exc_class_ptcl_idxs"):
				excluded = self.average["exc_class_ptcl_idxs"]

			all = included + excluded
			#all.sort()

			bdata = []
			data = []
			idx_included = []
			running_idx = 0
			from emimagemx import ApplyAttribute
			for val in included:
				bdata.append([self.particle_file,val,[ApplyAttribute("Img #",val)]])
				idx_included.append(running_idx)
				running_idx += 1

			idx_excluded = []
			for val in excluded:
				bdata.append([self.particle_file,val,[ApplyAttribute("Img #",val)]])
				idx_excluded.append(running_idx)
				running_idx += 1

			data = EMLightWeightParticleCache(bdata)

			first = False
			if self.particle_viewer == None:
				first = True
				self.particle_viewer = EMImageMXWidget(data=None,application=get_application())
				self.particle_viewer.set_mouse_mode("App" )
				QtCore.QObject.connect(self.particle_viewer,QtCore.SIGNAL("module_closed"),self.on_particle_mx_view_closed)
				QtCore.QObject.connect(self.particle_viewer,QtCore.SIGNAL("mx_image_selected"), self.particle_selected)
				get_application().show_specific(self.particle_viewer)


			self.check_images_in_memory()

			if self.sel== 0 or self.alignment_file == None:
				self.particle_viewer.set_data(data)
			else:

				for i,[name,idx,f] in enumerate(bdata):
					index = -1
					if self.classes.get_xsize() == 1:
						index = 0 # just assume it's the first one - this is potentially fatal assumption, but in obscure situations only
					else:
						for j in range(self.classes.get_xsize()):
							if int(self.classes.get(j,idx)) == self.class_idx:
								index = j
								break
					if index == -1:
						print "couldn't find"
						get_application().setOverrideCursor(Qt.ArrowCursor)
						return

					x = self.dx.get(index,idx)
					y = self.dy.get(index,idx)
					a = self.da.get(index,idx)
					m = self.dflip.get(index,idx)

					t = Transform({"type":"2d","alpha":a,"mirror":int(m)})
					t.set_trans(x,y)
					from emimagemx import ApplyTransform
					f.append(ApplyTransform(t))
					#data[i].transform(t)
				self.particle_viewer.set_data(data)


			if first:
				self.particle_viewer.updateGL()
				self.particle_viewer.optimally_resize()

			self.particle_viewer.clear_sets(False)
			self.particle_viewer.enable_set("Excluded",idx_excluded,True,False)
			self.particle_viewer.enable_set("Included",idx_included,False,False)
			self.particle_viewer.updateGL()

			get_application().setOverrideCursor(Qt.ArrowCursor)

			self.updateGL()

	def check_images_in_memory(self):
		if self.alignment_file != None:
			if self.dx == None:
				self.dx = EMData(self.alignment_file,2)
			if self.dy == None:
				self.dy = EMData(self.alignment_file,3)
			if self.da == None:
				self.da = EMData(self.alignment_file,4)
			if self.dflip == None:
				self.dflip  = EMData(self.alignment_file,5)
			if self.classes == None:
				self.classes  = EMData(self.alignment_file,0)
			if self.inclusions == None:
				self.inclusions  = EMData(self.alignment_file,1)


	def set_events_mode(self,mode):
		if not mode in self.events_mode_list:
			print "error, unknown events mode", mode
			return

		else:
			self.events_mode = mode

	def closeEvent(self,event):
		if self.inspector !=None: self.inspector.close()
		if self.proj_class_viewer !=None: self.proj_class_viewer.close()
		if self.proj_class_single !=None: self.proj_class_single.close()
		if self.particle_viewer != None: self.particle_viewer.close()
		get_application().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is


def set_included_1(e):
	e.set_attr("included",1)
	e.mxset = [1]

def set_included_0(e):
	e.set_attr("included",0)
	e.mxset = [0]



class EMAsymmetricUnitInspector(EMSymInspector):
	def __init__(self,target,enable_trace=False,enable_og=False) :
		EMSymInspector.__init__(self,target,enable_trace=enable_trace,enable_og=enable_og)

		if target.au_data != None and len(target.au_data) > 0:
			self.add_au_table()

	def add_au_table(self):

		self.au_tab= QtGui.QWidget()
		self.au_tab.vbl = QtGui.QVBoxLayout(self.au_tab)

		self.au_data = self.target().au_data
		combo_entries = self.au_data.keys()
		combo_entries.sort()
		combo_entries.reverse()
		self.combo = QtGui.QComboBox(self)
		for e in combo_entries:
			self.combo.addItem(e)

		self.connect(self.combo,QtCore.SIGNAL("currentIndexChanged(QString&)"),self.on_combo_change)
		self.connect(self.combo,QtCore.SIGNAL("currentIndexChanged(const QString&)"),self.on_combo_change)

		self.au_tab.vbl.addWidget(self.combo)
		self.refine_dir = combo_entries[0]

		self.list_widget = QtGui.QListWidget(None)

		self.list_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
		self.list_widget.setMouseTracking(True)
		QtCore.QObject.connect(self.list_widget,QtCore.SIGNAL("itemClicked(QListWidgetItem *)"),self.list_widget_item_clicked)

		self.update_classes_list(first_time=True)
		self.au_tab.vbl.addWidget(self.list_widget)
		self.tabwidget.insertTab(0,self.au_tab,"Refinement")
		self.tabwidget.setCurrentIndex(0)

	def on_combo_change(self,s):
		self.refine_dir = str(s)
		self.update_classes_list()

	def update_classes_list(self,first_time=False):
		selected_items = self.list_widget.selectedItems() # need to preserve the selection

		s_text = None
		if len(selected_items) == 1 :
			s_text = str(selected_items[0].text())
			if len(s_text) > 4: s_text = s_text[-4:]

		self.list_widget.clear()
		for i,vals in enumerate(self.au_data[self.refine_dir]):
			choice = vals[0]

			a = QtGui.QListWidgetItem(str(choice),self.list_widget)
			if first_time and i == 0:
				self.list_widget.setItemSelected(a,True)
			elif len(choice) > 4 and (choice[-4:] == s_text):
				self.list_widget.setItemSelected(a,True)

		selected_items = self.list_widget.selectedItems() # need to preserve the selection
		if len(selected_items) == 1:
			self.emit(QtCore.SIGNAL("au_selected"),self.refine_dir,str(selected_items[0].text()))

	def list_widget_item_clicked(self,item):
		self.emit(QtCore.SIGNAL("au_selected"),self.refine_dir,str(item.text()))

	def on_mouse_mode_clicked(self,bool):
		for button in self.mouse_mode_buttons:
			if button.isChecked():
				self.target().set_events_mode(str(button.text()))



if __name__ == '__main__':
	main()
