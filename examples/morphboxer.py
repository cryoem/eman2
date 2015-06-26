#!/usr/bin/env python

#
# Authors: James Michael Bell, 06/03/2015
# Copyright (c) 2015 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	2111-1307 USA
#

from EMAN2 import *
import sys
from emboxerbase import *
import subprocess as sp

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image> <image2>....
	
	Auto Boxing strategy making use of mathematical morphology. It is still in need of some work.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="micrographs",help="List the file to process with e2boxer here.", default="", guitype='filebox', browser="EMBoxesTable(withmodal=True,multiselect=True)",row=0,col=0,rowspan=1,colspan=3,mode="boxing,extraction")
	parser.add_header(name="boxerheader", help='Options below this label are specific to e2boxer', title="### e2boxer options ###", row=1, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--boxsize","-b",type=int,help="Box size in pixels",default=-1, guitype='intbox', row=2, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--xmin",type=int,default=0)
	parser.add_argument("--xmax",type=int,default=-1)
	parser.add_argument("--xstep",type=int,default=16)
	parser.add_argument("--ymin",type=int,default=0)
	parser.add_argument("--ymax",type=int,default=-1)
	parser.add_argument("--ystep",type=int,default=16)
	(options, args) = parser.parse_args()

	app = EMApp()
	boxer = AutoBoxer(args,options)
	boxer.show_interfaces()
	app.execute()
	
	#FilterTester(args[1]).display_processing()


class AutoBoxer(EMBoxerModule):
	
	"""
	Class to automatically pick particles from an cryo-EM micrograph
	"""
	
	def __init__(self,micrographs,options):
		super(AutoBoxer,self).__init__(micrographs,options.boxsize)
		self.box_list = MorphBoxList(self)
		
		self.add_tool(MorphBoxingTool)
		
		for i, mgph in enumerate(micrographs):
			self.set_current_file_by_idx(i)
			f = self.current_file()
			self.tile_micrograph(options)
			
	def tile_micrograph(self,options,type="morphological"):
		f = self.current_file()
		hdr = EMData(f,0,True).get_attr_dict()
		nx = hdr['nx']
		ny = hdr['ny']
		if options.xmin == 0: fx = int(options.boxsize)
		if options.ymin == 0: fy = int(options.boxsize)
		if options.xmax == -1: lx = int(nx-options.boxsize)
		if options.ymax == -1: ly = int(ny-options.boxsize)
		boxes = []
		for y in xrange(fy,ly,options.ystep):
			for x in xrange(fx,lx,options.xstep):
				boxes.append([x,y,type])
		self.add_boxes(boxes)
		self.nboxes = len(boxes)
	
	def morph_tiles(self,options):
		# apply filters in filter test to each box individually
		pass
	
	def classify_tiles(self,options):
		# perform a classification to distinguish particles from nonsense
		pass
	
	def remove_nonsense(self,options):
		# removes bad particles (maybe store in a separate file for easy reference?)
		pass
	
	def adjust_filters(self):
		# Just like the filtertool. Allow user to adjust all parameters in the morph panel.
		pass

class MorphBoxList(EMBoxList):

	def __init__(self,target):
		super(MorphBoxList,self).__init__(target)

	def get_particle_images(self,image_name,box_size):
		filtered = []
		for box in self.boxes:
			im = box.get_image(image_name,box_size,"normalize.edgemean")
			fl = self.filter(im)
			filtered.append(fl)
		return filtered #[box.get_image(image_name,box_size,"normalize.edgemean") for box in self.boxes]

	def filter(self,image):
		img = image.copy()
		img.process_inplace('filter.highpass.gauss',{'cutoff_abs':0.02})
		img.process_inplace('normalize.edgemean')
		#img.process_inplace('normalize.maxmin')
		img.process_inplace('math.sigma',{'value1':15.0,'value2':0.0})
		img.process_inplace('normalize.edgemean')
		#img.process_inplace('normalize.maxmin')
		img.process_inplace('filter.lowpass.gauss',{'cutoff_abs':0.088})
		#img.process_inplace('normalize.maxmin')
		#img.process_inplace('normalize.edgemean')
		img.process_inplace('threshold.belowtozero',{'minval':0.098})
		#img.process_inplace('normalize.maxmin')
		#img.process_inplace('normalize.edgemean')
		img.process_inplace('morph.majority',{'nmaj':1,'thresh':0.0})
		#img.process_inplace('normalize.maxmin')
		#img.process_inplace('normalize.edgemean')
		img.process_inplace('morph.object.density',{'thresh':0.0})
		#img.process_inplace('normalize.maxmin')
		#img.process_inplace('normalize.edgemean')
		#img.process_inplace('mask.addshells.multilevel',{'nshells':2})
		#img.process_inplace('normalize.maxmin')
		#img.process_inplace('normalize.edgemean')
		img.process_inplace('histogram.bin',{'nbins':3})
		img.process_inplace('normalize.maxmin')
		return img

class MorphBoxingTool(EMBoxingTool):

	BOX_TYPE = "morphological"
	EMBox.set_box_color(BOX_TYPE,[1,1,1])

	def __init__(self,target):
		super(MorphBoxingTool,self).__init__(target)
		self.target = weakref.ref(target)
		self.moving = None
		self.panel_object = None
		self.moving_data = None

	def get_widget(self):
		if self.panel_object == None:
			self.panel_object = MorphBoxingPanel(self)
		return self.panel_object.get_widget()

	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "white_box.png")

	def set_panel_object(self,panel): self.panel_object = panel
	def unique_name(self): return MorphBoxingTool.BOX_TYPE

	def set_current_file(self,file_name,active_tool=False):
		'''
		If the behavior of this Handler does not if the file changes, but the function needs to be supplied
		'''
		pass

	def get_2d_window(self): return self.target().get_2d_window()

	def mouse_down(self,event) :
		m = self.get_2d_window().scr_to_img((event.x(),event.y()))
		box_num = self.target().detect_box_collision(m)
		from PyQt4.QtCore import Qt
		if box_num == -1:
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			box_num = self.target().add_box(m[0],m[1],MorphBoxingTool.BOX_TYPE)
			if self.panel_object.auto_center_checkbox.isChecked():
				self.try_to_center_ref(box_num)

			self.moving=[m,box_num]
		else:
			box = self.target().get_box(box_num)
			if box.type == MorphBoxingTool.BOX_TYPE:
		 		if event.modifiers()&Qt.ShiftModifier :
					self.target().remove_box(box_num)
				else:
					# if we make it here than the we're moving a box
					self.moving=[m,box_num]
					#self.target().moving_box_established(box_num)
			else:
				raise EMUnknownBoxType,box.type

	def mouse_drag(self,event) :
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.target().detect_box_collision(m)
			if ( box_num != -1):
				box = self.target().get_box(box_num)
				if box.type ==  MorphBoxingTool.BOX_TYPE:
					self.target().remove_box(box_num)
				else:
					raise EMUnknownBoxType,box.type

		elif self.moving != None:
			oldm = self.moving[0]

			self.target().move_box(self.moving[1],m[0]-oldm[0],m[1]-oldm[1])
			self.moving[0] = m

	def mouse_up(self,event) :
		if self.moving != None:
			self.target().box_released(self.moving[1])
		self.moving=None

	def mouse_move(self,event): pass

	def clear_all(self):
		self.target().clear_boxes([MorphBoxingTool.BOX_TYPE],cache=True)

	def moving_ptcl_established(self,box_num,x,y):
		box = self.target().get_box(box_num)
		if box.type != MorphBoxingTool.BOX_TYPE:
			raise EMUnknownBoxType,box.type

		self.moving_data = [x,y,box_num]

	def move_ptcl(self,box_num,x,y,scale):
		if self.moving_data == None: return
		dx = self.moving_data[0] - x
		dy = y - self.moving_data[1]
		self.target().move_box(self.moving_data[2],dx,dy)

		self.moving_data = [x,y,self.moving_data[2]]

	def release_moving_ptcl(self,box_num,x,y):
		if self.moving_data == None: return
		self.target().box_placement_update_exclusion_image_n(box_num)
		self.moving_data = None

	def delete_ptcl(self,box_num):
		box = self.target().get_box(box_num)
		if box.type != MorphBoxingTool.BOX_TYPE:
			raise EMUnknownBoxType,box.type
		self.target().remove_box(box_num)

	def get_unique_box_types(self):
		return [MorphBoxingTool.BOX_TYPE]

	def boxes_erased(self,list_of_boxes):
		'''
		No need to act here for the morph boxing tool - everything is fine?
		'''
		pass

	#TODO: better code reuse, not copy and paste, here
	#COPIED FROM e2boxer's SwarmBoxer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	def xform_center_propagate(self,box,image_name,template,box_size):
		'''
		Centers a box that was generated in a shrunken image by getting the 'real particle' out of the large
		image on disk and doing a ccf with the template - then I just find the peak and use that to center
		@param box a list like [x,y,type,float] - only x and y are used
		@param image_name the name of the image we're operating on
		@param template the template correctly scaled to the have the same angstrom per pixel as the image (named image_name) stored on disk
		@param box_size the size of the box used to center
		Returns the dx and dy parameters, i.e. does not actually alter the box
		'''
	  	global BigImageCache
	  	image = BigImageCache.get_image_directly(image_name)

		xc = box[0]-box_size/2
		yc = box[1]-box_size/2
		r = Region(xc,yc,box_size,box_size)
		particle = image.get_clip(r)
		ccf  = particle.calc_ccf(template)
		trans = ccf.calc_max_location_wrap(particle.get_xsize()/2,particle.get_ysize()/2,0)
		dx = trans[0]
		dy = trans[1]
		return dx,dy

	def try_to_center_ref(self,box_num): #Modified from that in SwarmBoxer
		box = self.target().get_box(box_num)
		box_size = self.target().get_box_size()
		img_filename = self.target().current_file()
		ptcl = box.get_image(img_filename, box_size)
		centered_ptcl = ptcl.process("xform.centeracf")
		dx,dy = self.xform_center_propagate([box.x,box.y],img_filename,centered_ptcl,box_size)
		self.target().move_box(box_num, dx, dy)

		#self.mgph = EMData(micrograph)
		#self.metadata = self.mgph.get_attr_dict()
		#self.particles = {}
		#self._tilesize = 1.5*boxsize
		#self._tile()
		#self._filter()
	
# 	def _tile(self):
# 		# Use e2boxer style. generate coordinates: "{}\t{}\t{}\t{}".format(x,y,boxsize,boxsize)
# 		# use EMBoxList class from emboxerbase?
# 		raise(NotImplementedError)
# 		
# 	def _filter(self):
# 		# need data driven filtration, but this will do for now
# 		for tile in self._tiles:
# 			tile.process_inplace('filter.highpass.gauss',{'cutoff_abs':0.02})
# 			tile.process_inplace('math.sigma',{'value1':15.0,'value2':0.0})
# 			tile.process_inplace('filter.lowpass.gauss',{'cutoff_abs':0.088})
# 			tile.process_inplace('threshold.belowtozero',{'minval':0.098})
# 			tile.process_inplace('morph.majority',{'nmaj':1,'thresh':0.0})
# 			tile.process_inplace('morph.object.density',{'thresh':0.0})
# 			tile.process_inplace('mask.addshells.multilevel',{'nshells':2})
# 			tile.process_inplace('histogram.bin',{'nbins':2})

# 	def _search(self):
# 		"""
# 		Find possible particle coordinates
# 		"""
# 		raise(NotImplementedError)
# 	
# 	def _classify(self):
# 		"""
# 		Performs a classification to discriminate between particles and pure background
# 		"""
# 		raise(NotImplementedError)

class MorphBoxingPanel:
	
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.widget = None
	
	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")
			self.auto_center_checkbox = QtGui.QCheckBox("Auto-center")
			self.clear=QtGui.QPushButton("Clear")
			vbl.addWidget(self.auto_center_checkbox)
			vbl.addWidget(self.clear)
			QtCore.QObject.connect(self.clear, QtCore.SIGNAL("clicked(bool)"), self.clear_clicked)
		return self.widget
	
	def clear_clicked(self,val):
		self.target().clear_all()

class ErasingPanel: # copied for ideas for the morph panel

	def __init__(self,target,erase_radius=128):
		self.busy = True
		self.erase_radius = erase_radius
		self.target = weakref.ref(target)
		self.erase_rad_edit = None
		self.widget = None
		self.busy = False

	def set_erase_radius(self, erase_rad_edit):
		self.busy=True
		self.erase_radius = erase_rad_edit
		if self.erase_rad_edit != None: self.erase_rad_edit.setValue(erase_rad_edit)
		self.busy=False

	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")

			hbl = QtGui.QHBoxLayout()
			hbl.addWidget(QtGui.QLabel("Erase Radius:"))
			from valslider import ValSlider
			self.erase_rad_edit = ValSlider(None,(0.0,1000.0),"")
			self.erase_rad_edit.setValue(int(self.erase_radius))
			self.erase_rad_edit.setEnabled(True)
			hbl.addWidget(self.erase_rad_edit)

			self.unerase = QtGui.QCheckBox("Unerase")
			self.unerase.setChecked(False)

			vbl.addLayout(hbl)
			vbl.addWidget(self.unerase)
			QtCore.QObject.connect(self.erase_rad_edit,QtCore.SIGNAL("sliderReleased"),self.new_erase_radius) #"editingFinished()"
			QtCore.QObject.connect(self.unerase,QtCore.SIGNAL("clicked(bool)"),self.unerase_checked)

		return self.widget

	def new_erase_radius(self, erase_rad_edit):
		if self.busy: return
		self.target().set_erase_radius(erase_rad_edit)

	def unerase_checked(self,val):
		if self.busy: return
		self.target().toggle_unerase(val)


class FilterTester:
	
	"""
	Class to test a specific filter combination on an image.
	"""
	
	def __init__(self,fname):
		self.data = EMData(fname).process('normalize.edgemean')
		self.procedures = []
		self.processed = self.data.copy()
		self._append_procedure()
		self._filter()
	
	def _apply_processor(self,processor,param_dict=None,edgenorm=True,maxmin=True,include=True):
		if param_dict: self.processed.process_inplace(processor,param_dict)
		else: self.processed.process_inplace(processor)
		if maxmin: self.processed.process_inplace('normalize.maxmin')
		if edgenorm: self.processed.process_inplace('normalize.edgemean')
		if include: self._append_procedure()
	
	def _append_procedure(self): self.procedures.append(self.processed.copy())
	
	def display_processing(self): display(self.procedures)
	def display_result(self): display(self.processed)
	def display_initial_final(self): display((self.data,self.processed))
	
	def get_data(self): return self.data
	def get_processing(self): return self.procedures
	def get_processed(self): return self.processed

	def _filter(self):
		self._apply_processor('filter.highpass.gauss',{'cutoff_abs':0.02})
		self._apply_processor('math.sigma',{'value1':15.0,'value2':0.0})
		self._apply_processor('filter.lowpass.gauss',{'cutoff_abs':0.088})
		self._apply_processor('threshold.belowtozero',{'minval':0.098})
		self._apply_processor('morph.majority',{'nmaj':1,'thresh':0.0})
		self._apply_processor('morph.object.density',{'thresh':0.0})
		self._apply_processor('mask.addshells.multilevel',{'nshells':2})
		self._apply_processor('histogram.bin',{'nbins':2})


class MicrographSampler:

	def __init__(self,path,boxsize,nboxes=100):
		fname = os.path.basename(path)
		name,ext = fname.split('.')
		ndir = self.mkpth(os.path.abspath(path)[:-4])
		npath = ndir+'/'+fname
		shutil.copyfile(path,npath)
		self.path = npath
		self.proc_file = self.path[:-4] + '_proc.hdf'
		self.samp_file = self.path[:-4] + '_samp.hdf'
		self.img = EMData(self.path)
		self.hdr = self.img.get_attr_dict()
		self.boxsize = int(boxsize)
		self.boxes = []
		self.nboxes = 0
		self._generated = False
		self._processed = False
		self.gen_random_boxes(nboxes)

	def __repr__(self):
		return "Micrograph(%s,boxsize=%i)"%(self.path,self.boxsize)

	def get_box_at(self,x,y,asarray=False):
		region = Region(x-self.boxsize/2,y-self.boxsize/2,self.boxsize,self.boxsize)
		if asarray: return region.numpy()
		else: return region

	def add_box_at(self,x,y,meta={}):
		region = self.get_box_at(x,y)
		self.add_box(region,meta)

	def add_box(self,region,meta={}):
		if region not in self.boxes:
			self.boxes.append({'region':region,'meta':meta})

	def gen_random_boxes(self,nsamps):
		self.nboxes = self.nboxes + nsamps
		if not self._generated:
			self._random = True
			self._generated = True
		else: print("Appending %i regions for a total of %i"%(nsamps,self.nboxes))
		xlist = np.random.random_integers(int(self.boxsize/2),int(self.hdr['nx']-self.boxsize/2),nsamps)
		ylist = np.random.random_integers(int(self.boxsize/2),int(self.hdr['ny']-self.boxsize/2),nsamps)
		for x,y in zip(xlist,ylist):
			self.add_box_at(x,y)

	def gen_exhaustive_boxes(self,fx=0,fy=0,lx=-1,ly=-1,sx=1,sy=1):
		if not self._generated:
			self._exhaustive = True
			self._generated = True
			if fx == 0: fx = int(self.boxsize/2)
			if fy == 0: fy = int(self.boxsize/2)
			if lx == -1: lx = int(self.hdr['nx']-self.boxsize/2)
			if ly == -1: ly = int(self.hdr['ny']-self.boxsize/2)
			for y in xrange(fy,ly,sy):
				for x in xrange(fx,lx,sx):
					self.add_box_at(x,y)
			self.nboxes = len(self.boxes)
		else: print("Exhaustive regions have already been generated.")

	def display(self):
		display(self.img)

	def show_box_at(self,x,y):
		r = self.get_box_at(x,y)
		display(self.img.get_clip(r))

	def show_boxes(self,infile=None):
		if not infile: infile = self.proc_file
		sp.call(['e2display.py',infile])

	def get_box(self,box_number):
		return self.boxes[box_number]

	def get_box_coords(self,box_number):
		return self.boxes[box_number]['region'].get_origin()[:2]

	def process_boxes(self,processor,params=None,infile=None,outfile=None):
		if not self._generated:
			print("You must first generate regions randomly or exhaustively.")
			return
		if infile == None:
			if self._processed: infile = self.proc_file
			else:
				infile = self.samp_file
				self._processed = True
		# perform additional processing in place
		if not outfile: outfile = self.proc_file
		elif infile == self.proc_file: outfile = infile
		for imgnum in xrange(EMUtil.get_image_count_c(infile)):
			self.process_box(infile,imgnum,processor,params,outfile)

	@staticmethod
	def process_box(infile,imgnum,processor,params=None,outfile=None):
		if not outfile: outfile = infile # process in place
		box = EMData(infile,imgnum)
		if params: box.process_inplace(processor,params)
		else: box.process_inplace(processor)
		box.write_image_c(outfile,imgnum)

	def write_stack(self,outfile=None,imgperfile=100000):
		if not self._generated: print("You must first generate regions randomly or exhaustively.")
		if self._random:
			if not outfile: outfile = self.samp_file
			for i,b in enumerate(self.boxes):
				box = self.img.get_clip(b['region'])
				box.write_image_c(outfile,i)
		elif self._exhaustive:
			path = self.path[:-4]
			os.makedirs(path)
			for i,b in enumerate(self.boxes):
				box = self.img.get_clip(b['region'])
				outname = path + "/r%s.hdf"%(i/imgperfile)
				box.write_image_c(outname,i)

	def get_feature_matrix(self,fname,asdf=False):
		mat = []
		for i in xrange(self.nboxes):
			img = EMData(fname,i)
			mat.append(img.numpy().flatten())
		if asdf:
			try:
				import pandas as pd
				return pd.DataFrame(np.vstack(mat))
			except: return np.vstack(mat)
		else: return np.vstack(mat)

	@staticmethod
	def mkpth(path, stem=''):
		containing = os.path.dirname(os.path.realpath(path))
		contents = os.listdir(containing)
		basename = os.path.basename(path)
		if '_' not in basename: basename += '_00'
		while basename in contents:
			components=basename.split('_')
			if components[-1].isdigit(): components[-1] = str(int(components[-1])+1).zfill(2)
			else: components.append('00')
			basename = '_'.join(components)
		if basename not in contents:
			path = containing + '/' + basename
			os.makedirs(path)
		return path

# mgs = MicrographSampler(path='./micrographs/dh3962.hdf',boxsize=96,nboxes=10000)
# # add labeled background boxes
# mgs.add_box_at(2606,4506,meta={'type':"background",'intensity':"darker"})
# mgs.add_box_at(1724,5078,meta={'type':"background",'intensity':"lighter"})
# mgs.add_box_at(669,503,meta={'type':"background",'intensity':"light"})
# mgs.add_box_at(447,5688,meta={'type':"background",'intensity':"dark"})
# # add labeled particle boxes
# mgs.add_box_at(1938,3916,meta={'type':"particle",'view':"top"})
# mgs.add_box_at(2112,4417,meta={'type':"particle",'view':"side"})
# mgs.add_box_at(3882,3099,meta={'type':"particle",'view':"side"})
# mgs.add_box_at(2644,2445,meta={'type':"particle",'view':"side"})

if __name__ == "__main__":
	main()
