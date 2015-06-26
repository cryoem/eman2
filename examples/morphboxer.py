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
from emboxerbase import EMBoxerModule, EMBoxerTool

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image> <image2>....
	
	Auto Boxing strategy making use of mathematical morphology. It is still in need of some work.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="micrographs",help="List the file to process with e2boxer here.", default="", guitype='filebox', browser="EMBoxesTable(withmodal=True,multiselect=True)",row=0,col=0,rowspan=1,colspan=3,mode="boxing,extraction")
	parser.add_header(name="boxerheader", help='Options below this label are specific to e2boxer', title="### e2boxer options ###", row=1, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--boxsize","-b",type=int,help="Box size in pixels",default=-1, guitype='intbox', row=2, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	(options, args) = parser.parse_args()

	app = EMApp()
	
	# Check . 
	# The EMBoxerModule class (emboxerbase.py:1681) should be the foundation of this autoboxing strategy.
	boxer = MorphBoxer(args,options.boxsize)
	boxer.show_interfaces()
	
	app.execute()
	

	FilterTester(args[1]).display_processing()

class MorphBoxer(EMBoxerModule):
	
	"""
	Class to automatically pick particles from an cryo-EM micrograph
	"""
	
	def __init__(self,micrograph,boxsize):
		super(MorphBoxer,self).__init__(micrograph,boxsize)
		self.mgph = EMData(micrograph)
		self.metadata = self.mgph.get_attr_dict()
		self.particles = {}
		self._tilesize = 1.5*boxsize
		self._tile()
		self._filter()
	
	def _tile(self):
		# important! avoid edges.
		
		# Use e2boxer style. generate coordinates:
		# "{}\t{}\t{}\t{}".format(x,y,boxsize,boxsize)
		
		# use EMBoxList class from emboxerbase?
		
		raise(NotImplementedError)
		
	def _filter(self):
		# need data driven filtration, but this will do for now
		for tile in self._tiles:
			tile.process_inplace('filter.highpass.gauss',{'cutoff_abs':0.02})
			tile.process_inplace('math.sigma',{'value1':15.0,'value2':0.0})
			tile.process_inplace('filter.lowpass.gauss',{'cutoff_abs':0.088})
			tile.process_inplace('threshold.belowtozero',{'minval':0.098})
			tile.process_inplace('morph.majority',{'nmaj':1,'thresh':0.0})
			tile.process_inplace('morph.object.density',{'thresh':0.0})
			tile.process_inplace('mask.addshells.multilevel',{'nshells':2})
			tile.process_inplace('histogram.bin',{'nbins':2})

	def _search(self):
		"""
		Find possible particle coordinates
		"""
		raise(NotImplementedError)
	
	def morph_autobox(self):
		"""
		Perform automatic particle picking
		"""
		raise(NotImplementedError)
	
	def _classify(self):
		"""
		Performs a classification to discriminate between particles and pure background
		"""
		raise(NotImplementedError)


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

# mgs = MicrographSampler(path='./micrographs/dh3962.hdf',# 
#                         boxsize=96,
#                         nboxes=10000)
# 
# # add labeled background boxes
# mgs.add_box_at(2606,4506,meta={'type':"background",'intensity':"darker"})
# mgs.add_box_at(1724,5078,meta={'type':"background",'intensity':"lighter"})
# mgs.add_box_at(669,503,meta={'type':"background",'intensity':"light"})
# mgs.add_box_at(447,5688,meta={'type':"background",'intensity':"dark"})
# 
# # add labeled particle boxes
# mgs.add_box_at(1938,3916,meta={'type':"particle",'view':"top"})
# mgs.add_box_at(2112,4417,meta={'type':"particle",'view':"side"})
# mgs.add_box_at(3882,3099,meta={'type':"particle",'view':"side"})
# mgs.add_box_at(2644,2445,meta={'type':"particle",'view':"side"})

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


class MorphBoxingTool(EMBoxingTool):

	def __init__(self,target):
		super(MorphBoxingTool,self).__init__(target)

if __name__ == "__main__":
	main()
