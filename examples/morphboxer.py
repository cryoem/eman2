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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#

from EMAN2 import *
import sys


def main():
	FilterTester(sys.argv[1]).display_processing()


class EMAutoBoxer:
	
	"""
	Class to automatically pick particles from an cryo-EM micrograph
	"""
	
	def __init__(self,micrograph,boxsize):
		self.mgph = EMData(micrograph)
		self.metadata = self.mgph.get_attr_dict()
		self.particles = {}
		self._tilesize = 1.5*boxsize
		self._tile()
		self._filter()
	
	def _tile(self):
		# important! avoid edges.
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
	
	def autobox(self):
		"""
		Perform automatic particle picking
		"""
		raise(NotImplementedError)
	
	def _classify(self):
		"""
		Performs a classification to discriminate between particles and pure background
		"""
		raise(NotImplementedError)


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


if __name__ == "__main__":
	main()
