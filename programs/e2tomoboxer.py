#!/usr/bin/env python

#
# Author: David Woolford 04/16/2009 (woolford@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

from optparse import OptionParser
import sys
import os
from EMAN2 import file_exists, gimme_image_dimensions3D,EMANVERSION,EMData,Region,Transform
from pyemtbx.boxertools import Cache,get_idd_image_entry,set_idd_image_entry

tomo_db_name = "bdb:e2tomoboxercache#"

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
	
Manual particle selection from tomograms. This version is specifically aimed at cubic boxes
for tomographic analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: ref, grid, db, cmd",default=[])

	(options, args) = parser.parse_args()
	print args[0]
	a = TomogramZProjectionImage(args[0])
	a.get_image()
class TomogramZProjectionImage:
	'''
	Sums the tomograph (3D image) along z, stores and returns the result
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the file on disk - should be a 3D tomogram
		@exception RuntimeError raised if the file_name is not a file on disk
		@exception RuntimeError raised if the file on disk is not 3D (nz > 1)
		'''
		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)
		nx,ny,nz = gimme_image_dimensions3D(file_name)
		if nz <= 1: raise RuntimeError("The file must be 3D (nz>1)")
		self.file_name = file_name
		# this will be either None or an EMData object
		self.image = get_idd_image_entry(self.file_name,"tomogram_z_projection",tomo_db_name)
		print self.image
	def get_image_name(self):
		return self.image_name
	
	def get_image(self,use_alternate=False):
		if self.image == None:
			tmp = EMData(self.file_name)
			self.image = tmp.project("standard",Transform())
			self.image.write_image("result.mrc")
#			nx,ny,nz = gimme_image_dimensions3D(self.file_name)
#			self.image = EMData(nx,ny)
#			self.image.to_zero()
#			tmp = EMData()
#			for i in range(nz):
#				print i
#				tmp.read_image(self.file_name,0,False,Region(0,0,i,nx,ny,i+1))
#				tmp.set_size(nx,ny)
#				self.image = self.image + tmp
#			
			# cache to database
			set_idd_image_entry(self.file_name,"tomogram_z_projection",self.image,tomo_db_name)	
		return self.image


class EMBoxerTomoStorage:
	'''
	Initial storage object for a tomographic image
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the tomography file
		@exception RuntimeError raised if the file_name is not a file on disk
		'''
		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)
		self.file_name = file_name
		
if __name__ == "__main__":
	main()
		
		
		