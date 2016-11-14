#!/usr/bin/env python

#
# Author: Liwei Peng, 01/30/2005 (sludtke@bcm.edu)
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
import os.path
import math
import os

IMAGIC_IMG_EXT = "img"
IMAGIC_HED_EXT = "hed"

def get_imagic_filename_pair(infile):
    base = Util.sbasename(infile)
    hedfile = Util.change_filename_ext(base, IMAGIC_HED_EXT)
    imgfile = Util.change_filename_ext(base, IMAGIC_IMG_EXT)
    return (hedfile, imgfile)


def safe_unlink(filename):
	if os.access(filename, os.F_OK):
		os.unlink(filename)

def unlink_data_header_files(filename):
    base = Util.sbasename(filename)
    safe_unlink(base + TestUtil.EMDATA_HEADER_EXT)
    safe_unlink(base + TestUtil.EMDATA_DATA_EXT)

    ext = Util.get_filename_ext(base)
    if ext == IMAGIC_HED_EXT:
        datafile = Util.change_filename_ext(base, IMAGIC_IMG_EXT)
        #safe_unlink(datafile)
    elif ext == IMAGIC_IMG_EXT:
        hedfile = Util.change_filename_ext(base, IMAGIC_HED_EXT)
        #safe_unlink(hedfile)
        

def get_list(typename):
    l = []
    for i in range(3):
        if typename == "int":
            n = TestUtil.get_debug_int(i)
        elif typename == "float":
            n = TestUtil.get_debug_float(i)
        elif typename == "long":
            n1 = TestUtil.get_debug_int(i)
            n = long(n1)
        elif typename == "string":
            n = TestUtil.get_debug_string(i)
        elif typename == "transformarray":
            n = TestUtil.get_debug_transform(i)
        l.append(n)
    return l


def get_dict(typename):
	d = {}
	
	for i in range(3):
		s = TestUtil.get_debug_string(i)
		
		if typename == "int":
			n = TestUtil.get_debug_int(i)
		elif typename == "long":
			n1 = TestUtil.get_debug_int(i)
			n = long(n1)
		elif typename == "float":
			n = TestUtil.get_debug_float(i)
		elif typename == "string":
			n = TestUtil.get_debug_string(i)
#		elif typename == "emobject":
#			n = EMObject(TestUtil.get_debug_float(i))
		d[s] = n
		
	return d

def get_rbase(filename):
    filename = os.path.basename(filename)
    ext_i = filename.rfind(".")
    if ext_i >= 0:
        return filename[0:ext_i]
    else:
        return filename

emdata_counter = 0

def check_emdata(e, progname):
	global emdata_counter
	nx = e.get_xsize()
	ny = e.get_ysize()
	nz = e.get_zsize()

	progname = get_rbase(progname)
	
	if nx > 0 and ny > 0:
		emdata_counter = emdata_counter + 1
		filename = progname + "_" + str(emdata_counter) + ".mrc"
		e.write_image(filename)

		os.unlink(filename)
		
def check_emdata_list(elist, progname):
	for e in elist:
		check_emdata(e, progname)
		
def assertfloat(self, f1, f2):
    delta = 0.001
    diff = f1 - f2
    if math.fabs(diff) > delta:
        self.assertEqual(f1, f2)
