#!/bin/env python

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
        safe_unlink(datafile)
    elif ext == IMAGIC_IMG_EXT:
        hedfile = Util.change_filename_ext(base, IMAGIC_HED_EXT)
        safe_unlink(hedfile)
        

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
		elif typename == "emobject":
			n = EMObject(TestUtil.get_debug_float(i))
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
		
