#!/bin/env python

from EMAN2 import *
import os.path
import math

def assertfloat(f1, f2):
	diff = f1 - f2
	if math.fabs(diff) > 0.0001:
		assert(f1 == f2)

def assertlist(l1, t1):
	assert(l1[0] == t1[0] and l1[1] == t1[1] and l1[2] == t1[2])

def assertdict(d1, d2):
	assert(len(d1) == len(d2))
	k1 = d1.keys()
	k2 = d2.keys()
	k1.sort()
	k2.sort()
	assertlist(k1, k2)
	for k in k1:
		assert(d1[k] == d2[k])


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
	
def check_emdata_list(elist, progname):
	for e in elist:
		check_emdata(e, progname)
		
