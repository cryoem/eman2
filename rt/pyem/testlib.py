#!/bin/env python

from EMAN2 import *


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

emdata_counter = 0

def check_emdata(e):
	global emdata_counter
	nx = e.get_xsize()
	ny = e.get_ysize()
	nz = e.get_zsize()
	
	if nx > 0 and ny > 0:
		emdata_counter = emdata_counter + 1
		filename = "check_emdata_" + str(emdata_counter) + ".mrc"
		e.write_image(filename)
	
def check_emdata_list(elist):
	for e in elist:
		check_emdata(e)
		

