#!/bin/env python

from EMAN2 import *

def assertlist(l1, t1):
	assert(l1[0] == t1[0] and l1[1] == t1[1] and l1[2] == t1[2])

def get_nlist():
	nlist = []
	for i in range(3):
		n = TestUtil.get_debug_int(i)
		nlist.append(n)
	return nlist

def get_flist():
	flist = []
	for i in range(3):
		f = TestUtil.get_debug_float(i)
		flist.append(f)
	return flist

nlist = get_nlist()
flist = get_flist()

vec3i = TestUtil.test_Vec3i(nlist)
assertlist(vec3i.as_list(), nlist)

vec3f = TestUtil.test_Vec3f(flist)
assertlist(vec3f.as_list(), flist)

ip1 = TestUtil.test_IntPoint(nlist)
assertlist(ip1, nlist)

fp1 = TestUtil.test_FloatPoint(flist)
assertlist(fp1, flist)

is1 = TestUtil.test_IntSize(nlist)
assertlist(is1, nlist)

fs1 = TestUtil.test_FloatSize(flist)
assertlist(fs1, flist)

mapi1 = TestUtil.test


