#!/bin/env python

from EMAN2 import *
import sys

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
		s = str(i+1) + str(i+1)
		
		if typename == "int":
			#n = TestUtil.get_debug_int(i)
			n = i
		elif typename == "long":
			#n1 = TestUtil.get_debug_int(i)
			n1 = long(i)
			n = long(n1)
		elif typename == "float":
			#n = TestUtil.get_debug_float(i)
			n = float(i)
		elif typename == "string":
			#n = TestUtil.get_debug_string(i)
			n = s
		d[s] = n
		
	return d


def test_point_size():
	nlist = get_list("int")
	flist = get_list("float")

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

def test_map():

	smap = get_dict("string")
	print "before: ", smap
	smap2 = TestUtil.test_map_string(smap)
	print "after: ", smap
	print "result: ", smap2


	lmap = get_dict("long")
	print "before: ", lmap
	lmap2 = TestUtil.test_map_long(lmap)
	print "after: ", lmap
	print "result: ", lmap2

	fmap = get_dict("float")
	print "before: ", fmap
	fmap2 = TestUtil.test_map_float(fmap)
	#print "after: ", fmap
	print "result: ", fmap2
	sys.exit(1)

	imap = get_dict("int")
	print "before: ", imap
	imap2 = TestUtil.test_map_int(imap)
	print "after: ", imap
	print "result: ", imap2
	#assertdict(imap, imap2)

def test_vector():
	nlist = get_list("int")
	flist = get_list("float")
	llist = get_list("long")
	slist = get_list("string")

	nlist2 = TestUtil.test_vector_int(nlist)
	assertlist(nlist, nlist2)

	flist2 = TestUtil.test_vector_float(flist)
	assertlist(flist, flist2)
	
	llist2 = TestUtil.test_vector_long(llist)
	assertlist(llist, llist2)

	slist2 = TestUtil.test_vector_string(slist)
	assertlist(slist, slist2)


test_map()
test_vector()
test_point_size()








