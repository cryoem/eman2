#!/bin/env python

from EMAN2 import *
from testlib import *
import os
import sys

def test_emobject():
	num = TestUtil.get_debug_int(0)
	TestUtil.to_emobject({"int": num})

	fnum = TestUtil.get_debug_float(0)
	TestUtil.to_emobject({"float": fnum})

	lnum = long(num)
	TestUtil.to_emobject({"long": lnum})

	fl = get_list("float")
	TestUtil.to_emobject({"farray": fl})

	e = EMData()
	nx = TestUtil.get_debug_int(0)
	ny = TestUtil.get_debug_int(1)
	nz = TestUtil.get_debug_int(2)
	e.set_size(nx, ny, nz)
	TestUtil.to_emobject({"emdata": e})

	xyd = XYData()
	testfile = "xydata.txt"
	out = open(testfile, "wb")
	for f in fl:
		out.write(str(f) + " " + str(f) + "\n")
	out.close()

	xyd.read_file(testfile)
	TestUtil.to_emobject({"xydata" : xyd})
	os.unlink(testfile)



def test_Dict():
	edict = get_dict("emobject")
	edict2 = TestUtil.test_dict(edict)
	assertdict(edict, edict2)
	

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

	imap = get_dict("int")
	imap2 = TestUtil.test_map_int(imap)
	assertdict(imap, imap2)

	lmap = get_dict("long")
	lmap2 = TestUtil.test_map_long(lmap)
	assertdict(lmap, lmap2)

	fmap = get_dict("float")
	fmap2 = TestUtil.test_map_float(fmap)
	assertdict(fmap, fmap2)

	smap = get_dict("string")
	smap2 = TestUtil.test_map_string(smap)
	assertdict(smap, smap2)

	emobjectmap = get_dict("emobject")
	emobjectmap2 = TestUtil.test_map_emobject(emobjectmap)
	assertdict(emobjectmap, emobjectmap2)


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
	
	e1 = EMData()
	e1.read_image(TestUtil.get_debug_image("samesize1.mrc"))
	e2 = EMData()
	e2.set_size(10, 20, 5)
	e3 = EMData()

	elist = [e1, e2, e3]
	elist2 = TestUtil.test_vector_emdata(elist)
	check_emdata_list(elist, sys.argv[0])
	check_emdata_list(elist2, sys.argv[0])

	p1 = Pixel(1,2,3, 1.1)
	p2 = Pixel(4,5,6, 4.4)
	p3 = Pixel(7,8,9, 5.5)

	plist = [p1,p2,p3]
	plist2 = TestUtil.test_vector_pixel(plist)

	for i in range(len(plist)):
		assert(plist[i] == plist2[i])


test_Dict()
test_map()
test_vector()
test_point_size()
test_emobject()








