#!/usr/bin/python

from EMAN2 import *

e = EMData()
e.set_size(100, 100, 1)

try:
	e.filter("GoodMorning")
except RuntimeError, detail:
	print detail

fake_img = Util.get_debug_image("fake.mrc")
try:
	e.read_image(fake_img)
except RuntimeError, detail:
	print detail
	
