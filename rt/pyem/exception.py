#!/usr/bin/python

from EMAN2 import *
import os

e = EMData()
e.set_size(100, 100, 1)

try:
	e.filter("GoodMorning")
except RuntimeError, detail:
	print detail

fake_img = os.environ['HOME'] + "/images/fake.mrc"
try:
	e.read_image(fake_img)
except RuntimeError, detail:
	print detail
	
