#!/usr/bin/python

from EMAN2 import *

e = EMData()
e.set_size(100, 100, 1)
try:
	e.filter("GoodMorning")
except RuntimeError, detail:
	print type(detail)
except:
	print "nothing"
	
	
