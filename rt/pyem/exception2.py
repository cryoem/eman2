#!/usr/bin/python

from EMAN2 import *
import os
from pyemtbx.exceptions import *

e = EMData()
e.set_size(100, 100, 1)

try:
	e.read_image("notexistingfile.mrc")
	
except RuntimeError, runtime_err:
	err_type = exception_type(runtime_err)
	
	if err_type == "ImageFormatException":
		print "eman2", err_type
	else:
		print "RuntimeError  ", runtime_err
else:
	print "Expected a RuntimeError!"
	
