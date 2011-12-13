#!/usr/bin/env python

from EMAN2    import *
from sparx    import *
from random   import random, seed, randint
from sys      import argv, exit
from optparse import OptionParser
from EMAN2db import db_check_dict

dbkey  = sys.argv[1]
print dbkey
gbdbname = 'bdb:e2boxercache#gauss_box_DB'

gbdb = db_open_dict(gbdbname)
gbdb[dbkey] = {'gauss_width':1.0,'pixel_input':5.2,'pixel_output':5.2,'thr_low':4.0,'thr_hi':60.0,"invert_contrast":False,"use_variance":True,"boxsize":64}
