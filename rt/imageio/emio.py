#!/bin/env python

from EMAN2 import *

e = EMData()
e.read_image(Util.get_debug_image("search.dm3"))
e.write_image("test.em", 0, EM)
