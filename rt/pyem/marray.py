#!/bin/env python
from EMAN2 import *

e  = EMData()
e.read_image(TestUtil.get_debug_image("search.dm3"))
marray = e.get_view()

			
