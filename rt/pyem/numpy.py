#!/bin/env python

from EMAN2 import *
import os

e = EMData()
e.read_image(os.environ['HOME'] + "/images/search.dm3")

a1 = Wrapper.em2numpy(e)
print a1.info()
