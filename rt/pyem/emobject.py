#!/bin/env python

from EMAN2 import *
from testlib import *

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

