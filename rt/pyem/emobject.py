#!/bin/env python

from EMAN2 import *

num = 100
TestUtil.to_emobject({"int": num})

fnum = 1.234
TestUtil.to_emobject({"float": fnum})

lnum = 1000l
TestUtil.to_emobject({"long": lnum})


fl = [1.0, 2.0, 3.0]
nl = range(5)

TestUtil.to_emobject({"farray": fl})
TestUtil.to_emobject({"farray": nl})

e = EMData()
e.set_size(10,20,30)
TestUtil.to_emobject({"emdata": e})

