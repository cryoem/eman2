#!/bin/env python

from EMAN2 import *

num = 100
EMUtil.test_pyem_emobject({"int": num})

fnum = 1.234
EMUtil.test_pyem_emobject({"float": fnum})

lnum = 1000l
EMUtil.test_pyem_emobject({"long": lnum})


fl = [1.0, 2.0, 3.0]
nl = range(5)

EMUtil.test_pyem_emobject({"farray": fl})
EMUtil.test_pyem_emobject({"farray": nl})

e = EMData()
e.set_size(10,20,30)
EMUtil.test_pyem_emobject({"emdata": e})
