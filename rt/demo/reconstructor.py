#!/bin/env python

from EMAN2 import *
import math

e1 = EMData()
e1.read_image(TestUtil.get_debug_image("samesize1.mrc"))

e2 = EMData()
e2.read_image(TestUtil.get_debug_image("samesize2.mrc"))

r = Reconstructors.get("back_projection")
r.set_params({"size":100, "weight":1})
r.setup()
r.insert_slice(e1, Transform(EULER_EMAN, 0,0,0))
r.insert_slice(e2, Transform(EULER_EMAN, math.pi/2,0,0))

result = r.finish()
result.write_image("reconstructor.mrc")
