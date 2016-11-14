#!/usr/bin/env python

# This simple program just makes sure that scaling works in the refine aligner

from EMAN2 import *

r=test_image()
i1=test_image()
i1.process_inplace("normalize.edgemean")
i2=test_image()
i2.process_inplace("normalize.edgemean")

x1=Transform({"type":"2d","alpha":4.0,"tx":3.2,"ty":-4.5,"scale":.97})
x2=Transform({"type":"2d","alpha":4.0,"tx":3.2,"ty":-4.5,"scale":1.06})

i1.transform(x1)
i2.transform(x2)

print x1.inverse()
i1a=i1.align("rotate_translate_flip",r)
print "Image 1, initial :",i1a["xform.align2d"]
i1b=i1.align("refine",r,{"stepscale":0.02,"xform.align2d":i1a["xform.align2d"]},"ccc")
print "Image 1, final :",i1b["xform.align2d"]

print x2.inverse()
i2a=i2.align("rotate_translate_flip",r)
print "Image 2, initial :",i2a["xform.align2d"]
i2b=i2.align("refine",r,{"stepscale":0.02,"xform.align2d":i1a["xform.align2d"]},"ccc")
print "Image 2, final :",i2b["xform.align2d"]
