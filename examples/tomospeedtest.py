
from EMAN2 import *
import time
import random

out=file("logfile.txt","a")

for i in range(12,144):
	a=test_image_3d(size=(i,i,i))
	b=a.process("xform",{"transform":Transform({"type":"eman","alt":18.0,"az":240.0,"phi":47.0})})

	t=time.time()
	c=a.align("rotate_translate_3d",b)
	d=a.align("refine_3d",b,{"xform.align3d":c["xform.align3d"]})
	t2=time.time()

	print "%d\t%1.2f\t%s"%(i,t2-t,str(d["xform.align3d"]))

	out.write("%d\t%1.2f\n"%(i,t2-t))
	out.flush()

