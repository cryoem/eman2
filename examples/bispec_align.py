from EMAN2 import *
from sys import argv

# <file> <im1 #> <im2 #> [rfp=4 [size=32]]

a=EMData.read_images(argv[1],[int(argv[2]),int(argv[3])])

try: rfp=int(argv[4])
except: rfp=4

try: size=int(argv[5])
except: size=32

b=a[0].process("math.bispectrum.slice",{"rfp":rfp,"size":size})
#b.process_inplace("normalize.rows",{"unitlen":1})		# in testing, normalizing rows is definitely a bad thing
#display(b)

c=a[1].process("math.bispectrum.slice",{"rfp":rfp,"size":size})
#c.process_inplace("normalize.rows",{"unitlen":1})

a.append(a[1].process("xform.flip",{"axis":"x"}))
cf=a[2].process("math.bispectrum.slice",{"rfp":rfp,"size":size})
#cf.process_inplace("normalize.rows",{"unitlen":1})

#cc=b.calc_ccfx(c,0,-1,1,0,1)		# z score scaling also seems to work poorly
#ccf=b.calc_ccfx(cf,0,-1,1,0,1)
cc=b.calc_ccfx(c,0,-1,1)
ccf=b.calc_ccfx(cf,0,-1,1)
display((cc,ccf),1)

cc=b.calc_ccfx(c,0,-1,0)
ccf=b.calc_ccfx(cf,0,-1,0)
display((cc,ccf),0,1)

d=a[1].align("rotate_translate_tree",a[0],{"flip":1})
xf=d["xform.align2d"].inverse()
print("rtt\t{:0.2f} {}".format(xf.get_rotation("2d")["alpha"],xf.get_params("2d")["mirror"]))

e=a[1].align("rotate_translate_flip",a[0])
xf=e["xform.align2d"].inverse()
print("rt\t{:0.2f} {}".format(xf.get_rotation("2d")["alpha"],xf.get_params("2d")["mirror"]))

if cc["maximum"]>ccf["maximum"] :
	rot=cc.calc_max_index()*360.0/cc["nx"]
	print("bs\t{:0.2f}".format(rot))
	a[1].rotate(-rot,0,0)
	f=a[1].align("translational",a[0])
	display((a[0],f,d,e),1)
else:
	rot=ccf.calc_max_index()*360.0/cc["nx"]
	print("bs\t{:0.2f} 1".format(rot))
	a[2].rotate(-rot,0,0)
	f=a[2].align("translational",a[0])
	display((a[0],f,d,e),1)
	
