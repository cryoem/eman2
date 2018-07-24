from EMAN2 import *
from sys import argv

a=EMData.read_images(argv[1],[int(argv[2]),int(argv[3]))

b=a[0].process("math.bispectrum.slice",{"rfp":4,"size":64})
b.process_inplace("normalize.rows",{"unitlen":1})
#display(b)

c=a[1].process("math.bispectrum.slice",{"rfp":4,"size":64})
c.process_inplace("normalize.rows",{"unitlen":1})

cc=b.calc_ccfx(c,0,-1,1)
display(cc)

cc=b.calc_ccfx(c,0,-1,0)

display(cc)
