# This is a simple example showing how to generate a histogram from a text file
# specify the filename and column number
from EMAN2 import *
from numpy import *
from sys import argv

data=loadtxt(argv[1]).transpose()
coln=int(argv[2])
col=data[coln]

m=col.mean()
s=col.std()
col=col[abs(col-m)<s*4.0]

m=col.mean()
s=col.std()
col=col[abs(col-m)<s*4.0]

his=histogram(col,100)

try: os.mkdir("hist")
except: pass

out=file("hist/"+argv[1],"w")
for i in xrange(len(his[0])): out.write("%f\t%f\n"%(his[1][i],his[0][i]))


