# This is a simple example showing how to generate a histogram from a text file
# specify the filename and column number with an optional number of bins, column number 0 indexed
# Note that outliers are filtered out (>sigma*4 twice)
from EMAN2 import *
from numpy import *
from sys import argv,exit
try:
	import matplotlib
	matplotlib.use("AGG")
	import matplotlib.pyplot as plt
	pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]
except:
	print "Matplotlib not available, some output will not be generated"

if len(argv)<2 :
	print "usage:\nhistogram.py <txtfile> [col#=0] [nbins=100]"
	sys.exit(1)

data=loadtxt(argv[1])
if data.ndim!=1:
	data=data.transpose()
	try: coln=int(argv[2])
	except: coln=0
	col=data[coln]
else: col=data

m=col.mean()
s=col.std()
col=col[abs(col-m)<s*4.0]

m=col.mean()
s=col.std()
col=col[abs(col-m)<s*4.0]

lz=len(col[col<0])
gz=len(col[col>0])
print argv[1]
print "%1.2f (%d) less than zero"%(float(lz)/(lz+gz),lz)
print "%1.2f (%d) less than zero"%(float(gz)/(lz+gz),gz)

try: his=histogram(col,int(argv[3]))
except: his=histogram(col,100)

try: os.mkdir("hist")
except: pass

out=file("hist/"+argv[1],"w")
for i in xrange(len(his[0])): out.write("%f\t%f\n"%(his[1][i],his[0][i]))

fig = plt.figure()
ax = plt.axes([.15,.15,.8,.8])
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)

#plt.title("Convergence plot (not resolution)")
plt.xlabel("Similarity Difference",fontsize=24)
plt.ylabel("Number of Particles",fontsize=24)

plt.bar(his[1][:-1],his[0],his[1][1]-his[1][0])
plt.savefig("hist/"+argv[1]+".pdf")
