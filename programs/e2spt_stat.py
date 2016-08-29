#!/usr/bin/env python

# This is a simple example showing how to generate a histogram from a text file
# specify the filename and column number with an optional number of bins, column number 0 indexed
# Note that outliers are filtered out (>sigma*4 twice)
from EMAN2 import *
from numpy import *
from sys import argv,exit
try:
	import matplotlib
#	matplotlib.use("AGG")
	import matplotlib.pyplot as plt
	pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]
except:
	print "ERROR: Matplotlib not available, cannot generate histogram"
	sys.exit(1)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_stat.py [options] 
Note that this program is not part of the original e2spt hierarchy, but is part of an experimental refactoring.

This program will look in an spt_XX folder at particle_parms_xx.json and show a quality histogram of the particles as compared to the alignment reference. Smaller values are better.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path",type=str,default=None,help="Path to a folder containing current results (default = highest spt_XX)")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = auto",default=0)
	parser.add_argument("--bins",type=int,help="Number of bins to use in the histogram",default=100)
	parser.add_argument("--gui",action="store_true",help="If set will open an interactive plot with the results",default=False)
	parser.add_argument("--mode",type=str,default="score",help="Which variable to histogram, score, coverage, alt, az, phi, dx, dy, dz. default=score")
	parser.add_argument("--extract",action="store_true",help="If set, will convert the .json file to a .txt file suitable for plotting. No histogramming is involved, this is a per-particle conversion",default=False)
	#parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	#parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	#parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls.hdf) containing the aligned subtomograms.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="spt_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : 
			print "Error, cannot find any spt_XX folders"
			sys.exit(2)
		options.path = "spt_{:02d}".format(max(fls))

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print "Cannot find a {}/particle_parms* file".format(options.path)
			sys.exit(2)
		options.iter=max(fls)
		

	angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
	if options.mode in ("score","coverage"):
		data=[angs[a][options.mode] for a in angs.keys()]
	else:
		data=[angs[a]["xform.align3d"].get_params("eman")[options.mode] for a in angs.keys()]
	

	if options.extract:
		out=open("{}/particle_parms_{:02d}.txt".format(options.path,options.iter),"w")
		k=angs.keys()
		k.sort(key=lambda i:int(eval(i)[1]))
		for i in k:
			itm=angs[i]
			ort=itm["xform.align3d"].get_params("eman")
			out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(int(eval(i)[1]),itm["score"],itm["coverage"],ort["az"],ort["alt"],ort["phi"]))

	col=array(data)

	# This clears out any serious outliers
	m=col.mean()
	s=col.std()
	col=col[abs(col-m)<s*4.0]

	m=col.mean()
	s=col.std()
	col=col[abs(col-m)<s*4.0]

	if options.verbose:
		lz=len(col[col<0])
		gz=len(col[col>0])
		print "%1.2f (%d) less than zero"%(float(lz)/(lz+gz),lz)
		print "%1.2f (%d) less than zero"%(float(gz)/(lz+gz),gz)

	his=histogram(col,options.bins)

	#out=file("hist/"+argv[1],"w")
	#for i in xrange(len(his[0])): out.write("%f\t%f\n"%(his[1][i],his[0][i]))

	fig = plt.figure()
	ax = plt.axes([.15,.15,.8,.8])
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	#plt.title("Convergence plot (not resolution)")
	plt.title("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
	plt.xlabel("Score",fontsize=24)
	plt.ylabel("Number of Particles",fontsize=24)

	plt.bar(his[1][:-1],his[0],his[1][1]-his[1][0])
	if options.gui: plt.show()
	plt.savefig("{}/hist_score.pdf".format(options.path))

if __name__ == "__main__":
	main()
