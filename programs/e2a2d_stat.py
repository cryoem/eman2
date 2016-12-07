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
	parser.add_argument("--mode",type=str,default="score",help="Which variable to histogram, score, coverage, alpha, tx, ty. default=score")
	parser.add_argument("--cmp",type=str,help="A comparitor used to generate new quality scores. Will add an additional column with --extract, otherwise replaces score for this run.", default=None)
	parser.add_argument("--multicmp",action="store_true",help="If set will apply a range of different comparators to the aligned particle data and generate a multicolumn text file",default=False)
	parser.add_argument("--ref",type=str,help="A reference image to be used in conjunction with --cmp and --multicmp",default=None)
	parser.add_argument("--extract",action="store_true",help="If set, will convert the .json file to a .txt file suitable for plotting. No histogramming is involved, this is a per-particle conversion",default=False)
	#parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	#parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	#parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls.hdf) containing the aligned subtomograms.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="m2d_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : 
			print "Error, cannot find any m2d_XX folders"
			sys.exit(2)
		options.path = "m2d_{:02d}".format(max(fls))

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print "Cannot find a {}/particle_parms* file".format(options.path)
			sys.exit(2)
		options.iter=max(fls)
		
	angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))

	options.cmp=parsemodopt(options.cmp)
	
	# if cmp is specified, then we add a new "score.alt" element to each image
	if options.cmp[0]!=None:
		try:
			if "," in options.ref : refimg=EMData(options.ref.split(",")[0],int(options.ref.split(",")[1]))
			else: refimg=EMData(options.ref,0)
		except:
			print "ERROR: Unable to read reference image required with --cmp"
			sys.exit(1)
			
		t0=time.time()
		t1=t0
		N=len(angs.keys())
		for n,i in enumerate(angs.keys()):
			if options.verbose==1 and time.time()-t1>1:
				t1=time.time()
				frac=n/float(N)
				try:
					remain=int((time.time()-t0)/frac-(time.time()-t0))	# est remaining time in sec
					print "{:6d}/{:-6d}   time remaining: {}:{:02d}     \r".format(n,N,remain//60,remain%60),
				except:
					print "{:6d}/{:-6d}     \r".format(n,N),
				sys.stdout.flush()
			itm=angs.get(i,True)
			ort=itm["xform.align2d"]
			img=EMData(eval(i)[0],int(eval(i)[1])).process("xform",{"transform":ort})
			sim=img.cmp(options.cmp[0],refimg,options.cmp[1])
			itm["score.alt"]=sim
			angs.setval(i,itm,True)

	# if multicmp is specified, then we make a large multicolumn output file.
	if options.multicmp:
		try:
			if "," in options.ref : refimg=EMData(options.ref.split(",")[0],int(options.ref.split(",")[1]))
			else: refimg=EMData(options.ref,0)
		except:
			print "ERROR: Unable to read reference image required with --multicmp"
			sys.exit(1)
		
		mcmps=[ ("frc",{"minres":80,"maxres":20}),
				("frc",{"minres":20,"maxres":12}),
				("frc",{"minres":12,"maxres":5}),
				("phase",{"minres":80,"maxres":20}),
				("phase",{"minres":20,"maxres":12}),
				("phase",{"minres":12,"maxres":5}),
				("ccc",{}),
				("optsub",{"minres":80,"maxres":12})]
		
		t0=time.time()
		t1=t0
		N=len(angs.keys())
		out=open("{}/particle_multicmp_{:02d}.txt".format(options.path,options.iter),"w")
		for n,i in enumerate(angs.keys()):
			if options.verbose==1 and time.time()-t1>1:
				t1=time.time()
				frac=n/float(N)
				try:
					remain=int((time.time()-t0)/frac-(time.time()-t0))	# est remaining time in sec
					print "{:6d}/{:-6d}   time remaining: {}:{:02d}     \r".format(n,N,remain//60,remain%60),
				except:
					print "{:6d}/{:-6d}     \r".format(n,N),
				sys.stdout.flush()
			itm=angs.get(i,True)
			ort=itm["xform.align2d"]
			img=EMData(eval(i)[0],int(eval(i)[1])).process("xform",{"transform":ort})
			nx=img["nx"]
			apix=img["apix_x"]
			
			out.write("{}\t{}".format(img["ctf"].defocus,sum(img["ctf"].snr[int(2*nx*apix/100):int(2*nx*apix/20)])))
			for c in mcmps:
				sim=img.cmp(c[0],refimg,c[1])
				out.write("\t{}".format(sim))
			imgf=img.process("filter.highpass.gauss",{"cutoff_freq":1.0/80}).process("filter.lowpass.gauss",{"cutoff_freq":1.0/20})
			out.write("\t{}".format(imgf.cmp("ccc",refimg)))
			imgf=img.process("filter.highpass.gauss",{"cutoff_freq":1.0/20}).process("filter.lowpass.gauss",{"cutoff_freq":1.0/12})
			out.write("\t{}".format(imgf.cmp("ccc",refimg)))
			
			out.write("\n")
	
	if options.mode=="score":
		if options.cmp[0]!=None : data=[angs[a]["score.alt"] for a in angs.keys()]
		else : data=[angs[a][options.mode] for a in angs.keys()]
	else:
		data=[angs[a]["xform.align2d"].get_params("2d")[options.mode] for a in angs.keys()]

	if options.extract:
		out=open("{}/particle_parms_{:02d}.txt".format(options.path,options.iter),"w")
		k=angs.keys()
		k.sort(key=lambda i:int(eval(i)[1]))
		for i in k:
			itm=angs[i]
			ort=itm["xform.align2d"].get_params("2d")
			try: altscore=itm["score.alt"]
			except: altscore=0
			out.write("{}\t{}\t{}\t{}\t{}\t{}\t# {};{}\n".format(int(eval(i)[1]),itm["score"],altscore,ort["alpha"],ort["tx"],ort["ty"],eval(i)[1],eval(i)[0]))

	col=array(data)

	# This clears out any serious outliers
	m=col.mean()
	s=col.std()
	col=col[abs(col-m)<s*4.0]

	m=col.mean()
	s=col.std()
	col=col[abs(col-m)<s*4.0]
	savetxt("{}/hist_score.txt".format(options.path),col)
	print "Mean: {}\tSigma: {}".format(m,s)

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
