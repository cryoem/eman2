#!/usr/bin/env python

import os, sys, time

try:
	from optparse import OptionParser
except:
	print "You need Python version 2.3 or later. Check http://www.python.org for more information"
	sys.exit(-1)
try:
	import MLab, Numeric
except:
	print "You need to install Numeric python. Check http://www.pfdubois.com/numpy for more information"
	sys.exit(-1)
try:
	import EMAN2
except:
	print "You need to install EMAN2 python wrapper. Check http://ncmi.bcm.tmc.edu/homes/stevel/EMAN/doc for more infomation"
	sys.exit(-1)
try:
	import Pycluster
except:
	print "You need to install Pycluster package. Check http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/software.htm for more infomation"
	sys.exit(-1)


def main():
	(options,args) =  parse_command_line()
	inputs = args[:-1]
	output = args[-1]
	shrink = options.shrink
	maskrad = options.maskrad
	useamp = options.useamp
	outputamp = options.outputamp
	svalfile = options.svalfile
	eigenimgfile = options.eigenimgfile
	classesnum = options.classes
	svalnum = options.svals
	verbose = options.verbose
	
	d = EMAN2.EMData()
	d.read_image(inputs[0],0,1)
	nx = d.get_xsize()
	ny = d.get_ysize()
	if maskrad<=0: maskrad = ny/2-1

	num = 0
	for f in inputs:
		num+=EMAN2.EMUtil.get_image_count(f)

	mat = Numeric.zeros((num,nx*ny/(shrink*shrink)),typecode='d')
	ptcls = [()] * num
	count=0
	for f in inputs:
		fn, fext = os.path.splitext(f)
		if useamp and outputamp:
			ampfile = "%s.amp%s" % (fn, fext)
			try:
				os.remove(ampfile)
			except:
				pass
		d = EMAN2.EMData()
		for i in range(EMAN2.EMUtil.get_image_count(f)):
			ptcls[count] = (f, i)
			d.read_image(f,i)
			d.process_inplace("eman1.mask.ringmean",{"xc":d.get_xsize()/2.0, "yc":d.get_ysize()/2.0, "outer_radius":maskrad})
			if useamp:
				fft = d.do_fft()
				amp = fft.get_fft_amplitude()
				amp.process_inplace("eman1.math.log")
				img = amp
			else:
				img = d
			if shrink!=1: img.median_shrink(shrink)
			if useamp and outputamp: amp.write_image(ampfile,i)
			
			narray = EMAN2.Wrapper.em2numpy(img)[0]
			mat[count] = narray.flat
			if verbose:
				if useamp:
					print "Image %s %2d (%2d/%d): %dx%d -> %dx%d -> %dx%d" % (f, i, count, num, d.get_xsize(), d.get_ysize(), \
															fft.get_xsize(), fft.get_ysize(), img.get_xsize(), img.get_ysize() )
				else:
					print "Image %s %2d (%2d/%d): %dx%d -> %dx%d" % (f, i, count, num, d.get_xsize(), d.get_ysize(), \
															img.get_xsize(), img.get_ysize() )
			count+=1

	u, w, v = MLab.svd(mat)
	
	if verbose:
		print "mat shape =",mat.shape
		print "u shape =", u.shape
		print "w shape =", w.shape
		print "v shape=", v.shape
		print "w=", w
	
	if svalfile:
		wfile = open(svalfile,"w")
		for i in range(len(w)):
			wfile.write("%d\t%g\n" % (i, w[i]))
		wfile.close()
	if eigenimgfile:
		d = EMAN2.EMData()
		for i in range(len(w)):
			vi = v[i]
			vi = Numeric.reshape(vi,(1,nx/shrink,ny/shrink))
			EMAN2.Wrapper.numpy2em(vi,d)
			d.write_image(eigenimgfile,i)
		
	projval = Numeric.zeros((num,svalnum), typecode='d')
	for i in range(num):
		img = mat[i]
		for k in range(svalnum):
			vk = v[k]
			projval[i][k] = Numeric.dot(img, vk)
		
	#clusterid, centroids, error, nfound = Pycluster.kcluster(whiten, nclusters=3, transpose=0, npass=1, method='a', dist='e')
	clusterid, centroids, error, nfound = Pycluster.kcluster(projval, nclusters=classesnum, transpose=0, npass=1, method='a', dist='e')
	if verbose:
		print "clusterid =", clusterid
		print "centroids=", centroids
	
	outfp = open(output,"w")
	outfp.write("#LST\n")
	for i in range(len(ptcls)):
		outfp.write("%d\t%s\t%d\n" % (ptcls[i][1], ptcls[i][0], clusterid[i]))
	outfp.close()

def parse_command_line():
	usage="%prog <*.img>|<*.hdf> <output lst file> [options]"
	parser = OptionParser(usage=usage, version="%prog v1.0, November 2004. By Wen Jiang <wjiang@bcm.tmc.edu>")
	parser.add_option("--classes",dest="classes",type="int",metavar="<n>",help="number of classes to generate",default=1)
	parser.add_option("--svals",dest="svals",type="int",metavar="<n>",help="number of singular values/vectors to keep for classification",default=1)
	parser.add_option("--shrink",dest="shrink",type="int",metavar="<n>",help="shrink <n> times the power spectrum before classification",default=1)
	parser.add_option("--mask",dest="maskrad",type="int",metavar="<n>",help="mask radius in pixels",default=-1)
	parser.add_option("--useamp",dest="useamp",action="store_true",help="use the FFT amplitude images to classify instead of original images")
	parser.add_option("--outputamp",dest="outputamp",action="store_true",help="to output the computed FFT amplitude images")
	parser.add_option("--svalfile",dest="svalfile",type="string",metavar="<output file name>",help="file name for the singular values output",default="")
	parser.add_option("--eigenimgfile",dest="eigenimgfile",type="string",metavar="<output file name>",help="file name for the eigen image output",default="")
	parser.add_option("-v","--verbose",dest="verbose",action="store_true",help="set to verbose mode")
	
	(options, args)=parser.parse_args()

	if len(sys.argv)<2: 
		parser.print_help()
		sys.exit(-1)

	return (options,args)

if __name__== "__main__":
	main()
