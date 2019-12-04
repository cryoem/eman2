#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: Steven Ludtke, 12/02/2019 (sludtke42@gmail.com)
# Copyright (c) 2019 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA


from past.utils import old_div
from EMAN2 import *
from numpy import *
from sys import argv,exit
import os
import os.path
import time

def compress_files(jsd,i,fsps,options) :
	"""calls a subprocess to do actual compression"""
	cmd="e2compress.py {}".format(" ".join(fsps))
	for k,v in vars(options).items():
		if not k in ("threads","ppid","positionalargs") and v!=None: cmd+=" --{}={}".format(k,v)
#	print(cmd)
	launch_childprocess(cmd,True)
	jsd.put(i)
	return

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2compress.py [options] <file1> <file2> <file3> 
converts a list of image files to compressed hdf files. If the input file is also HDF, it will overwrite the file, otherwise it will change
the file extension. When read by EMAN2, compressed HDF files will be rescaled to their original (rounded) values, not the mapped integer
values stored in the file. If read using other software it is likely that the mapped integer values will be seen.

Default behavior is to perform 12 bit integer compression, which is sufficient for pretty much any CryoEM image file
or reconstruction. Raw movie frames may need only 2-4 bits and aligned averaged micrographs are likely to be fine with 4-6 bits, so it
is wise to specify the number of bits to use. Specifying 0 bits is a special case which will cause compression of the native floating point
format. In most cases this will result in only 10-30%% compression, whereas most files can be compressed by a factor of 5-10 with no impairment
of results.

Additionally, default behavior is to try and preserve the entire image range in the compressed file, and further, if the input file contains
integer values, to try and produce integer values on output. If there are a number of values which are exactly 0.0, this value will also be
preserved. However, in some cases, the bulk of the data will have very small values with a few large outliers. In such cases, truncating the
outliers will dramatically improve the information retention and/or compression of the output file.

The smaller the number of bits, the faster the compression, and the better the compression ratio. Noise compresses poorly, so eliminating bits
containing pure noise is benificial in multiple ways.

--compresslevel will not impact the quality of the stored images, but will impact compression size and time required. For example, when storing
a typical movie stack of 50 K2 frames using a single thread:
uncompressed   2848 MB   7.1 s
level 0         738      9.3     3.8x compression
level 1         193     15.7    14.7
level 2         184     16.9    15.5
level 3         175     24.7    16.3
level 4         165     21.2    17.3
level 5         158     32.8    18.0  (typical int-compressed tiff)
level 6         152     62.5    18.7
level 7         149     95.4    19.1
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--bits",type=int,help="Bits to retain in the output file, 0 or 2-16. 0 is a flag indicating the native floating point format.",default=12)
	parser.add_argument("--compresslevel",type=int,help="Compression level to use when writing. No impact on image quality, but large impact on speed. Default = 1",default=None)
	parser.add_argument("--range",type=str,help="Specify <minval>,<maxval> representing the largest and smallest values to be saved in the output file. Automatic if unspecified.",default=None)
	parser.add_argument("--sigrange",type=str,help="Specify <minsig>,<maxsig>, eg- 4,4 Number of standard deviations below and above the mean to retain in the output. Default is not to truncate. 4-5 is usually safe.",default=None)
	parser.add_argument("--outpath",type=str,help="Specify a destination folder for the compressed files. This will avoid overwriting existing files.", default=None);
# HDF5 not threadsafe so --threads implemented in an unusual way
	parser.add_argument("--threads",type=int,help="Compression requires significant CPU, this can significantly improve speed",default=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if options.range!=None and options.sigrange!=None:
		print("ERROR: only one of --range and --sigrange may be specified")

	logid=E2init(sys.argv,options.ppid)

	# if threading requested we just call ourselves multiple times with independent processes. Since HDF5 is not threadsafe
	# we cannot do proper internal threading, and must parallelize at the file level. Obnoxious,
	# but the only reasonable approach.
	if options.threads>1:
		import threading
		import queue
		jsd=queue.Queue(0)

		n=-1
		thrds=[(jsd,i,list(args)[i::options.threads],options) for i in range(options.threads)]

		# here we run the threads and save the results, no actual alignment done here
		print(len(thrds)," threads")
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) :
				while (threading.active_count()==options.threads+1) : time.sleep(.1)
				if options.verbose : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
				thrds[thrtolaunch]=threading.Thread(target=compress_files,args=thrds[thrtolaunch])
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(1)

			# no return other than the thread that finished
			while not jsd.empty():
				n=jsd.get()
				thrds[n].join()
				thrds[n]=None
	# this is the actual compression code, which only happens when we have a single thread
	else:
		for f in args:
			fname=os.path.split(f)[-1].rsplit(".",1)[0]
			tmpname="tmp_{}.hdf".format(fname)
			if options.outpath!=None : outpath=os.path.join(options.outpath,fname+".hdf")
			else : outpath=f
			try: os.unlink(tmpname)
			except: pass
			if outpath[-4:]!=".hdf" : outpath=outpath.rsplit(".",1)[0]+".hdf"
			
			N=EMUtil.get_image_count(f)
			if options.verbose : print("Processing {} with {} images".format(f,N))
			t0=time.time()
			for i in range(N):
				im=EMData(f,i)
				im["render_bits"]=options.bits
				if options.compresslevel!=None : im["render_compress_level"]=options.compresslevel
				if options.range!=None:
					rendermin,rendermax=options.range.split(",")
					im["render_min"]=float(rendermin)
					im["render_max"]=float(rendermax)
				elif options.sigrange!=None:
					renderminsig,rendermaxsig=options.sigrange.split(",")
					im["render_min"]=im["mean"]-im["sigma"]*float(renderminsig)
					im["render_max"]=im["mean"]+im["sigma"]*float(rendermaxsig)
				im.write_image(tmpname,i,IMAGE_UNKNOWN,0,None,EM_COMPRESSED)
				
			if options.verbose>1 : print("{:0.1f} s".format(time.time()-t0))
			try: os.unlink(outpath)
			except:pass
			os.rename(tmpname,outpath)
		
	E2end(logid)

		
if __name__ == "__main__":
	main()
