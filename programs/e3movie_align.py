#!/usr/bin/env python
#
# Author: Steven Ludtke  06/07/2023
# Copyright (c) 2023- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

from EMAN3 import *
from EMAN3tensor import *
import random

def main():

	usage="""e3movie.py <movie stack> ...

At the moment this program provides only an option for estimating the gain image from a large set of counting-mode images. Eventually this will include movie alignment.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--est_gain", type=str,help="specify output file for gain image. Estimates a gain image when given a set of many movies via hierarchical median estimation", default=None)
	parser.add_argument("--ccftest", type=str,help="specify output file for gain image. Estimates a gain image when given a set of many movies via hierarchical median estimation", default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	pid=E3init(argv)
	nmov=len(args)

	if options.est_gain is not None:

		sig=EMData()
		avgr=Averagers.get("mean",{"sigma":sig})
		for i in range(nmov):
			if (i%10==0): E3progress(pid,i/nmov)
			nimg=EMUtil.get_image_count(args[i])
			print(f"{args[i]}: {nimg}")
			if (nimg>500) : raise Exception("Can't deal with movies with >500 frames at present")
			for j in range(0,nimg,50):
				imgs=EMData.read_images(f"{args[i]}:{j}:{j+50}")
				for img in imgs: avgr.add_image(img)

		avg=avgr.finish()
		avg.write_image("average.hdf",0)
		sig.write_image("average.hdf",1)
		avg2=avg.copy()
		avg2.mult(sig.process("math.reciprocal"))
		avg2.write_image("average.hdf",2)

	if options.ccftest is not None:
		avg=EMData("average.hdf",0)
		avg.div(avg["mean"])
		#avg.add(-avg["mean"])
		nimg=EMUtil.get_image_count(args[0])
		imgs=EMDataStack2D(EMData.read_images(f"{args[i]}:0:{min(50,nimg)}"))
		for im in imgs:
			im.div(avg)
		ffts=imgs.do_fft()
		ccfs=imgs.calc_ccf(imgs[len(imgs)//2])
		ccfsr=ccfs.do_ift()

		ccfsr.write_images("ccfs.hdf")


	E3end(argv)

if __name__ == '__main__':
	main()
