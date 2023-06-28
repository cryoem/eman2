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

from EMAN2 import *
import random

def main():

	usage="""e3movie_align.py <movie stack> ...

At the moment this program provides only an option for estimating the gain image from a large set of counting-mode images. Eventually this will include movie alignment.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--est_gain", type=str,help="specify output file for gain image. Estimates a gain image when given a set of many movies via hierarchical median estimation", default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	nmov=len(args)

	if options.est_gain is not None:

		med1=Averagers.get("mean")
		for i in range(400):
			med2=Averagers.get("median")
			for j in range(25):
				msel=random.randint(0,nmov-1)	# pick a random movie
				nimg=EMUtil.get_image_count(args[msel])
				isel=random.randint(2,nimg-2)	# skip the first 2 frames since they may be unusual
				img=EMData(args[msel],isel)
				if options.verbose>1 : print(i,j,msel,isel)
				med2.add_image(img)
			img2=med2.finish()
			med1.add_image(img2)
		final=med1.finish()
		#final.mult(1.0/final["mean"])

		final.write_image(options.est_gain,0)

if __name__ == '__main__':
	main()
