#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
#
#

# $Id$


from EMAN2 import *
import sys
from optparse import OptionParser
import os.path


def main(args):
	sys.argv = args
	progname = os.path.basename(sys.argv[0])
	
	usage = progname + " options inputfile outputfile"
	parser = OptionParser(usage,version=EMANVERSION)
	
	parser.add_option("--nkeep", metavar="N", type="int", help="Number of particles to keep")
	parser.add_option("--sym", metavar="Cn", type="string", help="Symmetry to search for")

	parser.add_option("--mirror", metavar="outputfile", type="string", help="search for particles with mirror symmetry and write them out to outputfile.")
	parser.add_option("--rtp", action="store_true", help="make a rotational footprint")

	parser.add_option("--mask", metavar="rad", type="int", help="Mask radius")
	parser.add_option("--imask", metavar="rad", type="int", help="Inside mask radius")
	
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if len(args) != 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
        sys.exit(1)
		
	
	inputfile = args[0]
	outputfile = args[1]

	if sym[0] != 'c' and sym[0] != 'C':
		print "Error: invalid  symmetry ",  sym
		sys.exit(1)

	csym = int(sym[1:])

	nimg = EMUtil.get_image_count(inputfile)
	if nimg < 2:
		print "Error: Not enough images to sort"
		sys.exit(1)


	d1 = EMData()
	e1.read_image(inputfile, 0)
	if options.rfp:
		d2 = d1.make_rotational_footprint()
	else:
		d2 = d1

	msk = d2.copy()
	#msk.set_parent(None)  #no parent for EMAData any more 
	msk.to_one()
	#msk->applyMask(mask,4);
	#msk->applyMask(imask,5);

	imgsort = ImageSort(nimg)
	imgsortm = ImageSort(nimg)

	angle_step = 2.0 * math.pi / csym
	angle_max = 2.0 - 1.0 / csym * math.pi
	
	for i in range(n):
		d1.read_image(inputfile, i)

		d1.normalize()
		if options.rfp:
			d2 = d1.make_rotational_footprint()
		else:
			d2 = d1

		tmpimg = d2.copy()

		score = 0
		scorem = 0
		
		angle = angle_step
		while angle < angle_max:
			#tmpimg.setRAlign(angle,0,0);
			tmpimg.rotate_translate()
			tmpimg.mult(msk)
			#cmpscore=2.0-tmp->lcmp(d2)
			score += cmpscore

			angle += angle_step

		if options.mirror:
			tmpimg.set_parent(None)
			#tmpimg.vFlip()
			#tmpimg.rotAlign(d2)
			tmpimg.rotate_translate()
			#scorem = -score + (2.0 - tmpimg->lcmp(d2))
			imgsortm.set(i, scorem)
			print i, score/(csym-1), scorem/(csym-1)
		else:
			imgsort.set(i, score)
			print i, score/(csym-1)

	if options.mirror:
		imgsortm.sort()
		for i in range(options.nkeep):
			score1 = "%1.1f" % imgsortm.get_score(i)
			d1.write_lst(options.mirror, -1, inputfile, imgsortm.get_index(i), score1)
	
	imgsort.sort()
	for i in range(options.nkeep):
		score1 = "%1.1f" % imgsort.get_score(i)
		d1.write_lst(outputfile, -1, inputfile, imgsort.get_index(i), score1)


if __name__ == "__main__":
    main(sys.argv)
