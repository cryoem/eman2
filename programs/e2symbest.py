#!/bin/env python

# $Id$


from EMAN2 import *

from optparse import OptionParser
import os.path


def main():
	progname = os.path.basename(sys.argv[0])
	
	usage = progname + " options inputfile outputfile"
	parser = OptionParser(usage,version=EMANVERSION)
	
	parser.add_option("--nkeep", metavar="N", type="int", help="Number of particles to keep")
	parser.add_option("--sym", metavar="Cn", type="string", help="Symmetry to search for")

	parser.add_option("--mirror", type="string", help="search for particles with mirror symmetry")
	parser.add_option("--rtp", action="store_true", help="make a rotational footprint")

	parser.add_option("--mask", metavar="rad", type="int", help="Mask radius")
	parser.add_option("--imask", metavar="rad", type="int", help="Inside mask radius")

	(options, args) = parser.parse_args()
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

	# sort

	d1 = EMData()
	e1.read_image(inputfile, 0)
	if options.rfp:
		d2 = d1.make_rotational_footprint()
	else:
		d2 = d1

	msk = d2.copy()
	msk.set_parent(None)
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
		imgsort.sort()
		for i in range(options.nkeep):
			s = "%1.1f" % imgsort.score
			d1.write_image(options.mirror, -1, EMUtil::IMAGE_LST, inputfile, imgsort.get_index(i), s)
			
		

if __name__ == "__main__":
    main()
