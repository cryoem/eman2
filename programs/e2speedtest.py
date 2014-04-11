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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	2111-1307 USA
#
#

import os
import sys
import math
import time
from pprint import pprint
from random import uniform
from numpy import array

from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

This program runs a set of speed tests on the current machine

This program runs a set of speed tests in the current computer. It should
have at least 1 unloaded processor, and minimal i/o occuring when this
is run. It will give a single value which should be generally proportional
to how fast single particle refinement will run on a given computer. It is a 'real world' test,
in that it tests a variety of actual operations performed by a
refinement procedure. Don't compare values given by speed test in
different versions of EMAN, since the underlying routines may be
improved with time."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--slow",action="store_true",help="rtf_slow alignment",default=False)
	parser.add_argument("--best",action="store_true",help="rtf_best alignment",default=False)
	parser.add_argument("--low",action="store_true",help="low level test",default=False)
	parser.add_argument("--size",type=int,help="Size of particles, 192 default for comparisons",default=192)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	SIZE = options.size
	NTT = 105

	print 'This could take a few minutes. Please be patient.'

	print 'Initializing'

#	 pat = EMData();
#	 pat.set_size(SIZE, SIZE, 1)
#	 for i in xrange(-SIZE/2, SIZE/2):
#		 for j in xrange(-SIZE/2, SIZE/2):
#			 value = -3.0 * math.exp(-Util.square(math.fabs(i)+math.fabs(j)) / 10.0) + math.exp(-Util.square( (math.fabs(i)+math.fabs(j))/2.0 )/100.0) * 2.0 if abs(i)<2 else 1.0
#			 pat.set_value_at(i+SIZE/2, j+SIZE/2, value)
#	 pat.update()
	pat=test_image(size=(SIZE,SIZE))
	pat.process_inplace('normalize.circlemean')
	pat.process_inplace("mask.sharp", {"outer_radius":pat.get_xsize()/2.0})

	data = [None for i in xrange(NTT)]

	tmpl=test_image(0,size=(SIZE,SIZE))
	for i in xrange(NTT):
		data[i]=tmpl.process("xform",{"transform":Transform({"type":"2d","alpha":uniform(0,360.0),"tx":uniform(-4.0,4.0),"ty":uniform(-4.0,4.0)})})
		noise=test_image(1,size=(SIZE,SIZE))
		noise.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
		noise.mult(.2)
		data[i].add(noise)
		data[i].process_inplace('normalize.circlemean')
		data[i].process_inplace('mask.sharp', {'outer_radius':data[i].get_xsize()/2})

#		 if i < 5 :
#			 data[i].write_image('speed.hed', i, EMUtil.ImageType.IMAGE_IMAGIC)

	if options.low:
		print 'Low level tests starting. Please note that compiling with optimization may \
invalidate certain tests. Also note that overhead is difficult to compensate for, \
so in most cases it is not dealt with.'
		t1 = time.clock()
		for fj in xrange(500):
			for i in xrange(NTT/2):
				data[i].dot(data[i + NTT / 2])
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 1: %d, %d x %d dot()s in %1.1f sec	%1.1f/sec-> ~%1.2f mflops' % \
			( 500 * NTT / 2, SIZE, SIZE, ti, 500 * NTT / (2.0 * ti), SIZE * SIZE * 4 * 500.0 * NTT / (1000000.0 * ti))

		t1 = time.clock()
		for fj in xrange(500):
			for i in xrange(NTT/2):
				data[1].dot(data[12])
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 1a: %d, %d x %d optimized cached dot()s in %1.1f s %1.1f/s-> ~%1.2f mflops' % \
			(500 * NTT / 2, SIZE, SIZE, ti, 500 * NTT / (2.0 * ti), SIZE * SIZE * 4 * 500.0 * NTT / (1000000.0 * ti))

		t1 = time.clock()
		for j in xrange(500):
			for i in xrange(NTT/2):
				d = {}
				d["keepzero"] = 1
				data[i].cmp("sqeuclidean", data[i + NTT / 2], d)
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 2: %d, %d x %d lcmp()s in %1.1f sec	  %1.1f/sec' % \
			(500 * NTT / 2, SIZE, SIZE, ti, 500 * NTT / (2.0 * ti))

		t1 = time.clock()
		for j in xrange(100):
			for i in xrange(NTT/2):
				data[i].cmp("phase", data[i + NTT / 2], {})
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 3: %d, %d x %d pcmp()s in %1.1f sec	 %1.1f/sec' % \
			(100 * NTT / 2, SIZE, SIZE, ti, 100 * NTT / (2.0 * ti))

		t1 = time.clock()
		for j in xrange(100):
			for i in xrange(NTT/2):
				data[i].cmp("frc", data[i + NTT / 2], {})
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 4: %d, %d x %d fscmp()s in %1.1f sec  %1.1f/sec' % \
			(100 * NTT / 2, SIZE, SIZE, ti, 100 * NTT / (2.0 * ti))

		t1 = time.clock()
		for j in xrange(500):
			for i in xrange(NTT/2):
				data[i].process_inplace("math.absvalue")
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5a: %d, %d x %d fabs in %1.1f sec -> ~%1.2f mfabs/sec' %  \
			(500 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti))

		t1 = time.clock()
		for j in xrange(100):
			for i in xrange(NTT/2):
				data[i].process_inplace("math.sqrt")
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5b: %d, %d x %d sqrts in %1.1f sec -> ~%1.2f msqrt/sec' % \
			(100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti))

		t1 = time.clock()
		dat = data[0].get_2dview()
		for j in xrange(100):
			for i in xrange(NTT/2):
				for m in xrange(SIZE):
					for n in xrange(SIZE):
						math.sqrt(dat[m][n])
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5c: %d, %d x %d sqrts in %1.1f sec -> ~%1.2f msqrt/sec (cached)' % \
			(100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti))
		data[0].update()

		d = data[0].get_2dview()
		t1 = time.clock()
		for j in xrange(100):
			for i in xrange(NTT/2):
				for m in xrange(SIZE):
					for n in xrange(SIZE):
						math.cos(d[m][n])
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5d: %d, %d x %d cos in %1.1f sec -> ~%1.2f mcos/sec (cached)' % \
			(100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti))
		data[0].update()

		d = data[0].get_2dview()
		t1 = time.clock()
		for j in xrange(100):
			for i in xrange(NTT/2):
				for m in xrange(SIZE-1):
					for n in xrange(SIZE-1):
						math.hypot(d[m][n], d[m+1][n+1])
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5e: %d, %d x %d hypot in %1.1f sec -> ~%1.2f mhypot/sec (cached)' % \
			   (100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti))
		data[0].update()

		d = data[0].get_2dview()
		t1 = time.clock()
		for j in xrange(1000):
			for i in xrange(NTT/2):
				for m in xrange(SIZE-1):
					for n in xrange(SIZE-1):
						f = d[m][n] + d[m+1][n+1]
						f = f + f
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5f: %d, %d x %d mult in %1.1f sec -> ~%1.2f mmult/sec (cached)' % \
			   (100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 1000.0 * NTT / (1000000.0 * ti))
		data[0].update()

		d = data[0].get_2dview()
		t1 = time.clock()
		for j in xrange(500):
			for i in xrange(NTT/2):
				for m in xrange(SIZE-1):
					for n in xrange(SIZE-1):
						a = d[m][n] / d[m+1][n+1]
						a = a + a
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5g: %d, %d x %d div in %1.1f sec -> ~%1.2f mdiv/sec (cached)' % \
			   (100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti))
		data[0].update()

		d = data[0].get_2dview()
		t1 = time.clock()
		for j in xrange(500):
			for i in xrange(NTT/2):
				for m in xrange(SIZE):
					for n in xrange(SIZE):
						f = math.fabs(d[m][n])
						f = f + f
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5h: %d, %d x %d fabs in %1.1f sec -> ~%1.2f fabs/sec (cached)' % \
			   (100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti))
		data[0].update()

		d = data[0].get_2dview()
		t1 = time.clock()
		for j in xrange(500):
			for i in xrange(NTT/2):
				for m in xrange(SIZE-1):
					for n in xrange(SIZE-1):
						math.atan2(d[m][n], d[m+1][n+1])
						math.hypot(d[m][n], d[m+1][n+1])
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 5i: %d, %d x %d ri2ap in %1.1f sec -> ~%1.2f ri2ap/sec (cached)' % \
			   (100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti))
		data[0].update()

		t1 = time.clock()
		for i in xrange(NTT*100):
			cp = data[i%NTT].copy()
			d = {}
			d['n'] = 2
			cp.process_inplace('math.meanshrink',d)
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 6:	%1.1f sec %f meanshrink x 2/sec' % (ti, NTT * 100.0 / ti)

		dla = data[0].copy()
		t1 = time.clock()
		for i in xrange(NTT*1000):
			dla.do_fft_inplace()
			dla.update()
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 7:	%1.1f sec %f ffts/sec' % (ti, NTT * 1000 / ti)

		dla = data[0].copy()
		t1 = time.clock()
		for i in xrange(NTT*1000):
			dla.translate(-1,-3,0)
		t2 = time.clock()
		ti = t2-t1
		print 'Baseline 8:	%1.1f sec	%f translates/sec' % (ti, NTT * 1000 / ti)

		return

	tms=[]
	t1 = time.clock()
	for i in xrange(8):
		t11 = time.clock()
		for j in xrange(5, NTT):
			if options.best:
				tmp = data[i].align('rtf_best', data[j], {"flip":None, "maxshift":SIZE/8})
			elif options.slow:
				tmp = data[i].align('rtf_slow', data[j], {"flip":None, "maxshift":SIZE/8})
			else:
				tmp = data[i].align('rotate_translate_flip', data[j], {})
				data[i].del_attr("xform.align2d")
				tmp2 = data[i].align('refine',data[j],{"verbose":0,"xform.align2d":tmp.get_attr("xform.align2d")},"ccc",{})

			if j%5 == 0:
				sys.stdout.write('.')
				sys.stdout.flush()
		tms.append(time.clock()-t11)
		print
	t2 = time.clock()
	tms=array(tms)

	ti = t2 - t1

	if SIZE==256:
		print 'For comparison (Values approximate. Repeated runs will give some variation.)'
		print 'A 2011 MacBook Pro (2.2ghz core-i7) -----------------------------'
		print 'An Intel Xeon E5645 2.4Ghz SF -----------------------------------'
		print 'An Intel Core i5-2500 3.30GHz (depends on turbo) ----------------'
		print 'An Intel Xeon E5-2670 2.6Ghz SF ---------------------------------'
		print 'An Intel Xeon X5675 3.07Ghz SF ----------------------------------'
		print 'An Intel Core i7-3960X 3.3Ghz SF --------------------------------'

	print '\nYour machines speed factor = %1.4f +- %1.4f (%1.4f +- %1.5f sec)\n' % (2.3/tms.mean(),2.3/tms.mean()-2.3/(tms.mean()+tms.std()),tms.mean()/(NTT-5.0),tms.std()/(NTT-5.0))
#	print '\nThis represents %1.2f (RTFAlign+Refine)/sec\n' % (5.0 * (NTT - 5.0) / ti)

if __name__ == "__main__":
	main()
