#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Benjamin Bammes, 12/08/2008 (bammes@bcm.edu)
# Based on e2ctf.py by: Steven Ludtke, 10/29/2008 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the source
# code in this file under either license. However, note that the complete EMAN2
# and SPARX software packages have some GPL dependencies, so you are responsible
# for compliance with the licenses of these packages if you opt to use BSD
# licensing. The warranty disclaimer below holds in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing author
# citations must be preserved.
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

# e2ctfeman1.py  12/08/2008 Benjamin Bammes
# This is a program for determining CTF parameters using the EMAN 1 model

from EMAN2 import *
from EMAN2db import db_open_dict
from math import *
from optparse import OptionParser
from Simplex import Simplex
import glob, os

# Global variables
# Most variables are global in order to work with the Simplex algorithm
debug = False		  # Set to True to print DEBUG information
name = ""			  # Name of the image file whose CTF parameters are loaded
e1ctf = EMAN1Ctf()	# EMAN1 CTF object
e2ctf = EMAN2Ctf()	# EMAN2 CTF object
bg_1d = []			 # 1D background curve from EMAN2
im_1d = []			 # 1D power spectrum (particle + background) curve from EMAN2
firstmax = 1		   # Index of the first CTC maxima
firstzero = 1		  # Index of the first CTF zero
empsf = []			 # Empirical structure factor (if specified)

def main() :
	"""
	This will take EMAN2 CTF parameters and determine approximate EMAN1 CTF parameters, producing a
	standard EMAN1 ctfparm.txt file. Note the the CTF model in EMAN2 is more flexible than EMAN1, and
	it is thus impossible to make a completely accurate translation. However, it should produce a
	good starting point for further work in EMAN1.
	"""
	global debug, name, e1ctf
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <EMAN2 project directory>
	
Converts EMAN2 CTF model to the EMAN1 CTF model as best as possible.
Outputs CTF parameters to ctfparm.txt file in the current directory.

Always manually check/refine output, since the EMAN1 and EMAN2 CTF
models are not completely compatible.  Most of the time, amplitude
(and possibly envelop) are the only parameters that require manual
adjustment in the EMAN1 ctfit program."""
	parser = OptionParser(usage=usage, version=EMANVERSION)
	parser.add_option("--ac", type="float", help="Set amplitude contrast (percentage, default=10).", default=10.)
	parser.add_option("--sf", type="string", help="The name of a file containing a structure factor curve.", default="NULL")
	parser.add_option("--debug", action="store_true", default=False)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	if len(args) < 1 :
		parser.error("You must specify the EMAN2 project directory.")
	if options.debug :
		debug = True
	
	pid=E2init(sys.argv)
	db_parms = db_open_dict("bdb:e2ctf.parms")
	ptcls = db_parms.keys( )
	
	#ptcls = glob.glob( "particles/EMAN2DB/*_ptcls.bdb" )
	
	if len(ptcls) < 1 :
		parser.error("No particles found in bdb:e2ctf.parms" )
	
	newctflines = [ ]
	loadsf(options.sf)
	for i, f in enumerate(ptcls) :
		print "[%0.1f" % (100. * float(i) / len(ptcls)) + "%" + "] Converting %s..." % (f)
		loadctf(f, options.ac)
		convertnoise()
		convertenvelop()
		newctflines.append([ name, str(e1ctf.defocus) + "," + str(e1ctf.bfactor) + "," + str(e1ctf.amplitude) + "," + str(e1ctf.ampcont) + "," + str(e1ctf.noise0) + "," + str(e1ctf.noise1) + "," + str(e1ctf.noise2) + "," + str(e1ctf.noise3) + "," + str(e1ctf.voltage) + "," + str(e1ctf.cs) + "," + str(e1ctf.apix) + "," + str(options.sf) ])
		if debug :
			print newctflines[ - 1]
	write_ctfparm('ctfparm.txt', newctflines)

	E2end(pid)

def median(vals) :
	"""
	Calculates the median value of a list.
	"""
	sortedvals = sorted(vals)
	if len(sortedvals) > 1 :
		if len(sortedvals) % 2 == 1 :
			return sortedvals[(len(sortedvals) - 1) / 2]
		else:
			lower = sortedvals[len(sortedvals) / 2 - 1]
			upper = sortedvals[len(sortedvals) / 2]
			return float(lower + upper) / 2.
	else :
		return sortedvals[0]

def bgobjfunc(noisevars) :
	"""
	Objective function for estimated the EMAN1 noise model parameters using the Simplex algorithm.
	"""
	global e2ctf, bg_1d, im_1d, firstzero
	obj = 0.
	for i in range(1, len(bg_1d)) :
		s = e2ctf.dsbg * i
		if s > 0.2 :
			break
		e1bg = noisevars[0] * exp(- 1. * noisevars[1] * s - 2.467 * (noisevars[2] * s) ** 2 - noisevars[3] * sqrt(s))
		e2bg = bg_1d[i]
		if (i >= firstzero) and (e1bg <= im_1d[i]) :
			obj += (100000. * s * (e1bg - e2bg)) ** 2
		elif e1bg > im_1d[i] :
			obj += (100000. * s * (e1bg - e2bg)) ** 8
	return obj

def sf(s) :
	"""
	Returns the value of a generic structure factor curve at the specified
	frequency (s).
	"""
	global empsf
	if s < 0.004 :
		return 0.
	if s > 0.2934 :
		s = 0.2934
	if len(empsf) < 10 :
		return pow(10., 3.6717 - 364.58 * s + 15597 * s ** 2 - 4.0678e+05 * s ** 3 + 6.7098e+06 * s ** 4 - 7.0735e+07 * s ** 5 + 4.7839e+08 * s ** 6 - 2.0574e+09 * s ** 7 + 5.4288e+09 * s ** 8 - 8.0065e+09 * s ** 9 + 5.0518e+09 * s ** 10) ** 2
	else :
		besti = 0.
		bestdist = 1.
		for i in range(len(empsf)) :
			if abs(empsf[i][0] - s) < bestdist :
				besti = i
				bestdist = abs(empsf[i][0] - s)
		return empsf[besti][1]

def ctf(s) :
	"""
	Returns the sinusoidal CTF term.
	"""
	global e1ctf
	wl = 12.2639 / sqrt(e1ctf.voltage * 1000 + 0.97845 * e1ctf.voltage ** 2)
	g = 2 * pi * (10000000 * e1ctf.cs * wl ** 3 * s ** 4 / 4. - 10000 * e1ctf.defocus * wl * s ** 2 / 2.)
	return (sqrt(1. - e1ctf.ampcont ** 2) * sin(g) + e1ctf.ampcont * cos(g)) ** 2

def getfirstzero() :
	"""
	Calculates the first CTF zero and the first CTF maxima.
	"""
	global e2ctf, firstzero, firstmax
	ctflist = []
	for i in range(len(im_1d)) :
		ctflist.append(ctf(e2ctf.dsbg * i))
	for i in range(2, len(ctflist) - 1) :
		if (firstmax < 2) and (ctflist[i] > ctflist[i - 1]) and (ctflist[i] >= ctflist[i + 1]) :
			firstmax = i
		if (firstzero < 2) and (ctflist[i] < ctflist[i - 1]) and (ctflist[i] <= ctflist[i + 1]) :
			firstzero = i
			break
	if firstmax < 2 :
		firstmax = firstzero - 2

def loadsf(filename) :
	"""
	Loads an EMAN1 structure factor file.
	"""
	global empsf
	sf = []
	if os.path.exists(filename) :
		sffile = open(filename, 'r')
		for line in sffile :
			splitline = map(str.strip, line.split('\t'))
			if len(splitline) > 1 :
				if (len(splitline[0]) > 0) and (len(splitline[1]) > 0) :
					empsf.append(map(float, splitline))
		sffile.close()

def loadctf(filename, ampcont=10.) :
	"""
	Clears the EMAN1 and EMAN2 CTF models (e1ctf and e2ctf), the loads the
	EMAN2 CTF information from the specified filename.
	"""
	global name, e1ctf, e2ctf, bg_1d, im_1d
	name = get_file_tag(filename)
	db_parms = db_open_dict("bdb:e2ctf.parms")
	[ ctfstring, im_1d, bg_1d, im_2d, bg_2d ] = db_parms[name]
	e2ctf = EMAN2Ctf()
	e2ctf.from_string(ctfstring)
	bg_1d = e2ctf.background
	e1ctf = EMAN1Ctf()
	e1ctf.defocus = e2ctf.defocus * - 1.
	e1ctf.voltage = e2ctf.voltage
	e1ctf.cs = e2ctf.cs
	e1ctf.apix = e2ctf.apix
	e1ctf.ampcont = ampcont / 100.
	getfirstzero()

def convertnoise() :
	"""
	Converts the EMAN2 noise model to the EMAN1 noise model by using the
	Simplex non-linear optimization algorithm to minimize the RMSD between
	the EMAN1 noise curve and the EMAN2 background curve.
	"""
	global bg_1d, e2ctf
	noisevars = [ 8, 10., 0.5, 0. ]
	incr = [ 0.1, 0.1, 0.01, 0.01 ]
	nlopt = Simplex(bgobjfunc, noisevars, incr)
	noisevars = nlopt.minimize(maxiters=1000, monitor=0)[0]
	incr = [ 0.001, 0.001, 0.0001, 0.0001 ]
	nlopt = Simplex(bgobjfunc, noisevars, incr)
	noisevars = nlopt.minimize(maxiters=10000, monitor=0)[0]
	e1ctf.noise0 = noisevars[3]
	e1ctf.noise1 = noisevars[1]
	e1ctf.noise2 = noisevars[0]
	e1ctf.noise3 = noisevars[2]
	
def convertenvelop(scale=1.) :
	"""
	Converts the EMAN2 envelop parameters to the EMAN1 envelop parameter. The
	scale input variable allows for an additional user-specified adjustment
	during the conversion.
	"""
	global e1ctf, e2ctf, bg_1d, im_1d, firstmax
	e1ctf.bfactor = (e2ctf.bfactor / 4.) * scale
	amps = []
	for i in range(2, firstmax) :
		s = e2ctf.dsbg * i
		if (ctf(s) > 0.1) and ((im_1d[i] - bg_1d[i]) > 0) :
			amps.append((im_1d[i] - bg_1d[i]) / (sf(s) * exp(- 2. * e1ctf.bfactor * s ** 2) ** 2 * ctf(s)))
	if len(amps) > 0 :
		e1ctf.amplitude = median(amps)
	else :
		e1ctf.amplitude = 0.1

def write_ctfparm(outputfile, newctflines) :
	"""
	Update the CTF parameters file with newctflines.
	"""
	ctflines = []
	if os.path.exists(outputfile) :
		ctffile = open(outputfile, 'r')
		for line in ctffile :
			splitline = map(str.strip, line.split('\t'))
			ctflines.append(splitline)
		ctffile.close()
	for newline in newctflines :
		foundline = False
		for line in ctflines :
			if line[0] == newline[0] :
				foundline = True
				line[1] = newline[1]
				break
		if not foundline :
			ctflines.append([ newline[0], newline[1] ])
	ctflines.sort()
	if os.path.exists(outputfile) :
		os.remove(outputfile)
	ctffile = open(outputfile, 'w')
	for line in ctflines :
		ctffile.write('\t'.join(map(str, line)) + '\n')
	ctffile.close()


if __name__ == "__main__" :
	main()
