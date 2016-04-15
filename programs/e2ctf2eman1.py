#!/usr/bin/env python

#
# Author: Benjamin Bammes, 06/03/2008 (bammes@bcm.edu)
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

# e2ctfeman1.py  (Last modified: 04/12/2010 Benjamin Bammes)
# This is a program for determining CTF parameters using the EMAN 1 model


from EMAN2 import *
from math import *
import os


def main( ) :
	"""
	Converts the EMAN2 CTF model to EMAN1 parameters and writes output to ctfparm.txt.
	"""
	
	# Setup command line options parser
	progname = os.path.basename( sys.argv[0] )
	usage = """prog [options]
	
Converts EMAN2 CTF model to the EMAN1 CTF model as well as possible.
Outputs CTF parameters to ctfparm.txt file in the current directory.

To use, run this program in the directory containing your EMAN2DB directory.
It is recommended that you specify a structure factor file (--sf=<file>). 

Always manually check/refine output, since the EMAN1 and EMAN2 CTF
models are not completely compatible."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument( "--first", type = int, help = "First image to include in the plot.", default = -1 )
	parser.add_argument( "--last", type = int, help = "Last image to include in the plot.", default = -1 )
	parser.add_argument( "--sf", type = str, help = "The name of a file containing a structure factor curve.", default = "NULL" )
	parser.add_argument( "--noisemin", type = float, help = "Minimum resolution to examine for determining the noise parameters.", default = -1. )
	parser.add_argument( "--noisemax", type = float, help = "Maximum resolution to examine for determining the noise parameters.", default = -1. )
	parser.add_argument( "--ac", type = float, help = "Set amplitude contrast (percentage, default=10).", default = 10. )
	parser.add_argument( "--bf", type = float, help = "Set constant B-factor (as defined in EMAN1) for all images.", default = 0. )
	parser.add_argument( "--df", action = "store_true", help = "Calculate defocus from entire CCD frame.", default = False )
	parser.add_argument( "--dfmin", type = float, help = "Set minimum possible defocus value (positive is underfocus).", default = 0.5 )
	parser.add_argument( "--dfmax", type = float, help = "Set maximum possible defocus value (positive is underfocus).", default = 5. )
	parser.add_argument( "--dfval", type = float, help = "Set constant defocus for all images (positive is underfocus).", default = 0. )
	parser.add_argument( "--ctfcoverage", action = "store_true", help = "Create a map showing the integrated SNR for the combined data.", default = False )
	parser.add_argument( "--debug", action = "store_true", help = "Show debugging messages.", default = False )
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	# Get and check command line arguments
	( options, args ) = parser.parse_args( )
	pid = E2init(sys.argv, options.ppid)
	
	# Get a list of all images
	ptcls=["particles/"+i for i in os.listdir("particles") if "__" not in i and i[0]!="." and ".hed" not in i ]
	ptcls.sort()
	if len( ptcls ) < 1 :
		parser.error( "No particles found in the particles folder" )
	
	# Convert CTF parameters for all images
	sf = [ ]
	oldds = 0.
	newctflines = [ ]
	intsnr = [ ]
	for i, f in enumerate( ptcls ) :
		if options.first > -1 :
			if i < options.first :
				continue
		if options.last > options.first :
			if i > options.last :
				continue
		print "[%0.1f" % ( 100. * float( i ) / len( ptcls ) ) + "%" + "] Converting %s..." % f
		
		# Load CTF parameters from EMAN2
		if options.debug :
			print "  Loading EMAN2 CTF parameters..."
		( e2ctf, ps, bg ) = loadctf( f )
		ds = e2ctf.dsbg
		
		# Load structure factor file
		if options.debug :
			print "  Loading structure factor file..."
		if ds != oldds :
			sf = loadsf( options.sf, ds, len( bg ) )
			oldds = ds
		
		# Create EMAN1 CTF object to store converted CTF parameters
		if options.debug :
			print "  Creating EMAN1 CTF object..."
		e1ctf = EMAN1Ctf( )
#		e1ctf.defocus = e2ctf.defocus * -1.		# no longer necessary
		e1ctf.voltage = e2ctf.voltage
		e1ctf.cs = e2ctf.cs
		e1ctf.apix = e2ctf.apix
#		e1ctf.ampcont = options.ac / 100.		# no longer necessary
		if options.bf > 0. :
			e1ctf.bfactor = options.bf
#		else :						# no longer necessary
#			e1ctf.bfactor = e2ctf.bfactor / 4.
		
		# Convert noise model
		if options.debug :
			print "  Converting noise model..."
		noise_min = [ 0., -50., 40., 0. ]
		noise_max = [ 30., 0., 140., 10. ]
		( e1ctf.noise0, e1ctf.noise1, e1ctf.noise2, e1ctf.noise3 ) = getnoise( bg, ps, ds, noise_min, noise_max, options.noisemin, options.noisemax )
		
		# Determine defocus
		if options.df :
			if options.debug :
				print "  Calculating defocus..."
			e1ctf.defocus = getdefocus( f, e1ctf, ps, -1. * options.dfmin, -1. * options.dfmax )
		elif options.dfval <> 0. :
			if options.debug :
				print "  Setting defocus..."
			e1ctf.defocus = -1. * options.dfval
		
		# Determine amplitude
		if options.debug :
			print "  Calculating amplitude..."
		e1ctf.amplitude = getamp( ps, sf, e1ctf, e2ctf, ds )
		
		# Add SNR to integrated SNR
		if options.ctfcoverage :
			if options.debug :
				print "  Adding SNR to integrated SNR..."
			e1bg = calcnoise( e1ctf.noise0, e1ctf.noise1, e1ctf.noise2, e1ctf.noise3, ds, bg )
			intsnr = addsnr( intsnr, ps, e1bg )
		
		# Save parameters
		if options.debug :
			print "  Saving result..."
		newctflines.append( [ base_name( f ), str( e1ctf.defocus ) + "," + str( e1ctf.bfactor ) + "," + str( e1ctf.amplitude ) + "," + str( e1ctf.ampcont ) + "," + str( e1ctf.noise0 ) + "," + str( e1ctf.noise1 ) + "," + str( e1ctf.noise2 ) + "," + str( e1ctf.noise3 ) + "," + str( e1ctf.voltage ) + "," + str( e1ctf.cs ) + "," + str( e1ctf.apix ) + "," + str( options.sf ) ] )
	
	# Save results
	if options.debug :
		print "Writing ctfparm.txt file..."
	write_ctfparm( "ctfparm.txt", newctflines )
	
	# Save integrated SNR
	if options.ctfcoverage :
		if options.debug :
			print "Writing ctfcoverage.png file..."
		write_snrmap( "ctfcoverage.png", intsnr, ds )
	
	# End program
	E2end( pid )


def loadsf( filename, ds, pssize ) :
	"""
	Loads an EMAN1 structure factor file or computes a generic structure factor.
	
	Parameters:
	  filename      Filename of EMAN1 structure factor file
	  ds            Size of each pixel in frequency space
	  pssize        Length of the power spectrum curve
	
	Returns:
	  sf            List of structure factor values
	"""
	
	# Load structure factor file if it exists
	empsf = [ ]
	if os.path.exists( filename ) :
		sffile = open( filename, 'r' )
		for line in sffile :
			splitline = map( str.strip, line.split( '\t' ) )
			if len( splitline ) > 1 :
				if ( len( splitline[0] ) > 0 ) and ( len( splitline[1] ) > 0 ) :
					empsf.append( map( float, splitline ) )
		sffile.close( )
		empsf.sort( )
	
	# Interpolate structure factor if loaded from file
	sf = [ ]
	if len( empsf ) > 0 :
		for i in range( pssize ) :
			s = i * ds
			notfound = True
			for j in range( 1, len( empsf ) ) :
				if empsf[j][0] >= s and empsf[j-1][0] < s :
					sf.append( ( ( empsf[j][0] - s ) * empsf[j-1][1] + ( s - empsf[j-1][0] ) * empsf[j][1] ) / ( empsf[j][0] - empsf[j-1][0] ) )
					notfound = False
					break
			if notfound :
				sf.append( genericsf( ds * i ) )
	
	# Get generic structure factor if no structure factor file
	else :
		for i in range( pssize ) :
			sf.append( genericsf( ds * i ) )
	
	# Return result
	return sf


def genericsf( s ) :
	"""
	Returns the value of a generic structure factor curve at the specified
	frequency.
	
	Parameters:
	  s             Frequency
	"""
	
	# Return generic structure factor
	return pow( 10., 3.6717 - 364.58 * s + 15597 * s ** 2 - 4.0678e+05 * s ** 3 + 6.7098e+06 * s ** 4 - 7.0735e+07 * s ** 5 + 4.7839e+08 * s ** 6 - 2.0574e+09 * s ** 7 + 5.4288e+09 * s ** 8 - 8.0065e+09 * s ** 9 + 5.0518e+09 * s ** 10 ) ** 2


def loadctf( filename ) :
	"""
	Loads the CTF, power spectrum, and background curve from EMAN2.
	
	Parameters:
	  filename      Filename of EMAN2 image to load
	
	Returns:
	  e2ctf         EMAN2Ctf object containing EMAN2 CTF parameters
	  ps            1D power spectrum
	  bg            1D background
	"""
	
	js_parms=js_open_dict(info_name(filename))
	ctf=js_parms["ctf"]
	return ctf[:3]

def calcnoise( n0, n1, n2, n3, ds, e2bg ) :
	"""
	Calculates the EMAN1 noise curve.
	
	Parameters:
	  n0            Noise parameter 0
	  n1            Noise parameter 1
	  n2            Noise parameter 2
	  n3            Noise parameter 3
	  ds            Size of each pixel in frequency space
	  e2bg          1D EMAN2 background curve
	
	Returns:
	  e1bg          1D EMAN1 background curve
	"""
	
	# Calcalate noise curve
	e1bg = [ ]
	for i in range( len( e2bg ) ) :
		s = ds * i
		e1bg.append( n2 * exp( -1. * n1 * s - 2.467 * ( n3 * s ) ** 2 - n0 * sqrt( s ) ) )
	
	# Return result
	return e1bg


def getnoise( bg, ps, ds, parm_min, parm_max, s_min, s_max ) :
	"""
	Converts the EMAN2 noise model to the EMAN1 noise model by minimizing
	the RMSD between the EMAN1 noise curve and the EMAN2 background curve.
	
	Parameters:
	  bg            1D background curve
	  ps            1D power spectrum curve
	  ds            Size of each pixel in frequency space
	  parm_min      List of minimum possible parameter values
	  parm_max      List of maximum possible parameter values
	
	Returns:
	  n0            Noise parameter 0
	  n1            Noise parameter 1
	  n2            Noise parameter 2
	  n3            Noise parameter 3
	"""
	
	# Set boundaries for computing noise curve
	if s_min > 0. :
		i0 = min( int( ( 1. / s_min ) / ds ), len( ps ) - 2 )
	else :
		i0 = 1
	if s_max > s_min :
		i1 = min( int( ( 1. / s_max ) / ds ), len( ps ) - 2 )
	else :
		i1 = min( int( 0.20 / ds ), len( ps ) - 2 )
	
	# Set coarseness of search
	c = 10
	
	# Cycle through all possible parameter values coarsely then finely
	parm = [ ( parm_max[0] + parm_min[0] ) / 2., ( parm_max[1] + parm_min[1] ) / 2., ( parm_max[2] + parm_min[2] ) / 2., ( parm_max[3] + parm_min[3] ) / 2. ]  # Result = [ n0, n1, n2, n3 ]
	for dp in range( 3 ) :
		sd = [ ]  # List of [ sd, n0, n1, n2, n3 ]
		for d2 in range( int( c + 1 ) ) :
			n2 = parm_min[2] + d2 * ( parm_max[2] - parm_min[2] ) / float( c )
			for d1 in range( int( c + 1 ) ) :
				n1 = parm_min[1] + d1 * ( parm_max[1] - parm_min[1] ) / float( c )
				for d3 in range( int( c + 1 ) ) :
					n3 = parm_min[3] + d3 * ( parm_max[3] - parm_min[3] ) / float( c )
					for d0 in range( int( c + 1 ) ) :
						n0 = parm_min[0] + d0 * ( parm_max[0] - parm_min[0] ) / float( c )
						validcurve = True
						sdval = 0.
						for i in range( i0, i1 ) :
							s = ds * i
							bgmodel = n2 * exp( -1. * n1 * s - 2.467 * ( n3 * s ) ** 2 - n0 * sqrt( s ) )
							if bgmodel > ps[i] :
								validcurve = False
								break
							else :
								sdval += s ** 2 * ( bg[i] - bgmodel ) ** 2
						if validcurve :
							sd.append( [ sdval, n0, n1, n2, n3 ] )
		if len( sd ) > 0 :
			sd.sort( )
			parm = [ sd[0][1], sd[0][2], sd[0][3], sd[0][4] ]
			dparm = ( parm_max[0] - parm_min[0] ) / ( 1.5 * float( c ) )
			parm_min[0] = parm[0] - dparm
			parm_max[0] = parm[0] + dparm
			dparm = ( parm_max[1] - parm_min[1] ) / ( 1.5 * float( c ) )
			parm_min[1] = parm[1] - dparm
			parm_max[1] = parm[1] + dparm
			dparm = ( parm_max[2] - parm_min[2] ) / ( 1.5 * float( c ) )
			parm_min[2] = parm[2] - dparm
			parm_max[2] = parm[2] + dparm
			dparm = ( parm_max[3] - parm_min[3] ) / ( 1.5 * float( c ) )
			parm_min[3] = parm[3] - dparm
			parm_max[3] = parm[3] + dparm
		else :
			break
	
	# Return result
	return ( parm[0], parm[1], parm[2], parm[3] )


def getdefocus( f, e1ctf, ptclps, parm_min, parm_max ) :
	"""
	Calculates the defocus based on the entire micrograph if available.
	
	Parameters:
	  f             Name of particles file
	  e1ctf         EMAN1Ctf object with correct parameters
	  ptclsize      Particle size (pixels)
	  parm_min      Minimum possible parameter value
	  parm_max      Maximum possible parameter value
	
	Returns:
	  defocus       Defocus
	"""
	
	# Get original defocus
	df = -1. * e1ctf.defocus
	parm_min *= -1.
	parm_max *= -1.
	ptclsize = len( ptclps )
	
	# Try to get entire CCD frame
	ccd_db = db_open_dict( "bdb:raw_data#" + f[:-6] )
	ccd = ccd_db.get( 0 )
	db_close_dict( "bdb:raw_data#" + f[:-6] )
	
	# Return zero defocus if CCD does not exist
	hasccd = True
	if not ccd :
		hasccd = False
	
	# Calculate smoothed 1D power signal for entire frame
	if hasccd :
		ft = ccd.do_fft( )
		ps = ft.calc_radial_dist( int( ft.get_xsize( ) / 2. ), 0., 1., True )
		ds = 0.5 / ( e1ctf.apix * len( ps ) )
		smkernel = [ ]
		smksize = int( len( ps ) / 100. ) + 1
		for i in range( -1 * smksize, smksize + 1 ) :
			smkernel.append( exp( -1. * ( float( i ) / float( smksize ) ) ** 2. ) )
		smps = smooth( ps, smkernel )
		ipm = min( int( 0.01 / ds ), len( ps ) )
		sig = [ ]
		for i in range( len( smps ) ) :
			ipm0 = max( i - ipm, 0 )
			ipm1 = min( i + ipm, len( ps ) )
			sig.append( 100. * smps[i] * len( smps[ipm0:ipm1] ) / sum( smps[ipm0:ipm1] ) )
	
	# Calculate smoothed 1D power signal for particles
	smkernel = [ ]
	smksize = int( len( ptclps ) / 100. ) + 1
	for i in range( -1 * smksize, smksize + 1 ) :
		smkernel.append( exp( -1. * ( float( i ) / float( smksize ) ) ** 2. ) )
	ptclsmps = smooth( ptclps, smkernel )
	ipm = min( int( 0.01 / ds ), len( ptclps ) )
	sig = [ ]
	for i in range( len( ptclsmps ) ) :
		ipm0 = max( i - ipm, 0 )
		ipm1 = min( i + ipm, len( ptclps ) )
		ptclsig.append( 100. * ptclsmps[i] * len( ptclsmps[ipm0:ipm1] ) / sum( ptclsmps[ipm0:ipm1] ) )
	
	# Create model
	ctf = EMAN2Ctf( )
	ctf.from_dict( { "defocus" : df, "voltage" : e1ctf.voltage, "bfactor" : 0.0, "cs" : e1ctf.cs, "ampcont" : e1ctf.ampcont, "apix" : e1ctf.apix, "dsbg" : ds } )
	
	# Cycle through all possible parameter values
	for iteration, courseness in enumerate( [ 0.1, 0.01, 0.002 ] ) :
		
		# Use whole frame for first iteration, then particle images
		ps_use = ps
		smps_use = smps
		sig_use = sig
		if iteration > 0 :
			ps_use = ptclps
			smps_use = ptclsmps
			sig_use = ptclsig
		
		# Set boundaries for computing noise curve
		i0 = max( int( 0.02 / ds ), 3 )
		i1 = min( int( 0.2 / ds ), len( ps_use ) - 3 )
		
		defocusvals = [ ]  # [ zeroavg, df, firstzero, secondzero ]
		ctf.defocus = parm_min
		while ctf.defocus < parm_max :
			ctfmodel = ctf.compute_1d( 2 * len( ps_use ), ds, Ctf.CtfType.CTF_AMP )
			ctfmodel = [ v * v for v in ctfmodel ]
			firstzero = 0
			secondzero = 0
			zerosum = 0.
			zeronum = 0.
			for i in range( i0, i1 ) :
				if ctfmodel[i-1] >= ctfmodel[i] and ctfmodel[i+1] > ctfmodel[i] :
					if firstzero < 1 :
						firstzero = i
					elif secondzero < 1 :
						secondzero = i
					zerosum += sig_use[i]
					zeronum += 1.
			if zeronum > 0. :
				defocusvals.append( [ zerosum / zeronum, ctf.defocus, firstzero, secondzero ] )
			ctf.defocus += courseness
		if len( defocusvals ) > 0 :
			defocusvals.sort( )
			if iteration > 0 :
				df = defocusvals[0][1]
			else :
				for d in defocusvals :
					if d[2] > 0 and d[3] > d[2] :
						drange = d[3] - d[2]
						smbetween = smooth( sig_use, [ 1. ] * int( drange / 2. ) )
						validdefocus = True
						dstart = d[2] + 1
						dstop = d[3]
						dipthresh = ( max( smbetween[dstart:dstop] ) - max( smbetween[dstart], smbetween[dstop] ) ) / 4.
						for i in range( dstart, dstop ) :
							if smbetween[i-1] >= smbetween[i] and smbetween[i+1] > smbetween[i] :
								nearestmax = smbetween[i]
								findnearplus = True
								findnearminus = True
								for j in range( 1, min( i - dstart, dstop - i ) ) :
									if smbetween[i+j-1] < smbetween[i+j] and smbetween[i+j+1] <= smbetween[i+j] and findnearplus :
										nearestmax = max( nearestmax, smbetween[i+j] )
										findnearplus = False
									if smbetween[i-j-1] < smbetween[i-j] and smbetween[i-j+1] <= smbetween[i-j] and findnearminus :
										nearestmax = max( nearestmax, smbetween[i-j] )
										findnearminus = False
								if ( nearestmax - smbetween[i] ) > dipthresh :
									validdefocus = False
									break
						if validdefocus :
							df = d[1]
							break
			parm_min = df - courseness
			parm_max = df + courseness
		else :
			break
	
	# Return result
	return -1. * ds


def getamp( ps, sf, e1ctf, e2ctf, ds ) :
	"""
	Calculates the EMAN1 CTF amplitude by .
	
	Parameters:
	  ps            1D power spectrum curve
	  sf            1D structure factor curve
	  e1ctf         EMAN1Ctf object containing correct noise parameters
	  e2ctf         EMAN2Ctf object
	  ds            Size of each pixel in frequency space
	
	Returns:
	  amp           Amplitude
	"""
	
	# Compute empirical structure factor from power spectrum and CTF model
	amplist = [ ]
	ctf = e2ctf.compute_1d( len( ps ) * 2, ds, Ctf.CtfType.CTF_AMP )
	for i in range( 3, 7 ) :
		e1bg = e1ctf.noise2 * exp( -1. * e1ctf.noise1 * ( ds * i ) - 2.467 * ( e1ctf.noise3 * ( ds * i ) ) ** 2 - e1ctf.noise0 * sqrt( ( ds * i ) ) )
		ampval = sqrt( max( ( ps[i] - e1bg ) / max( sf[i] * ctf[i] ** 2, 0.0001 ), 0. ) )
		if ampval > 0. and ampval < 20. :
			amplist.append( ampval )
	if len( amplist ) > 2 :
		amplist.sort( )
		return sum( amplist[1:-1] ) / len( amplist[1:-1] )
	elif len( amplist ) > 0 :
		return sum( amplist[:] ) / len( amplist[:] )
	else :
		return 0.


def addsnr( intsnr, ps, bg ) :
	"""
	Adds the SNR to the integrated SNR.
	"""
	
	# Create new integrated SNR if necessary
	if len( intsnr ) < len( ps ) :
		intsnr = [ [ 0. for i in range( 100 ) ] for j in range( len( ps ) ) ]
	
	# Calculate SNR
	snr = [ ( ps[i] - bg[i] ) / max( bg[i], 0.0001 ) for i in range( len( ps ) ) ]
	
	# Add SNR to intsnr
	for i in range( len( snr ) ) :
		for j in range( int( round( max( min( snr[i], 0.5 ), 0. ) * 200. ) ) ) :
			intsnr[i][j] += 1
	
	# Return result
	return intsnr


def write_ctfparm( outputfile, newctflines ) :
	"""
	Update the CTF parameters file with new CTF file lines.
	
	Parameters:
	  outputfile    Filename of output ctfparm.txt file
	  newctflines   New CTF file lines to save
	"""
	
	# Get existing CTF file lines
	ctflines = [ ]
	if os.path.exists( outputfile ) :
		ctffile = open( outputfile, 'r' )
		for line in ctffile :
			splitline = map( str.strip, line.split( '\t' ) )
			ctflines.append( splitline )
		ctffile.close( )
	
	# Update lines as necessary
	for newline in newctflines :
		foundline = False
		for line in ctflines :
			if line[0] == newline[0] :
				foundline = True
				line[1] = newline[1]
				break
		if not foundline :
			ctflines.append( [ newline[0], newline[1] ] )
	ctflines.sort( )
	
	# Save new CTF file
	if os.path.exists( outputfile ) :
		os.remove( outputfile )
	ctffile = open( outputfile, 'w' )
	for line in ctflines :
		ctffile.write( '\t'.join( map( str, line ) ) + '\n' )
	ctffile.close( )


def write_snrmap( outputfile, intsnr, ds ) :
	"""
	Update the CTF parameters file with new CTF file lines.
	
	Parameters:
	  outputfile    Filename of output image file
	  intsnr   Integrated SNR
	  ds   Frequency sampling
	"""
	
	# Import plotting libaries
	import matplotlib
	import numpy as np
	import matplotlib.cm as cm
	import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt

	# Setup plotting parameters
	matplotlib.rcParams['xtick.direction'] = 'out'
	matplotlib.rcParams['ytick.direction'] = 'out'

	# Create plotting variables
	x = np.arange( 0., ds * len( intsnr ), ds )
	y = np.arange( 0., 0.5, 0.005 )
	X,Y = np.meshgrid( x, y )
	Z = np.transpose( np.array( intsnr ) )

	# Plot result
	plt.figure( )
	CS = plt.contourf( X, Y, Z, np.arange( 0., np.max( intsnr ) + np.max( intsnr ) / 20., np.max( intsnr ) / 20. ), antialiased = True )
	CB = plt.colorbar( CS, shrink = 0.8, format = '%i' )
	plt.cool( )
	plt.title( 'Image Coverage of Frequency Space' )
	plt.xlabel( 'Frequency (1/A)' )
	plt.ylabel( 'Signal-to-Noise Ratio' )
	plt.savefig( outputfile )



def smooth( seq, kernel ) :
	"""
	Smooths the given sequence using the specified kernel.
	"""
	
	# Create indices range for kernel
	if ( len( kernel ) - 1 ) % 2 :
		kernel.insert( 0, 0. )
	maxj = int( ( len( kernel ) - 1 ) / 2 )
	krange = range( -1 * maxj, maxj + 1 )
	
	# Smooth the sequence
	sm = [ ]
	for i in range( len( seq ) ) :
		num = 0.
		denom = 0.
		for k,j in enumerate( krange ) :
			if i + j > 0 and i + j < len( seq ):
				num += seq[i+j] * kernel[k]
				denom += kernel[k]
		sm.append( num / denom )
	
	# Return result
	return sm


if __name__ == "__main__" :
	main()
