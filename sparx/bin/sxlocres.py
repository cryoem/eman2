#!/usr/bin/env python
#
# Author: Pawel A.Penczek and Edward H. Egelman 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
# Copyright (c) 2008-Forever The University of Virginia
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#
import global_def
from   global_def import *
from EMAN2 import *
from sparx import *

def main():
	import os
	import sys
	from optparse import OptionParser
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + """ firstvolume  secondvolume maskfile outputfile --wn --step --cutoff  --radius  --fsc

	    Compute local resolution in real space within are outlined by the maskfile and within regions wnxwnxwn

	"""
	parser = OptionParser(usage,version=SPARXVERSION)
	
	parser.add_option("--wn",		type="int",		default=7, 			help="Size of window wiyhin which local real-space FSC is computed")
	parser.add_option("--step",     type="float",	default= 1.0,       help="Shell step in Fourier size in pixels (default 1.0)")   
	parser.add_option("--cutoff",   type="float",	default= 0.5,       help="resolution cut-off for FSC (default 0.5)")
	parser.add_option("--radius",	type="int",		default=-1, 		help="if there is no maskfile, sphere with r=radius will be used, by default the radius is nx/2-wn")
	parser.add_option("--fsc",      type="string",	default= None,      help="overall FSC curve (migh be truncated) (default no curve)")

	(options, args) = parser.parse_args(arglist[1:])
	
	if len(args) <3 or len(args) > 4:
		print "See usage " + usage
		sys.exit()
	
	vi = get_im(args[0])
	ui = get_im(args[1])

	nn = vi.get_xsize()
	wn = int(options.wn)

	kern = model_blank(wn,wn,wn,1.0)

	
	if len(args) == 3:
		m = model_circle((nn-wn)//2,nn,nn,nn)
		outvol = args[2]
		
	elif len(args) == 4:
		m = get_im(args[2])
		outvol = args[3]

	mc = model_blank(nn,nn,nn,1.0)-m

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()


	vf = fft(vi)
	uf = fft(ui)

	lp = int(nn/2/options.step+0.5)
	step = 0.5/lp

	freqvol = model_blank(nn,nn,nn)
	resolut = []
	for i in xrange(lp):
		print lp,i
		fl = step*i
		fh = fl+step
		v = fft(filt_tophatb( vf, fl, fh))
		u = fft(filt_tophatb( uf, fl, fh))
		tmp1 = Util.muln_img(v,v)
		tmp2 = Util.muln_img(u,u)

		do = Util.infomask(square_root(Util.muln_img(tmp1,tmp2)),m,True)[0]


		tmp3 = Util.muln_img(u,v)
		dp = Util.infomask(tmp3,m,True)[0]

		resolut.append([i,(fl+fh)/2.0, dp/do])


		tmp1 = rsconvolution(tmp1, kern)
		tmp2 = rsconvolution(tmp2, kern)
		tmp3 = rsconvolution(tmp3, kern)

		Util.mul_img(tmp1,tmp2)

		tmp1 = square_root(tmp1)

		Util.mul_img(tmp1,m)
		Util.add_img(tmp1,mc)


		Util.mul_img(tmp3,m)
		Util.add_img(tmp3,mc)

		Util.div_img(tmp3,tmp1)

		Util.mul_img(tmp3,m)
		freq=(fl+fh)/2.0
		bailout = True
		for x in xrange(nn):
			for y in xrange(nn):
				for z in xrange(nn):
					if(m.get_value_at(x,y,z) > 0.5):
						if(freqvol.get_value_at(x,y,z) == 0.0):
							if(tmp3.get_value_at(x,y,z) <0.5):
								freqvol.set_value_at(x,y,z,freq)
							else:
								bailout = False
		if(bailout):  break

	freqvol.write_image(outvol)
	if(options.fsc != None): write_text_row(resolut, options.fsc)

if __name__ == "__main__":
	main()
