#!/usr/bin/env python

#
# Author: Steven Ludtke, 06/27/2007 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

# e2resolution.py  06/27/2007	Steven Ludtke
# This program computes a similarity matrix between two sets of images

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <volume> <mask> <output>
	This will compute an FSC resolution curve from a single structure by
	comparing the power spectrum of the noise region to the power spectrum
	of the reconstructed particle. The 'mask' input should mask out the particle,
	and should be a 'soft' mask. Sharp edges on the mask may ruin the results.
	proc3d (EMAN1) with the automask2 option can produce appropriate masks. The
	3rd parameter (in automask2) should be about 10% of the box size. This will not
	work with particles that have been tightly masked already. If doing an EMAN1
	reconstruction with the amask= option, you must also use the refmaskali option
	(in the refine command). Note that the resultant curve is guaranteed to go
	to zero at near Nyquist because it assumes that the SNR is zero near Nyquist.
	If your data is undersampled, the resulting curve may also be inaccurate.
	Models should also not be heavily low-pass filtered
	"""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--apix", "-A", type="float", help="A/voxel", default=-1.0)
	#parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	#parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
	
	parser.add_option("--path", default=None, type="string",help="The name the e2refine directory that contains the reconstruction data. If specified will place curves generated in bdb:path#convergence.results")
	
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input and output files required")
	
	E2n=E2init(sys.argv)
	
	#options.align=parsemodopt(options.align)

	print "read models"
	data=EMData(args[0],0)
	mask=EMData(args[1],0)

#	mask.to_one()
#	mask.process_inplace("mask.gaussian",{"outer_radius":mask.get_xsize()/6,"inner_radius":mask.get_xsize()/6})
	mask.write_image("msk.mrc",0)

	# inverted mask
	maski=mask.copy()
	maski*=-1.0
	maski+=1.0
	maski.process_inplace("mask.gaussian",{"outer_radius":mask.get_xsize()/6,"inner_radius":mask.get_xsize()/3})
	maski.write_image("msk2.mrc",0)

	noise=data.copy()
	noise*=maski

	data*=mask

	print "compute FFT"
	dataf=data.do_fft()
	noisef=noise.do_fft()

	print "compute power 1"
	datapow=dataf.calc_radial_dist(dataf.get_ysize()/2-1,1,1,1)
	print "compute power 2"
	noisepow=noisef.calc_radial_dist(noisef.get_ysize()/2-1,1,1,1)

	x=range(1,len(datapow)+1)
	if options.apix>0:
		x=[i/(len(datapow)*options.apix*2.0) for i in x]
	else:
		x=[i/(len(datapow)*data["apix_x"]*2.0) for i in x]

	# normalize noise near Nyquist
	s=0
	sn=0
	for i in range(int(len(noisepow)*.9),len(noisepow)-1):
		if datapow[i]<datapow[i+1] or noisepow[i]<noisepow[i+1] : continue
		s+=datapow[i]/noisepow[i]
		sn+=1.0
	if sn==0 :
		print "Warning, strange normalization"
		s=datapow[int(len(noisepow)*.9)]/noisepow[int(len(noisepow)*.9)]
	else: s/=sn

	noisepow=[i*s for i in noisepow]
#
#	# normalize based on volume
#	datapow=[v/mask["mean"] for v in datapow]
#	noisepow=[v/maski["mean"] for v in noisepow]

	# compute signal to noise ratio
	snr=[]
	for i in range(len(datapow)):
		try: snr.append((datapow[i]-noisepow[i])/noisepow[i])
		except: snr.append(0)
	
	# convert to FSC
	fsc=[i/(2.0+i) for i in snr]

	out=file(args[2],"w")
	for i in range(len(fsc)): out.write("%f\t%f\n"%(x[i],fsc[i]))
	out.close()
	
	out=file(args[2]+".dat","w")
	for i in range(len(fsc)): out.write("%f\t%f\n"%(x[i],datapow[i]))
	out.close()

	out=file(args[2]+".noi","w")
	for i in range(len(noisepow)): out.write("%f\t%f\n"%(x[i],noisepow[i]))
	out.close()
	
	if options.path != None:
		s = get_file_tag(args[0])
		db = db_open_dict("bdb:"+options.path+"#convergence.results")	
		db[s+"res_fsc"] = [x,fsc] # warning, changing this naming convention will disrupt forms in the workflow (make them fail)
		db[s+"res_datapow"] = [x,datapow]
		db[s+"res_noisepow"] = [x,noisepow]
		db_close_dict("bdb:"+options.path+"#convergence.results")
		

	E2end(E2n)
	

				
if __name__ == "__main__":
    main()
