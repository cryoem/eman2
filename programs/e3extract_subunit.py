#!/usr/bin/env python
#
# Author: Steven Ludtke, 02/23/25 (sludtke@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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

from builtins import range
import os
import sys
from EMAN3 import *
from EMAN3jax import *



def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <particle .lst file, 2-D> <gaussian .txt> <mask>

The goal of this program is to extract one region, defined by a 3-D mask, from each 2-D particle so additionl refinement can be performed on just
the isolated region. This idea dates back to EMAN2 in the mid-2000s, and was later adopted by Relion, etc. It is not a perfect solution, as, if
the reason for doing it is conformational variability, and clearly if the structure is variable it cannot be subtracted perfectly. Nonetheless, it
can lessen the influence of the subtracted region on overall alignment/refinement and produce improved results.

To use this, you must already have run a 3-D refinement on the full particle with Gaussian reconstruction. Inputs are the .lst file from the refinement run
containing particles and orientations, the corresponding Gaussian output model and a 3-D mask used to define the region you wish to remain in the output
particles. Everything outside this mask will be subtracted. The output will be a new .hdf file with the particle data for the full stack and a corresponding
.lst file with the correct orientation.

The program proceeds in several steps:
- split the gaussian model into two models, one inside the mask and one outside the mask
- make projections of the inside and outside (and combined) models in the correct orientation for each particle
- for each particle:
-  mask the particle with a mask computed from the compbined projection
-  compute a filter to optimally match the combined projection to the particle
-  apply this filter to the outside projection, and subtract from the masked particle
-  apply the inside mask to the subtracted particle
-  center the particle based on the center of mass of the inside projection
-  optionally clip the particle down to a smaller box size
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output",type=str,default=None,help="Specify output filename (.hdf), which will contain nptcl*sym particles")
	parser.add_argument("--sym", dest = "sym", default="c1",help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos. If specified each input 'particle' will become N extracted subparticles.")
	parser.add_argument("--newbox", default=-1, type=int, help="If set, will reduce the box size of the output particles around the center of each particle")
	parser.add_argument("--threads", default=4, type=int, help="Number of threads to run in parallel on a single computer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logid=E3init(sys.argv, options.ppid)

	if options.output is None : output=args[0].replace(".lst","_sub.hdf")
	else: output=options.output
	if ":" in output:
		outbits=int(output.split(":")[-1])
		output=output.split(":")[0]
	else: outbits=6
	outlst=output.replace(".lst",".hdf")

	# The Gaussian model from a file
	gauss=Gaussians(args[1])
	if options.sym!="c1": gauss.replicate_sym(options.sym)
	maskincl=to_jax(EMData(args[2]))
	maskexcl=1.0-maskincl
	gaussincl=gauss.mask(maskincl)		# Gaussians for portion we want to keep
	gaussexcl=gauss.mask(maskexcl)		# Gaussians for portion we want to eliminate
	print (f"{len(gauss)} gaussians -> {len(gaussincl)},{len(gaussexcl)}")

	# Open the particle .lst file and get some basic info
	lsxin=LSXFile(args[0])
	lsxout=LSXFile(output.replace(".hdf",".lst"))
	hdr=lsxin.read_image(0,True)
	nx=hdr["nx"]
	N=len(lsxin)
#	N=100

	lpfilt=jax_gaussfilt_2d(nx,0.125)	# modest lowpass filter to smooth out Gaussians for mask generation

	nblk=int(1.0e9/(nx*nx*4))	# target 1G of ram at a time for the particles

	sym=Symmetries.get(options.sym)
	tlast=0
	for sn in range(sym.get_nsym()):
		sxf=sym.get_sym(sn)
		for i in range(0,N,nblk):
			tlast=print_progress(tlast,f"{i}({sn})",sn*N+i,sym.get_nsym()*N)
			ptcl=EMStack2D(EMData.read_images(args[0],range(i,min(i+nblk,N))))
			orts=ptcl.orientations_withxf(sxf)
			ortsxf=orts[0].transforms(orts[1])
#			proj=gauss.project_simple(orts[0],nx,orts[1]/nx)
			projincl=gaussincl.project_simple(orts[0],nx,orts[1]/nx)
			projexcl=gaussexcl.project_simple(orts[0],nx,orts[1]/nx)
			proj=EMStack2D(projincl.jax+projexcl.jax)

			projfullmask=(jax_ift2d(proj.do_fft().jax*lpfilt)>0.1).astype(float)
			projinclmask=(jax_ift2d(jax_fft2d(projincl.jax)*lpfilt)>0.1).astype(float)

			ptcl.set_data(ptcl.jax*projfullmask)	# mask out periphery using projection of full mask to make subtraction more accurate

			out=[]
			for im,pr,prexcl in zip(ptcl.emdata,proj.emdata,projexcl.emdata):
				#out.append(im)
				out.append(im.process("math.sub.optimal",{"ref":pr,"actual":prexcl,"return_fft":0}))
				# im.write_image("dbg_im.hdf:6",-1)
				# pr.write_image("dbg_pr.hdf:6",-1)
				# prm.write_image("dbg_prm.hdf:6",-1)

			out=EMStack2D(out)

			# Mathematically, this is kind of a stupid and inefficient way to perform recentering
			out.set_data(out.jax*projinclmask)			# apply projection of individual subunit mask after subtraction
			for im,pr in zip(out.emdata,projincl.emdata):
				pr.process_inplace("xform.centerofmass",{"int_shift_only":1})	# we compute the center of mass from the gaussian projection which shouldn't be noisy
				x,y=pr["xform.align2d"].get_trans_2d()
				im.translate(x,y,0)							# then apply it to the image

			if options.newbox>0: out=out.center_clip(options.newbox)

			out.write_images(output,outbits,i+sn*N)

			# new lst entries with correct metadata
			for j in range(i,min(i+nblk,N)):
				nxt,xt,dct=lsxin[j]
				newort=ortsxf[j-i]*sxf.inverse()
				newort.set_trans(0,0,0)		# because we center above
				dct["xform.projection"]

#				dct["xform.projection"]=sxf*dct["xform.projection"]
				lsxout[j+sn*N]=j,output,dct

	E3end(logid)






if __name__ == "__main__":
	main()
