#!/usr/bin/env python
#
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Please do not copy or modify this file without written consent of the author.
# Copyright (c) 2000-2019 The University of Texas - Houston Medical School
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
from builtins import range
import EMAN2_cppwrap
import filter
import fundamentals
import global_def
import morphology
import mpi
import optparse
import numpy
import os
import statistics
import sys
import utilities

#Transforms the local resolution file from frequency units to angstroms.
def makeAngRes(freqvol, nx, ny, nz, pxSize, freq_to_real=True):
	if (pxSize == 1.0):
		print("Using a value of 1 for the pixel size. Are you sure this is correct?")

	outAngResVol = utilities.model_blank(nx,ny,nz)
	data_in = freqvol.get_3dview()
	data_out = outAngResVol.get_3dview()

	if freq_to_real:
		mask = data_in > 0.0
		data_out[mask] = pxSize / data_in[mask]
	else:
		mask = data_in >= 2 * pxSize
		data_out[mask] = pxSize / data_in[mask]

	return outAngResVol


def output_volume(freqvol, resolut, apix, outvol, fsc, out_ang_res, nx, ny, nz, res_overall):
	outvol_ang = os.path.splitext(outvol)[0] + "_ang.hdf"
	outvol_shifted = os.path.splitext(outvol)[0] + "_shift.hdf"
	outvol_shifted_ang = os.path.splitext(outvol_shifted)[0] + "_ang.hdf"

	freqvol.write_image(outvol)
	if(out_ang_res):
		outAngResVol = makeAngRes(freqvol, nx, ny, nz, apix)
		outAngResVol.write_image(outvol_ang)

	if res_overall !=-1.0:
		for ifreq in range(len(resolut)):
			if resolut[ifreq][0] > res_overall:
				 break
		for jfreq in range(ifreq, len(resolut)):
			resolut[jfreq][1] = 0.0	

		data_freqvol = freqvol.get_3dview()
		mask = data_freqvol != 0
		percentile_25 = numpy.percentile(data_freqvol[mask], 25)
		percentile_75 = numpy.percentile(data_freqvol[mask], 75)
		iqr = percentile_75 - percentile_25
		mask_low_pass = data_freqvol > percentile_75 + 1.5*iqr
		mask_high_pass = data_freqvol < percentile_25 - 1.5*iqr
		mean_real = 1 / float(numpy.mean(data_freqvol[mask & mask_low_pass & mask_high_pass]))
		overall_res_real = 1 / float(res_overall)
		#mean_ang = options.apix / float(EMAN2_cppwrap.Util.infomask(freqvol, m, True)[0])

		volume_out_real = makeAngRes(freqvol, nx, ny, nz, 1)
		volume_out_real += (overall_res_real - mean_real)
		volume_out = makeAngRes(volume_out_real, nx, ny, nz, 1, False)
		volume_out.write_image(outvol_shifted)
		if out_ang_res:
			outAngResVol = makeAngRes(freqvol, nx, ny, nz, apix)
			outAngResVol.write_image(outvol_shifted_ang)

	if(fsc != None): utilities.write_text_row(resolut, fsc)

def main():
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + """ firstvolume  secondvolume  maskfile  outputfile  --wn  --step  --cutoff  --radius  --fsc  --res_overall  --out_ang_res  --apix  --MPI

	Compute local resolution in real space within area outlined by the maskfile and within regions wn x wn x wn
	"""
	parser = optparse.OptionParser(usage,version=global_def.SPARXVERSION)
	
	parser.add_option("--wn",           type="int",           default=7,      help="Size of window within which local real-space FSC is computed. (default 7)")
	parser.add_option("--step",         type="float",         default= 1.0,   help="Shell step in Fourier size in pixels. (default 1.0)")   
	parser.add_option("--cutoff",       type="float",         default= 0.5,   help="Resolution cut-off for FSC. (default 0.5)")
	parser.add_option("--radius",       type="int",           default=-1,     help="If there is no maskfile, sphere with r=radius will be used. By default, the radius is nx/2-wn (default -1)")
	parser.add_option("--fsc",          type="string",        default= None,  help="Save overall FSC curve (might be truncated). By default, the program does not save the FSC curve. (default none)")
	parser.add_option("--res_overall",  type="float",         default= -1.0,  help="Overall resolution at the cutoff level estimated by the user [abs units]. (default None)")
	parser.add_option("--out_ang_res",  action="store_true",  default=False,  help="Additionally creates a local resolution file in Angstroms. (default False)")
	parser.add_option("--apix",         type="float",         default= 1.0,   help="Pixel size in Angstrom. Effective only with --out_ang_res options. (default 1.0)")
	parser.add_option("--MPI",          action="store_true",  default=False,  help="Use MPI version.")

	(options, args) = parser.parse_args(arglist[1:])

	if len(args) <3 or len(args) > 4:
		print("See usage " + usage)
		sys.exit()

	if global_def.CACHE_DISABLE:
		utilities.disable_bdb_cache()

	res_overall = options.res_overall

	if options.MPI:
		sys.argv = mpi.mpi_init(len(sys.argv),sys.argv)		

		number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
		myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
		main_node = 0
		global_def.MPI = True
		cutoff = options.cutoff

		nk = int(options.wn)

		if(myid == main_node):
			#print sys.argv
			vi = utilities.get_im(sys.argv[1])
			ui = utilities.get_im(sys.argv[2])
			
			nx = vi.get_xsize()
			ny = vi.get_ysize()
			nz = vi.get_zsize()
			dis = [nx, ny, nz]
		else:
			dis = [0,0,0,0]

		global_def.BATCH = True

		dis = utilities.bcast_list_to_all(dis, myid, source_node = main_node)

		if(myid != main_node):
			nx = int(dis[0])
			ny = int(dis[1])
			nz = int(dis[2])

			vi = utilities.model_blank(nx,ny,nz)
			ui = utilities.model_blank(nx,ny,nz)


		if len(args) == 3:
			m = utilities.model_circle((min(nx,ny,nz)-nk)//2,nx,ny,nz)
			outvol = args[2]
		
		elif len(args) == 4:
			if(myid == main_node):
				m = morphology.binarize(utilities.get_im(args[2]), 0.5)
			else:
				m = utilities.model_blank(nx, ny, nz)
			outvol = args[3]
		utilities.bcast_EMData_to_all(m, myid, main_node)

		"""Multiline Comment0"""
		freqvol, resolut = statistics.locres(vi, ui, m, nk, cutoff, options.step, myid, main_node, number_of_proc)

		if(myid == 0):
			# Remove outliers based on the Interquartile range
			output_volume(freqvol, resolut, options.apix, outvol, options.fsc, options.out_ang_res, nx, ny, nz, res_overall)
		mpi.mpi_finalize()

	else:
		cutoff = options.cutoff
		vi = utilities.get_im(args[0])
		ui = utilities.get_im(args[1])

		nn = vi.get_xsize()
		nk = int(options.wn)
	
		if len(args) == 3:
			m = utilities.model_circle((nn-nk)//2,nn,nn,nn)
			outvol = args[2]
		
		elif len(args) == 4:
			m = morphology.binarize(utilities.get_im(args[2]), 0.5)
			outvol = args[3]

		mc = utilities.model_blank(nn,nn,nn,1.0)-m

		vf = fundamentals.fft(vi)
		uf = fundamentals.fft(ui)
		"""Multiline Comment1"""
		lp = int(nn/2/options.step+0.5)
		step = 0.5/lp

		freqvol = utilities.model_blank(nn,nn,nn)
		resolut = []
		for i in range(1,lp):
			fl = step*i
			fh = fl+step
			#print(lp,i,step,fl,fh)
			v = fundamentals.fft(filter.filt_tophatb( vf, fl, fh))
			u = fundamentals.fft(filter.filt_tophatb( uf, fl, fh))
			tmp1 = EMAN2_cppwrap.Util.muln_img(v,v)
			tmp2 = EMAN2_cppwrap.Util.muln_img(u,u)

			do = EMAN2_cppwrap.Util.infomask(morphology.square_root(morphology.threshold(EMAN2_cppwrap.Util.muln_img(tmp1,tmp2))),m,True)[0]


			tmp3 = EMAN2_cppwrap.Util.muln_img(u,v)
			dp = EMAN2_cppwrap.Util.infomask(tmp3,m,True)[0]
			resolut.append([i,(fl+fh)/2.0, dp/do])

			tmp1 = EMAN2_cppwrap.Util.box_convolution(tmp1, nk)
			tmp2 = EMAN2_cppwrap.Util.box_convolution(tmp2, nk)
			tmp3 = EMAN2_cppwrap.Util.box_convolution(tmp3, nk)

			EMAN2_cppwrap.Util.mul_img(tmp1,tmp2)

			tmp1 = morphology.square_root(morphology.threshold(tmp1))

			EMAN2_cppwrap.Util.mul_img(tmp1,m)
			EMAN2_cppwrap.Util.add_img(tmp1,mc)

			EMAN2_cppwrap.Util.mul_img(tmp3,m)
			EMAN2_cppwrap.Util.add_img(tmp3,mc)

			EMAN2_cppwrap.Util.div_img(tmp3,tmp1)

			EMAN2_cppwrap.Util.mul_img(tmp3,m)
			freq=(fl+fh)/2.0
			bailout = True
			for x in range(nn):
				for y in range(nn):
					for z in range(nn):
						if(m.get_value_at(x,y,z) > 0.5):
							if(freqvol.get_value_at(x,y,z) == 0.0):
								if(tmp3.get_value_at(x,y,z) < cutoff):
									freqvol.set_value_at(x,y,z, freq)
									bailout = False
								else:
									bailout = False
			if(bailout):  break
		#print(len(resolut))
		# remove outliers
		output_volume(freqvol, resolut, options.apix, outvol, options.fsc, options.out_ang_res, nx, ny, nz, res_overall)

if __name__ == "__main__":
	main()
