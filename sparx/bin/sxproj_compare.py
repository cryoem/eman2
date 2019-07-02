#!/usr/bin/env python
import os
import EMAN2
import EMAN2_cppwrap
#from EMAN2 import EMUtil, EMArgumentParser, EMANVERSION
from applications import header, project3d
from utilities import get_im, write_header, model_circle, read_text_row
from statistics import ccc
from fundamentals import rops_table, fft
from projection import prep_vol, prgl
from math import sqrt
from filter import filt_table

# TO DO:
#	resize the class-averages and re-projections if they have different sizes?

USAGE = """
Write to the same directory as the input images:
sxproj_compare.py <input_imgs> <input_volume>

Write to a specific directory:
sxproj_compare.py <input_imgs> <input_volume> --outdir <output_directory>

Supply a file with projection angles if not present in the header of the input images:
sxproj_compare.py <input_imgs> <input_volume> --angles <angles_file>

Some input images may have been assigned projection angles, so include an image selection file:
sxproj_compare.py <input_imgs> <input_volume> --angles <angles_file> --select <img_selection_file>

Choose interpolation method for re-projection of input volume:
sxproj_compare.py <input_imgs> <input_volume> --prjmethod <interpolation_method>
Valid are options are: trilinear (default), gridding, and nn (nearest neighbor)

Automatically open a montage of output images:
sxproj_compare.py <input_imgs> <input_volume> --display
"""
	
def runcheck(classavgstack, reconfile, outdir, inangles=None, selectdoc=None, prjmethod='trilinear', displayYN=False, 
			 projstack='proj.hdf', outangles='angles.txt', outstack='comp-proj-reproj.hdf', normstack='comp-proj-reproj-norm.hdf'):
	
	print("\n%s, Modified 2018-12-07\n" % __file__)
	
	# Check if inputs exist
	check(classavgstack)
	check(reconfile)
	
	# Create directory if it doesn't exist
	if not os.path.isdir(outdir):
		os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
		print("mkdir -p %s" % outdir)

	# Expand path for outputs
	projstack = os.path.join(outdir, projstack)
	outangles = os.path.join(outdir, outangles)
	outstack  = os.path.join(outdir, outstack)
	normstack = os.path.join(outdir, normstack)
	
	# Get number of images
	nimg0 = EMAN2_cppwrap.EMUtil.get_image_count(classavgstack)
	recon = EMAN2_cppwrap.EMData(reconfile)
	nx = recon.get_xsize()
	
	# In case class averages include discarded images, apply selection file
	if selectdoc:
		goodavgs, extension = os.path.splitext(classavgstack)
		newclasses = goodavgs + "_kept" + extension
		
		# e2proc2d appends to existing files, so rename existing output
		if os.path.exists(newclasses):
			renamefile = newclasses + '.bak'
			os.rename(newclasses, renamefile)
			print("mv %s %s" % (newclasses, renamefile))
		
		cmd7="e2proc2d.py %s %s --list=%s" % (classavgstack, newclasses, selectdoc)
		print(cmd7)
		os.system(cmd7)
		
		# Update class-averages
		classavgstack = newclasses
	
	# Import Euler angles
	if inangles:
		cmd6 = "sxheader.py %s --params=xform.projection --import=%s" % (classavgstack, inangles)
		print(cmd6)
		header(classavgstack, 'xform.projection', fimport=inangles)
	
	try:
		header(classavgstack, 'xform.projection', fexport=outangles)
		cmd1 = "sxheader.py %s --params=xform.projection --export=%s" % (classavgstack, outangles) 
		print(cmd1)
	except RuntimeError:
		print("\nERROR!! No projection angles found in class-average stack header!\n")
		print('Usage:', USAGE)
		exit()
	
	#cmd2="sxproject3d.py %s %s --angles=%s" % (recon, projstack, outangles)
	#print(cmd2)
	#os.system(cmd2)
	
	#  Here if you want to be fancy, there should be an option to chose the projection method,
	#  the mechanism can be copied from sxproject3d.py  PAP
	if prjmethod=='trilinear':
		method_num = 1
	elif prjmethod=='gridding':
		method_num = -1
	elif prjmethod=='nn':
		method_num = 0
	else:
		print("\nERROR!! Valid projection methods are: trilinear (default), gridding, and nn (nearest neighbor).")
		print('Usage:', USAGE)
		exit()
	
	#project3d(recon, stack=projstack, listagls=outangles)
	recon = prep_vol(recon, npad = 2, interpolation_method = 1)

	result=[]
	#  Here you need actual radius to compute proper ccc's, but if you do, you have to deal with translations, PAP
	mask = model_circle(nx//2-2,nx,nx)
	
	# Number of images may have changed
	nimg1   = EMAN2_cppwrap.EMUtil.get_image_count(classavgstack)
	outangles = read_text_row(outangles)
	for imgnum in range(nimg1):
		# get class average
		classimg = get_im(classavgstack, imgnum)
		
		# compute re-projection
		prjimg = prgl(recon, outangles[imgnum], 1, False)
		
		# calculate 1D power spectra
		rops_dst = rops_table(classimg*mask)  
		rops_src = rops_table(prjimg)
		
		#  Set power spectrum of reprojection to the data.
		#  Since data has an envelope, it would make more sense to set data to reconstruction,
		#  but to do it one would have to know the actual resolution of the data. 
		#  you can check sxprocess.py --adjpw to see how this is done properly  PAP
		table = [0.0]*len(rops_dst)  # initialize table
		for j in range( len(rops_dst) ):
			table[j] = sqrt( rops_dst[j]/rops_src[j] )
		prjimg = fft(filt_table(prjimg, table))  # match FFT amplitdes of re-projection and class average

		cccoeff = ccc(prjimg, classimg, mask)
		#print(imgnum, cccoeff)
		classimg.set_attr_dict({'cross-corr':cccoeff})
		prjimg.set_attr_dict({'cross-corr':cccoeff})
		prjimg.write_image(outstack,2*imgnum)
		classimg.write_image(outstack, 2*imgnum+1)
		result.append(cccoeff)
	del outangles
	meanccc = sum(result)/nimg1
	print("Average CCC is %s" % meanccc)

	nimg2 = EMAN2_cppwrap.EMUtil.get_image_count(outstack)
	
	for imgnum in xrange(nimg2):
		if (imgnum % 2 ==0):
			prjimg = get_im(outstack,imgnum)
			meanccc1 = prjimg.get_attr_default('mean-cross-corr', -1.0)
			prjimg.set_attr_dict({'mean-cross-corr':meanccc})
			write_header(outstack,prjimg,imgnum)
		if (imgnum % 100) == 0:
			print(imgnum)
	
	# e2proc2d appends to existing files, so delete existing output
	if os.path.exists(normstack):
		os.remove(normstack)
		print("rm %s" % normstack)
		


	#  Why would you want to do it?  If you do, it should have been done during ccc calculations,
	#  otherwise what is see is not corresponding to actual data, thus misleading.  PAP
	#cmd5="e2proc2d.py %s %s --process=normalize" % (outstack, normstack)
	#print(cmd5)
	#os.system(cmd5)
	
	# Optionally pop up e2display
	if displayYN:
		cmd8 = "e2display.py %s" % outstack
		print(cmd8)
		os.system(cmd8)
	
	print("Done!")
	
def check(file):
	if not os.path.exists(file):
		print("ERROR!! %s doesn't exist!\n" % file)
		exit()

	
if __name__ == "__main__":
	# Command-line arguments
	parser = EMAN2_cppwrap.EMArgumentParser(usage=USAGE,version=EMAN2.EMANVERSION)
	parser.add_argument('classavgs', help='Input class averages')
	parser.add_argument('vol3d', help='Input 3D reconstruction')
	parser.add_argument('--outdir', "-o", type=str, help='Output directory')
	parser.add_argument('--angles', "-a", type=str, help='Angles files, which will be imported into the header of the class-average stack')
	parser.add_argument('--select', "-s", type=str, help='Selection file for included classes. RVIPER may exclude some images from the reconstruction.')
	parser.add_argument('--prjmethod', "-p", type=str, default='trilinear', help='Projection method: trilinear (default), gridding, nn (nearest neighbor)')
	parser.add_argument('--display', "-d", action="store_true", help='Automatically open montage in e2display')
	
	(options, args) = parser.parse_args()
	#print args, options  # (Everything is in options.)
	#exit()
	
	# If output directory not specified, write to same directory as class averages
	if not options.outdir:
		outdir = os.path.dirname(os.path.realpath(options.classavgs))
	else:
		outdir = options.outdir
	#print("outdir: %s" % outdir)

	runcheck(options.classavgs, options.vol3d, outdir, 
		  inangles=options.angles, selectdoc=options.select, prjmethod=options.prjmethod, displayYN=options.display)
