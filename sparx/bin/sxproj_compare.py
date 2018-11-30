#!/usr/bin/env python
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT from EMAN2 import EMUtil, EMArgumentParser, EMANVERSION
pass#IMPORTIMPORTIMPORT from applications import header, project3d
pass#IMPORTIMPORTIMPORT from utilities import get_im, write_header
pass#IMPORTIMPORTIMPORT from statistics import ccc

# TO DO:
#	resize the class-averages and re-projections if they have different sizes?

import EMAN2
import EMAN2_cppwrap
import EMAN2_meta
import applications
import os
import statistics
import utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import EMAN2_meta
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import utilities
def runcheck(classavgstack, recon, outdir, inangles=None, selectdoc=None, displayYN=False, 
			 projstack='proj.hdf', outangles='angles.txt', outstack='comp-proj-reproj.hdf', normstack='comp-proj-reproj-norm.hdf'):
	
	print
	
	# Check if inputs exist
	check(classavgstack)
	check(recon)
	
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
	#print("nimg0: %s" % nimg0)
	
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
		print cmd7
		os.system(cmd7)
		
		# Update class-averages
		classavgstack = newclasses
	
	# Import Euler angles
	if inangles:
		cmd6="sxheader.py %s --params=xform.projection --import=%s" % (classavgstack, inangles)
		print cmd6
		applications.header(classavgstack, 'xform.projection', fimport=inangles)
	
	cmd1="sxheader.py %s --params=xform.projection --export=%s" % (classavgstack, outangles) 
	print cmd1
	#os.system(cmd1)
	try:
		applications.header(classavgstack, 'xform.projection', fexport=outangles)
	except RuntimeError:
		print("\nERROR!! No projection angles found in class-average stack header!\n")
		exit()
	
	cmd2="sxproject3d.py %s %s --angles=%s" % (recon, projstack, outangles)
	print cmd2
	#os.system(cmd2)
	applications.project3d(recon, stack=projstack, listagls=outangles)
	
	imgcounter = 0  # montage will have double the number of images as number of class-averages
	result=[]
	
	# Number of images may have changed
	nimg1   = EMAN2_cppwrap.EMUtil.get_image_count(classavgstack)
	
	for imgnum in xrange(nimg1):
		#print imgnum
		classimg = utilities.get_im(classavgstack, imgnum)
		ccc1 = classimg.get_attr_default('cross-corr', -1.0)
		prjimg = utilities.get_im(projstack,imgnum)
		ccc1 = prjimg.get_attr_default('cross-corr', -1.0)
		cccoeff = statistics.ccc(prjimg,classimg)
		#print imgnum, cccoeff
		classimg.set_attr_dict({'cross-corr':cccoeff})
		prjimg.set_attr_dict({'cross-corr':cccoeff})
		prjimg.write_image(outstack,imgcounter)
		imgcounter += 1
		classimg.write_image(outstack, imgcounter)
		imgcounter += 1
		result.append(cccoeff)
	result1 = sum(result)
	#print result1

	nimg2   = EMAN2_cppwrap.EMUtil.get_image_count(outstack)
	meanccc = result1/nimg1
	print("Mean CCC is %s" % meanccc)
	
	for imgnum in xrange(nimg2):
		if (imgnum % 2 ==0):
			prjimg = utilities.get_im(outstack,imgnum)
			meanccc1 = prjimg.get_attr_default('mean-cross-corr', -1.0)
			prjimg.set_attr_dict({'mean-cross-corr':meanccc})
			utilities.write_header(outstack,prjimg,imgnum)
		if (imgnum % 100) == 0:
			print imgnum
	
	# e2proc2d appends to existing files, so delete existing output
	if os.path.exists(normstack):
		os.remove(normstack)
		print("rm %s" % normstack)
		
	cmd5="e2proc2d.py %s %s --process=normalize" % (outstack, normstack)
	print cmd5
	os.system(cmd5)
	
	# Optionally pop up e2display
	if displayYN:
		cmd8 = "e2display.py %s" % normstack
		print cmd8
		os.system(cmd8)
	
	print("Done!")
	
def check(file):
	if not os.path.exists(file):
		print("ERROR!! %s doesn't exist!\n" % file)
		exit()

	
if __name__ == "__main__":
	usage = """
	sxproj_compare.py <input_imgs> <input_volume>
	sxproj_compare.py <input_imgs> <input_volume> --outdir <output_directory>
	sxproj_compare.py <input_imgs> <input_volume> --outdir <output_directory> --angles <angles_file>
	sxproj_compare.py <input_imgs> <input_volume> --outdir <output_directory> --angles <angles_file> --select <img_selection_file>
	sxproj_compare.py <input_imgs> <input_volume> --outdir <output_directory> --angles <angles_file> --select <img_selection_file> --display
	"""
	
	# Command-line arguments
	parser = EMAN2.EMArgumentParser(usage=usage,version=EMAN2_meta.EMANVERSION)
	parser.add_argument('classavgs', help='Input class averages')
	parser.add_argument('vol3d', help='Input 3D reconstruction')
	parser.add_argument('--outdir', "-o", type=str, help='Output directory')
	parser.add_argument('--angles', "-a", type=str, help='Angles files, which will be imported into the header of the class-average stack')
	parser.add_argument('--select', "-s", type=str, help='Selection file for included classes. RVIPER may exclude some images from the reconstruction.')
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

	runcheck(options.classavgs, options.vol3d, outdir, inangles=options.angles, selectdoc=options.select, displayYN=options.display)
