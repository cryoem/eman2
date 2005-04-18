# Helper functions from Pawel Penczek
# Please do not alter this file without permision from the author.

from EMAN2 import *

def readImage(filename):
	"""Read an image from the disk.

	Usage: myimage = readImage("path/to/image")
	"""
	image = EMData()
	image.read_image(filename)
	return image

def descriptive_statistics(image):
	"""Calculate the descriptive statistics of an image.

	Usage: [mean, sigma, xmin, xmax, nx, ny, nz =] descriptive_statistics(image object)
	       or
	       [mean, sigma, xmin, xmax, nx, ny, nz =] descriptive_statistics(image filename)

	Purpose: calculate basic statistical characteristics of an image.
	"""
	try:
		mean = image.get_attr("mean")
		sigma = image.get_attr("sigma")
		imin = image.get_attr("minimum")
		imax = image.get_attr("maximum")
		nx = image.get_xsize()
		ny = image.get_ysize()
		nz = image.get_zsize()
	except:
		# hopefully the "image" is actually a filename
		e = EMData()
		e.read_image(image)
		mean = e.get_attr("mean")
		sigma = e.get_attr("sigma")
		imin = e.get_attr("minimum")
		imax = e.get_attr("maximum")
		nx = e.get_xsize()
		ny = e.get_ysize()
		nz = e.get_zsize()
	print "Image size: nx = %i, ny = %i, nz = %i" % (nx, ny, nz)
	print "avg = %g, std dev = %g, min = %g, max = %g" % (mean, sigma, imin, imax)
	return mean,sigma,imin,imax, nx, ny, nz
	

def printImage(image):
	"""Print the data in an image to standard out.

	Usage: printImage(image)
	   or
	       printImage("path/to/image")
	"""
	try:
		nx = image.get_xsize()
	except: # homefully a filename
		image = readImage(image)
		nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	for iz in xrange(nz):
		print "(z = %d slice)" % (iz)
		line = []
		for ix in xrange(nx):
			for iy in xrange(ny):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((iy + 1) % 5 == 0):
					line.append("\n   ")
			line.append("\n")
		print "".join(line)




