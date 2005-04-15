# Helper functions from Pawel Penczek
# Please do not alter this file without permision from the author.

from EMAN2 import *

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
		ny = image.get_xsize()
		nz = image.get_xsize()
	except:
		# hopefully the "image" is actually a filename
		e = EMData()
		e.read_image(image)
		mean = e.get_attr("mean")
		sigma = e.get_attr("sigma")
		imin = e.get_attr("minimum")
		imax = e.get_attr("maximum")
		nx = e.get_xsize()
		ny = e.get_xsize()
		nz = e.get_xsize()
	print "Image size: nx = %i, ny = %i, nz = %i" % (nx, ny, nz)
	print "avg = %g, std dev = %g, min = %g, max = %g" % (mean, sigma, imin, imax)
	return mean,sigma,imin,imax, nx, ny, nz
	
