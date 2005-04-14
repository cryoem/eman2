# Helper functions from Pawel Penczek
# Please do not alter this file without permision from the author.

from EMAN2 import *

def fs(image):
	try:
		mean = image.get_attr("mean")
		sigma = image.get_attr("sigma")
		imin = image.get_attr("minimum")
		imax = image.get_attr("maximum")
	except:
		# hopefully the "image" is actually a filename
		e = EMData()
		e.read_image(image)
		mean = e.get_attr("mean")
		sigma = e.get_attr("sigma")
		imin = e.get_attr("minimum")
		imax = e.get_attr("maximum")
	print "avg = %g, std dev = %g, min = %g, max = %g" % (mean, sigma, imin, imax)
	
	
