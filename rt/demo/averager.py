#!/usr/bin/python
# usage: averager.py inputfile outputfile

from EMAN2 import *
import sys

images = EMData.read_images_by_index(sys.argv[1])
interation_averager = AveragerFactory.instance().get("Image")
avg_image = interation_averager.average(images)
avg_image.write_image(sys.argv[2], 0, IMAGIC)
