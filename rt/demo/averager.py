#!/usr/bin/python
# usage: averager.py inputfile outputfile

from EMAN2 import *
import sys


images = EMData.read_images(sys.argv[1])
averager = Averagers.get("Image")
avg_image = averager.average(images)
avg_image.write_image(sys.argv[2])
