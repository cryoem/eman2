from libpyAligner2 import *
from libpyAverager2 import *
from libpyCmp2 import *
from libpyFilter2 import *
from libpyReconstructor2 import * 
from libpyProjector2 import *
from libpyEMObject2 import * 
from libpyEMData2 import *
from libpyGeometry2 import *
from libpyTransform2 import *
from libpyUtils2 import * 
from libpyPointArray2 import *
from libpyTypeConverter2 import *
from bisect import bisect_left
from pyemtbx.imagetypes import *
from pyemtbx.box import *


EMANVERSION="EMAN2 v1.90"

Vec3f.__str__=lambda x:"Vec3f"+str(x.as_list())
Transform.__str__=lambda x:"Transform(\t%1.4g\t%1.4g\t%1.4g\n\t\t%1.4g\t%1.4g\t%1.4g\n\t\t%1.4g\t%1.4g\t%1.4g)\nPretrans:%s\nPosttrans:%s"%(x.at(0,0),x.at(0,1),x.at(0,2),x.at(1,0),x.at(1,1),x.at(1,2),x.at(2,0),x.at(2,1),x.at(2,2),str(x.get_pretrans()),str(x.get_posttrans()))

def display(img):
	from os import system
	img.write_image("/tmp/img.mrc")
	system("v2 /tmp/img.mrc")
