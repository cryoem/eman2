from libpyEMData2 import *
from libpyFactory2 import *
from libpyGeometry2 import *
from libpyUtils2 import *
from libpyPointArray2 import *
from bisect import bisect_left

EMANVERSION="EMAN2 v1.90"

MRC = EMUtil.ImageType.IMAGE_MRC
SPIDER = EMUtil.ImageType.IMAGE_SPIDER
SINGLE_SPIDER = EMUtil.ImageType.IMAGE_SINGLE_SPIDER
IMAGIC = EMUtil.ImageType.IMAGE_IMAGIC
HDF = EMUtil.ImageType.IMAGE_HDF
DM3 = EMUtil.ImageType.IMAGE_DM3
TIFF = EMUtil.ImageType.IMAGE_TIFF
PGM = EMUtil.ImageType.IMAGE_PGM
LST = EMUtil.ImageType.IMAGE_LST
PIF = EMUtil.ImageType.IMAGE_PIF
VTK = EMUtil.ImageType.IMAGE_VTK
PNG = EMUtil.ImageType.IMAGE_PNG
SAL = EMUtil.ImageType.IMAGE_SAL
ICOS = EMUtil.ImageType.IMAGE_ICOS
EMIM = EMUtil.ImageType.IMAGE_EMIM
GATAN2 = EMUtil.ImageType.IMAGE_GATAN2
AMIRA = EMUtil.ImageType.IMAGE_AMIRA
XPLOR = EMUtil.ImageType.IMAGE_XPLOR
EM = EMUtil.ImageType.IMAGE_EM
IMAGE_UNKNOWN = EMUtil.ImageType.IMAGE_UNKNOWN

def parse_filter_params(filterparams):
    params = filterparams.split(":")
    filtername = params[0]

    if len(params) == 1:
        return (filtername, None)
    else:
        d = Dict()
        for param in params[1:]:
            key_values = param.split("=")
            d[key_values[0]] = EMObject(key_values[1])
        return (filtername, d)


def get_optionlist(argv):
    optionlist = []
    for arg1 in argv:
        if arg1[0] == "-":
            argname = arg1.split("=")
            optionlist.append(argname[0].lstrip("-"))
    return optionlist

# These are particle box sizes which have only 2,3,5 and 7 as prime factors, 
# and are all divisible by 4
good_box_sizes=[24, 28, 32, 36, 40, 48, 56, 60, 64, 72, 80, 84, 96, 100, 108, 112, 120, 128, 
	140, 144, 160, 168, 180, 192, 196, 200, 216, 224, 240, 252, 256, 280, 288, 300, 320, 324, 
	336, 360, 384, 392, 400, 420, 432, 448, 480, 500, 504, 512, 540, 560, 576, 588, 600, 640, 
	648, 672, 700, 720, 756, 768, 784, 800, 840, 864, 896, 900, 960, 972, 980, 1000, 1008, 1024, 
	1080, 1120, 1152, 1176, 1200, 1260, 1280, 1296, 1344, 1372, 1400, 1440, 1500, 1512, 1536, 
	1568, 1600, 1620, 1680, 1728, 1764, 1792, 1800, 1920, 1944, 1960, 2000, 2016, 2048, 2100, 
	2160, 2240, 2268, 2304, 2352, 2400, 2500, 2520, 2560, 2592, 2688, 2700, 2744, 2800, 2880, 
	2916, 2940, 3000, 3024, 3072, 3136, 3200, 3240, 3360, 3456, 3500, 3528, 3584, 3600, 3780, 
	3840, 3888, 3920, 4000, 4032, 4096]

def good_boxsize(val):
	"This will find the next largest 'good' boxsize"
	if val > 4096: return int(2**ceil(log(val)/log(2.0)))
	return good_box_sizes[bisect_left(good_box_sizes,val)]
