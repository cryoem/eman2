from libpyEMData2 import *
from libpyFactory2 import *
from libpyGeometry2 import *
from libpyUtils2 import *

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
