#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from optparse import OptionParser
import pprint
from EMAN2 import *

def get_data_type_string(datatype):
    dtstring = {
        0 : "UNKNOWN",
        1 : "CHAR",
        2 : "UNSIGNED CHAR",
        3 : "SHORT",
        4 : "UNSIGNED SHORT",
        5 : "INT",
        6 : "UNSIGNED INT",
        7 : "FLOAT",
        8 : "DOUBLE",
        9 : "SHORT_COMPLEX",
        10 : "USHORT_COMPLEX",
        11 : "FLOAT_COMPLEX"
    }
    return dtstring.get(datatype, 'UNKNOWN')

def main():
    progname = os.path.basename(sys.argv[0])
    usage = """%prog [options] imagefile
This program will print out some information about the image.
    """
    
    parser = OptionParser(usage=usage,version=EMANVERSION)
    
    parser.add_option("-H", "--header", action="store_true",help="show all header information",default=False)
    parser.add_option("-v", "--verbose", type="int", help="set verbosity level. N=0,1,2,3. large N means more verbose.")
    parser.add_option("-s", "--stat", action="store_true",help="show statistical information about the image(s).",default=False)
    
    (options, args) = parser.parse_args()
    
    if len(args)<1:
        print usage
        parser.error("Specify image file")
    
    if options.header:
        show_all_header = True
    else:
        show_all_header = False
    if options.stat:
        stat = True
    else:
        stat = False
    
    imagefile = args[0]
    nimg = EMUtil.get_image_count(imagefile)
    imgtype = EMUtil.get_image_type(imagefile)
    imgtypename = EMUtil.get_imagetype_name(imgtype)
    image_index = 0
    if imgtype == EMUtil.ImageType.IMAGE_SPIDER and not stat: 
        image_index = -1
    print '%20s: %d' % ('Number of Images', nimg)
    print '%20s: %s' % ('Image Format',imgtypename)
    
    d = EMData()
    if not stat:
        d.read_image(imagefile, image_index, True)
    else:
        d.read_image(imagefile, image_index, False)
    
    print '%20s: %d x %d x %d' % ("Image Dimensions", d.get_xsize(), d.get_ysize(), d.get_zsize())
    print '%20s: %s' % ("Image Data Type", get_data_type_string(d.get_attr("datatype")))
    if options.stat: 
        print '\nmean=%1.3g sigma=%1.3g skewness=%1.3g kurtosis=%1.3g' % \
            ( d.get_attr("mean"), d.get_attr("sigma"), d.get_attr("skewness"), d.get_attr("kurtosis"))
    
    ctf = d.get_ctf()
    if ctf is not None:
        print '\nCTF: %s' % (ctf.to_string())
    
    if options.header:
        dict = d.get_attr_dict()
        print '\nDetailed Header Information:'
        pprint.pprint(dict)

if __name__ == "__main__":
    main()