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

from libpyUtils2 import *
from libpyTransform2 import *

IMAGE_MRC = EMUtil.ImageType.IMAGE_MRC
IMAGE_SPIDER = EMUtil.ImageType.IMAGE_SPIDER
IMAGE_SINGLE_SPIDER = EMUtil.ImageType.IMAGE_SINGLE_SPIDER
IMAGE_IMAGIC = EMUtil.ImageType.IMAGE_IMAGIC
IMAGE_HDF = EMUtil.ImageType.IMAGE_HDF
IMAGE_DM3 = EMUtil.ImageType.IMAGE_DM3
IMAGE_DM4 = EMUtil.ImageType.IMAGE_DM4
IMAGE_TIFF = EMUtil.ImageType.IMAGE_TIFF
IMAGE_PGM = EMUtil.ImageType.IMAGE_PGM
IMAGE_LST = EMUtil.ImageType.IMAGE_LST
IMAGE_PIF = EMUtil.ImageType.IMAGE_PIF
IMAGE_VTK = EMUtil.ImageType.IMAGE_VTK
IMAGE_PNG = EMUtil.ImageType.IMAGE_PNG
IMAGE_SAL = EMUtil.ImageType.IMAGE_SAL
IMAGE_ICOS = EMUtil.ImageType.IMAGE_ICOS
IMAGE_EMIM = EMUtil.ImageType.IMAGE_EMIM
IMAGE_GATAN2 = EMUtil.ImageType.IMAGE_GATAN2
IMAGE_AMIRA = EMUtil.ImageType.IMAGE_AMIRA
IMAGE_XPLOR = EMUtil.ImageType.IMAGE_XPLOR
IMAGE_EM = EMUtil.ImageType.IMAGE_EM
IMAGE_V4L = EMUtil.ImageType.IMAGE_V4L
IMAGE_UNKNOWN = EMUtil.ImageType.IMAGE_UNKNOWN


#EULER_EMAN = Transform3D.EulerType.EMAN
#EULER_IMAGIC = Transform3D.EulerType.IMAGIC
#EULER_SPIN = Transform3D.EulerType.SPIN
#EULER_QUATERNION = Transform3D.EulerType.QUATERNION
#EULER_SGIROT = Transform3D.EulerType.SGIROT
#EULER_SPIDER = Transform3D.EulerType.SPIDER
#EULER_MRC = Transform3D.EulerType.MRC
#EULER_XYZ = Transform3D.EulerType.XYZ
#EULER = { "mrc":Transform3D.EulerType.MRC, "eman":Transform3D.EulerType.EMAN,
# "imagic":Transform3D.EulerType.IMAGIC,"spin":Transform3D.EulerType.SPIN,
# "quaternion":Transform3D.EulerType.QUATERNION,"sgirot":Transform3D.EulerType.SGIROT,
# "spider":Transform3D.EulerType.SPIDER,"xyz":Transform3D.EulerType.XYZ }

EM_UNKNOWN = EMUtil.EMDataType.EM_UNKNOWN
EM_CHAR = EMUtil.EMDataType.EM_CHAR
EM_UCHAR = EMUtil.EMDataType.EM_UCHAR
EM_SHORT = EMUtil.EMDataType.EM_SHORT
EM_USHORT = EMUtil.EMDataType.EM_USHORT
EM_INT = EMUtil.EMDataType.EM_INT
EM_UINT = EMUtil.EMDataType.EM_UINT
EM_FLOAT = EMUtil.EMDataType.EM_FLOAT
EM_DOUBLE = EMUtil.EMDataType.EM_DOUBLE
EM_SHORT_COMPLEX = EMUtil.EMDataType.EM_SHORT_COMPLEX
EM_USHORT_COMPLEX = EMUtil.EMDataType.EM_USHORT_COMPLEX
EM_FLOAT_COMPLEX = EMUtil.EMDataType.EM_FLOAT_COMPLEX
