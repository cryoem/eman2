/**
 * $Id$
 */
 
/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */
 
#ifndef eman__ioioio_h__
#define eman__ioioio_h__ 1

// NOTE: for some reason unknown, hdfio.h must be behind pngio.h. 

#include "mrcio.h"

#ifdef EM_PNG
	#include "pngio.h"
#endif	//EM_PNG

#ifdef EM_JPEG
	#include "jpegio.h"
#endif

#include "dm3io.h"

#ifdef EM_TIFF
	#include "tifio.h"
#endif	//EM_TIFF

#include "pifio.h"
#include "vtkio.h"
#include "sspiderio.h"
#include "pgmio.h"
#include "icosio.h"
#include "lstio.h"
#include "lstfastio.h"

#ifdef EM_HDF5
	#include "hdfio.h"
	#include "hdfio2.h"
#endif	//EM_HDF5

#include "salio.h"
#include "amiraio.h"
#include "xplorio.h"
#include "gatan2io.h"
#include "emio.h"
#include "imagicio.h"
#include "imagicio2.h"
#include "fitsio.h"
#include "df3io.h"
#include "omapio.h"
#include "situsio.h"

#ifdef ENABLE_V4L2
	#include "v4l2io.h"
#endif

#endif	//eman__ioioio_h__
