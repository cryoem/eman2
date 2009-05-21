/**
 * $Id$
 */

/*
 * Author: David Woolford, 11/06/2007 (woolford@bcm.edu)
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

#ifdef EMAN2_USING_OPENGL

#include "glutil.h"
#include "emdata.h"

#ifndef _WIN32
	#ifndef GL_GLEXT_PROTOTYPES
		#define GL_GLEXT_PROTOTYPES
	#endif	//GL_GLEXT_PROTOTYPES
#endif	//_WIN32

#ifdef __APPLE__
	#include "OpenGL/gl.h"
	#include "OpenGL/glu.h"
	#include "OpenGL/glext.h"
#else // WIN32, LINUX
	#include "GL/gl.h"
	#include "GL/glu.h"
	#include "GL/glext.h"
#endif	//__APPLE__

using namespace EMAN;

unsigned int GLUtil::gen_glu_mipmaps(const EMData* const emdata)
{
	if ( emdata->get_data() == 0 ) throw NullPointerException("Error, attempt to create an OpenGL mipmap without internally stored data");
	ENTERFUNC;

	unsigned int tex_name;
	glGenTextures(1, &tex_name);

	if ( emdata->ny == 1 && emdata->nz == 1 )
	{
		glBindTexture(GL_TEXTURE_1D, tex_name);
		gluBuild1DMipmaps(GL_TEXTURE_1D, GL_LUMINANCE, emdata->nx, GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
	} else if (emdata->nz == 1) {
		glBindTexture(GL_TEXTURE_2D, tex_name);
		gluBuild2DMipmaps(GL_TEXTURE_2D, GL_LUMINANCE, emdata->nx, emdata->ny, GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
	}
	else {
#ifdef	_WIN32
		//There is no gluBuild3DMipmaps() function in glu.h for VS2003 and VS2005
		printf("3D OpenGL mipmaps are not available on this platform.\n");
#else
		glBindTexture(GL_TEXTURE_3D, tex_name);
		gluBuild3DMipmaps(GL_TEXTURE_3D, GL_LUMINANCE, emdata->nx, emdata->ny, emdata->nz, GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
#endif	//_WIN32
	}

	EXITFUNC;
	return tex_name;
}

unsigned int GLUtil::gen_gl_texture(const EMData* const emdata)
{
	if ( emdata->get_data() == 0 ) throw NullPointerException("Error, attempt to create an OpenGL texture without internally stored data");
	ENTERFUNC;

	unsigned int tex_name;
	glGenTextures(1, &tex_name);

	if ( emdata->ny == 1 && emdata->nz == 1 )
	{
		glBindTexture(GL_TEXTURE_1D, tex_name);
		glTexImage1D(GL_TEXTURE_1D,0,GL_LUMINANCE,emdata->nx,0,GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
	} else if (emdata->nz == 1) {
		glBindTexture(GL_TEXTURE_2D, tex_name);
		glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE,emdata->nx,emdata->ny,0,GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
	}
	else {
		glBindTexture(GL_TEXTURE_3D, tex_name);
#ifdef _WIN32
	PFNGLTEXIMAGE3DPROC glTexImage3D;
#endif
		glTexImage3D(GL_TEXTURE_3D,0, GL_LUMINANCE,emdata->nx,emdata->ny,emdata->nz,0,GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
	}

	EXITFUNC;
	return tex_name;
}

unsigned int GLUtil::render_amp8_gl_texture(EMData* emdata, int x0, int y0, int ixsize, int iysize, int bpl, float scale, int mingray, int maxgray,	float render_min, float render_max,float gamma,int flags) {

	string pixels = emdata->render_amp8(x0, y0, ixsize,iysize, bpl, scale, mingray, maxgray, render_min, render_max, gamma,flags);

	unsigned int tex_name;
	glGenTextures(1, &tex_name);

	glBindTexture(GL_TEXTURE_2D,tex_name);
	glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE,ixsize,iysize,0,GL_LUMINANCE, GL_UNSIGNED_BYTE, pixels.c_str() );

	return tex_name;
}

// undef GL_GLEXT_PROTOTYPES
#ifdef GL_GLEXT_PROTOTYPES
#undef GL_GLEXT_PROTOTYPES
#endif

#endif // EMAN2_USING_OPENGL
