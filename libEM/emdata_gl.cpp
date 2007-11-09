
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

#include "emdata.h"
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES
#endif
#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glext.h"
using namespace EMAN;

unsigned int EMData::gen_glu_mipmaps() const
{
	if ( get_data() == 0 ) throw NullPointerException("Error, attempt to create an OpenGL mipmap without internally stored data");
	ENTERFUNC;
	
	unsigned int tex_name;
	glGenTextures(1, &tex_name);
	
	if ( ny == 1 and nz == 1 )
	{
		glBindTexture(GL_TEXTURE_1D, tex_name);
		gluBuild1DMipmaps(GL_TEXTURE_1D, GL_LUMINANCE, nx, GL_LUMINANCE, GL_FLOAT, (void*)get_data());
	} else if (nz == 1) {
		glBindTexture(GL_TEXTURE_2D, tex_name);
		gluBuild2DMipmaps(GL_TEXTURE_2D, GL_LUMINANCE, nx, ny, GL_LUMINANCE, GL_FLOAT, (void*)get_data());
	}
	else {
		glBindTexture(GL_TEXTURE_3D, tex_name);
		gluBuild3DMipmaps(GL_TEXTURE_3D, GL_LUMINANCE, nx, ny, nz, GL_LUMINANCE, GL_FLOAT, (void*)get_data());
	}
	
	return tex_name;
}

unsigned int EMData::gen_gl_texture() const
{
	if ( get_data() == 0 ) throw NullPointerException("Error, attempt to create an OpenGL mipmap without internally stored data");
	ENTERFUNC;
	
	unsigned int tex_name;
	glGenTextures(1, &tex_name);
		
	if ( ny == 1 and nz == 1 )
	{
		glBindTexture(GL_TEXTURE_1D, tex_name);
		glTexImage1D(GL_TEXTURE_1D,0,GL_LUMINANCE,nx,0,GL_LUMINANCE, GL_FLOAT, (void*)get_data());
	} else if (nz == 1) {
		glBindTexture(GL_TEXTURE_2D, tex_name);
		glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE,nx,ny,0,GL_LUMINANCE, GL_FLOAT, (void*)get_data());
	}
	else {
		glBindTexture(GL_TEXTURE_3D, tex_name);
		glTexImage3D(GL_TEXTURE_3D,0, GL_LUMINANCE,nx,ny,nz,0,GL_LUMINANCE, GL_FLOAT, (void*)get_data());
	}
	
	return tex_name;
}

#endif // EMAN2_USING_OPENGL
