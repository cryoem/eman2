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
#ifndef glutil_h__
#define glutil_h__

#include <vector>
#include "vec3.h"
#include "transform.h"

#ifdef __APPLE__
	#include "OpenGL/gl.h"
#else // WIN32, LINUX
	#include "GL/gl.h"
#endif	//__APPLE__

using std::vector;

namespace EMAN
{
	class EMData;
	class MarchingCubes;

	class GLUtil {
	public:
		/** create an OpenGL mipmap set
		 * @return the texture id used in the call to glBindTextures
		 */
		static unsigned int gen_glu_mipmaps(const EMData* const emdata);

		/** create an OpenGL texture
		 * @param emdata: The image data that will be used to create the texture
		 * @param format: glTexImage() sets its internalFormat and format parameters to this value. Can be GL_COLOR_INDEX, GL_RED, GL_GREEN, GL_BLUE, GL_ALPHA, GL_RGB, GL_BGR, GL_RGBA, GL_BGRA, GL_LUMINANCE, or GL_LUMINANCE_ALPHA
		 * @return the texture id used in the call to glBindTextures
		 */
		static unsigned int gen_gl_texture(const EMData* const emdata, GLenum format = GL_LUMINANCE);

		/** create an OpenGL texture using render_amp8
		 * @return the texture id used in the call to glBindTextures
		 */
		static unsigned int render_amp8_gl_texture(EMData* emdata, int x0, int y0, int ixsize, int iysize, int bpl, float scale, int mingray, int maxgray,	float render_min, float render_max,float gamma,int flags);

		/** Determine the intersection of .... just ask David Woolford
		 *
		 */
		static int nearest_projected_points(const vector<float>& model_matrix, const vector<float>& proj_matrix, const vector<int>& view_matrix, const vector<Vec3f>& points, const float mouse_x, const float mouse_y,const float& nearnes);
		static void colored_rectangle(const vector<float>& data,const float& alpha, const bool center_point=false);
		static void mx_bbox(const vector<float>& data, const vector<float>& text_color, const vector<float>& bg_color);

		/** Render the image into an 8-bit image. 2D images only.
		 * flags provide a way to do unusual things with this function, such
		 * as calculating a histogram of the rendered area.
		 *
		 * @param x	origin of the area to render
		 * @param y
		 * @param xsize	size of the area to render in output pixels
		 * @param ysize
		 * @param bpl	bytes per line, if asrgb remember *3
		 * @param scale	scale factor for rendering
		 * @param min_gray	minimum gray value to render (0-255)
		 * @param max_gray	maximum gray value to render (0-255)
		 * @param min_render	float image density corresponding to min_gray
		 * @param max_render	float image density corresponding to max_gray
		 * @param gamma
		 * @param flags	1-RGB (24 bit) rendering,2-add a 256 int greyscale histogram to the end of the image array,4-invert y axis,8-render 32 bit 0xffRRGGBB,16-Color display of complex images
		 * @exception ImageDimensionException If the image is not 2D.
		 */
		static std::string render_amp8(EMData * emdata, int x, int y, int xsize, int ysize,
						 int bpl, float scale, int min_gray, int max_gray,
						 float min_render, float max_render,float gamma,int flags);

		/** Get an isosurface display list
		* Traverses the tree, marches the cubes, and renders a display list using the associated vertices and normals
		* Uses OpenGL arrays for maximum performance
		* DEPRECATED
		* @return an OpenGL display list number
		*/
		static unsigned long get_isosurface_dl(MarchingCubes* mc, unsigned int tex_id = 0, bool surface_face_z = false, bool recontour = true);
		
		/** Render a isosurface using buffer objects, this uses non-deprecated methods and improves performance */
		static void render_using_VBOs(MarchingCubes* mc, unsigned int tex_id = 0, bool surface_face_z = false);
		
		/** Recountour isosurface, for use with VBOs */
		static void contour_isosurface(MarchingCubes* mc);
		
		/** Load a EMAN style transform to open GL w/o having to go through python
		* Calls glLoadTransposeMatrix rather than glLoadMatrix to convert between C/C++/Python row-major format and openGL's Column major format
		* @param xform The Transform to apply
		**/
		static void glLoadMatrix(const Transform& xform);
		
		/** Mult a EMAN style transform to open GL w/o having to go through python
		* Calls glMultTransposeMatrix rather than glMultMatrix to convert between C/C++/Python row-major format and openGL's Column major format
		* @param xform The Transform to apply
		**/
		static void glMultMatrix(const Transform& xform);
		
		/** Draw a bounding box. This is done on the C side to simplify vertex arrays */
		static void glDrawBoundingBox(float width, float height, float depth);
		
		static void glDrawDisk(float radius, int spokes);
		
	private:
		//This is a buffer for the bounding box
		static GLuint buffer[2];
		
	};
}

#endif	//glutil_h__
