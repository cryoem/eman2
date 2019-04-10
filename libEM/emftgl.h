/*
 * Author: David Woolford, 7/16/2008 (woolford@bcm.edu)
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
#ifndef emftgl_h__
#define emftgl_h__

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <FTGL/ftgl.h>
#include <FTGL/FTFont.h>
namespace EMAN {


/** EMFTGL is an interface for rendering fonts in EMAN2 using FTGL

 * The EMFTGL has an instance of an EMFTGLFontManager which caches FTFonts.
 * Internally, everytime the EMFTGL is asked to render or obtain bounding boxes,
 * it asks its EMFTGLFontManager for an FTFont pointer usng the the current state
 * of all member variables. The EMFTGLFontManager may already have the correct
 * FTFont, or it may have to construct it (and store it for later use, if necessary).
 *
 * The EMFTGL class is defined in terms of 5 things, them being
 * the font size, whether or not display lists are being used (in FTGL),
 * the font file itself (a .ttf file), the mode of font rendering (TEXTURE,
 * BITMAP etc), and the depth (which is only applicable when the font mode
 * is EXTRUDE in terms of FTGL). These five parameter act as an index when asking
 * the EMFTGLFontManager for the FTFont
 *
 * The EMFTGLFontManager is intentionally not static - this is because in EMAN2
 * it is possible to be rendering accross multiple OpenGL contexts. When the EMFTGL
 * destructor is called the associated EMFTGLFontsManager is destroyed all with all
 * of its previously stored FTFont pointers.
 *
 * @author David Woolford
 * @date July 2008
*/
class EMFTGL
{
	public:
		/** Constructor is OS dependent - it attempts to reference an Operating System specific fonts file
		 * by default
		 */
#ifdef __APPLE__
		EMFTGL() : font_file_name("/Library/Fonts/Arial.ttf" ), face_size(32), depth(32), use_display_lists(true), mode(TEXTURE) {};
#else
#ifdef WIN32
		EMFTGL() : font_file_name("" ), face_size(32), depth(32), use_display_lists(true), mode(TEXTURE) {};
#else
		EMFTGL() : font_file_name("/usr/share/fonts/dejavu/DejaVuSerif.ttf"), face_size(32), depth(32), use_display_lists(true), mode(TEXTURE) {};
#endif


#endif
		/** Destructor
		 */
		~EMFTGL() {};


		/** FontModes correspond to the different FTFont types - this correspondence is
		 * EXTRUDE - FTGLExtrdFont, PIXMAP - FTGLPixmapFont
		 * TEXTURE - FTGLTextureFont, BITMAP - FTGLBitmapFont
		 * OUTLINE - FTGLOutlineFont, POLYGON - FTGLPolygonFont
		 */
		enum FontMode {
			EXTRUDE,
			PIXMAP,
			TEXTURE,
			BITMAP,
			OUTLINE,
			POLYGON
		};

		/** Render string performs OpenGL rendering of the string using the
		 * current state variables of this class
		 * @param message the string to render
		 */
		void render_string(const string& message);

		/** Obtains the bounding box of the given message, as would be rendered
		 * using the current state variables of this class
		 * @param message the string that will be used to calculate the bounding box
		 * @return a vector of length 6 containing the values [x1,y1,z1,x2,y2,z2] which define the 2 points of the bounding box extremes
		 */
		vector<float> bounding_box(const string& message);

		/** Set the font file name - should be a .ttf file
		 * @param file_name the font file name - should be a .ttf file
		 */
		void set_font_file_name(const string& file_name) { font_file_name = file_name; }

		/** Set the face size of the rendered text
		 * @param size the face size of the rendered text
		 */
		void set_face_size(const unsigned int size) { face_size = size; }

		/** Set the depth of the rendered text - only useful if this->mode = EXTRUDE
		 */
		void set_depth(const unsigned int d ) { depth = d; }

		/** Set whether or not the font render should be using display lists
		 */
		void set_using_display_lists(const bool b) { use_display_lists = b; }

		/** Set the font mode
		 */
		void set_font_mode(const FontMode m ) { mode = m; }
		/**Get the name of the current font file in use
		 */
		string get_font_file_name() { return font_file_name; }
		/** Get the currently used face size
		 */
		unsigned int get_face_size() {return face_size; }
		/** Get the currently used depth
		 */
		unsigned int get_depth() { return depth; }
		/** Get whether or not font renderer is currently using display lists
		 */
		bool get_using_display_lists() { return use_display_lists; }
		/** Get the current font mode
		 */
		FontMode get_font_mode() { return mode; }



	private:
		/** Disallow copy construction */
// 		EMFTGL(const EMFTGL& );  FIXME solve this issue
		/** Disallow Assignment */
// 		EMFTGL& operator=(const EMFTGL& );

		string font_file_name;
		unsigned int face_size;
		unsigned int depth;
		bool use_display_lists;

		FontMode mode;

		/** A class for encapsulatiing a particular instance of an FTFont (pointer)
		 * Each FTFont is characterised by 5 parameters, them being the font mode,
		 * the font file name, the face size, whether or not display lists are being
		 * used, and depth (which is redundant, except when the font mode is EXTRUDE).
		 * @author David Woolford
		 * @date July 2008
		 */
		class EMFTGLFontInstance
		{
			public:
				/** Constructor - must supply the 5 important parameters
				 */
				EMFTGLFontInstance(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int d, const bool use_dl);
				~EMFTGLFontInstance();

				/** Checks to see if the argument params match the internally stored equivalents
				 */
				bool params_match(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int depth, const bool use_dl);

				/** Get the pointer to the font
				 */
				FTFont* get_font() { return font; }
			private:
				/** Disallow copy construction */
// 				EMFTGLFontInstance(const EMFTGLFontInstance& );
				/** Disallow Assignment */
// 				EMFTGLFontInstance& operator=(const EMFTGLFontInstance& );

				EMFTGL::FontMode font_mode;
				string font_file_name;
				unsigned int face_size;
				unsigned int depth;
				bool use_display_lists;
				FTFont* font;
		};

		/** A class for managing multiple instances of EMFTGLFontInstances, in particular for
		 * caching them, for constructing news ones if they don't exist, and for returning
		 * appropriate instances.
		 * @author David Woolford
		 * @date July 2008
		 */
		class EMFTGLManager
		{
			public:
				EMFTGLManager();
// 		static EMFTGLManager& instance();

				~EMFTGLManager();

				/** Get a font with the associated parameters
				 */
				FTFont* get_font(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int d, const bool use_dl);
			private:
				/** Disallow copy construction */
// 				EMFTGLManager(const EMFTGLManager& );
				/** Disallow Assignment */
// 				EMFTGLManager& operator=(const EMFTGLManager& );

				vector<EMFTGLFontInstance*> font_instances;
		};

		EMFTGLManager fm;
};




} // namespace EMAN

#endif //emftgl_h__
