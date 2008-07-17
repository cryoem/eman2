
#ifndef emftgl_h__
#define emftgl_h__

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <FTGL/FTGL.h>
#include <FTGL/FTGLExtrdFont.h>
#include <FTGL/FTGLPixmapFont.h>
#include <FTGL/FTGLTextureFont.h>
namespace EMAN {

class EMFTGL
{
	public:
#ifdef LINUX
		EMFTGL() : font_file_name("/usr/share/fonts/liberation/LiberationSans-Regular.ttf" ), face_size(32), depth(32), use_display_lists(true), mode(TEXTURE) {};
#endif
#ifdef __APPLE__
		EMFTGL() : font_file_name("/Library/Fonts/Arial.ttf" ), face_size(32), depth(32), use_display_lists(true), mode(TEXTURE) {};
#endif 
#ifdef WIN32
		EMFTGL() : font_file_name("" ), face_size(32), depth(32), use_display_lists(true), mode(TEXTURE) {};
#endif
		~EMFTGL() {};
		
		void render_string(const string& message);
		vector<float> bounding_box(const string& message);
				
		void set_font_file_name(const string& file_name) { font_file_name = file_name; }
		void set_face_size(const unsigned int size) { face_size = size; }
		void set_depth(const unsigned int d ) { depth = d; }
		void set_using_display_lists(const bool b) { use_display_lists = b; }
		
		string get_font_file_name() { return font_file_name; }
		unsigned int get_face_size() {return face_size; }
		unsigned int get_depth() { return depth; }
		bool get_using_display_lists() { return use_display_lists; }
		
		enum FontMode {
			EXTRUDE,
			PIXMAP,
			TEXTURE
		};
		
	private:
		
		string font_file_name;
		unsigned int face_size;
		unsigned int depth;
		bool use_display_lists;
		
		FontMode mode;
	
		class EMFTGLFontInstance
		{
			public:
				EMFTGLFontInstance(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int d, const bool use_dl);
				~EMFTGLFontInstance();
		
				bool params_match(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int depth, const bool use_dl);
		
				FTFont* get_font() { return font; }
			private:
				EMFTGL::FontMode font_mode;
				string font_file_name;
				unsigned int face_size;
				unsigned int depth;
				bool use_display_lists;
				FTFont* font;
		};

		class EMFTGLManager
		{
			public:
				EMFTGLManager();
// 		static EMFTGLManager& instance();
	
				~EMFTGLManager();
		
				FTFont* get_font(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int d, const bool use_dl);
			private:
		
		
				vector<EMFTGLFontInstance*> font_instances;
		};
		
		EMFTGLManager fm;
};
	



} // namespace EMAN

#endif //emftgl_h__