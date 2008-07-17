
#ifdef EMAN2_USING_FTGL

#include "emftgl.h"
using namespace EMAN;

#include "exception.h"
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <algorithm>

//static init
// string  EMFTGL::font_file_name = "/usr/share/fonts/liberation/LiberationSans-Regular.ttf";
// unsigned int EMFTGL::face_size = 32;
// unsigned int EMFTGL::depth = 32;
// EMFTGL::FontMode EMFTGL::mode = EMFTGL::TEXTURE;
// EMFTGL::EMFTGLManager EMFTGL::fm;

void EMFTGL::render_string(const string& message) {
	
	FTFont* font = fm.get_font(mode,font_file_name,face_size,depth,true);
	
	if (font == 0) {
		cerr << "Couldn't open font, no action taken. Current font is " << font_file_name << endl;
		return;
	}
	
	font->Render(message.c_str());
}

vector<float> EMFTGL::bounding_box(const string& message)
{
	FTFont* font = fm.get_font(mode,font_file_name,face_size,depth,true);
	if (font == 0) {
		cerr << "Couldn't open font, no action taken. Current font is " << font_file_name << endl;
		return vector<float>(); 
	}
	vector<float> bounds(6);
	font->BBox(message.c_str(),bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
	return bounds;
}


// EMFTGL::EMFTGLManager& EMFTGL::EMFTGLManager::instance() {
// 	static EMFTGLManager fm;
// 	return fm;
// }

EMFTGL::EMFTGLManager::EMFTGLManager() : font_instances() { }
	

EMFTGL::EMFTGLManager::~EMFTGLManager() {
	
	for (vector<EMFTGLFontInstance*>::iterator it = font_instances.begin(); it != font_instances.end(); ++it ) {
		delete (*it);
		*it = 0;
	}
	font_instances.clear();
}

FTFont* EMFTGL::EMFTGLManager::get_font(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int depth, const bool use_dl)
{
	for (vector<EMFTGLFontInstance*>::const_iterator it = font_instances.begin(); it != font_instances.end(); ++it ) {
		if ((*it)->params_match(mode,file_name,face_size,depth,use_dl)) {
			cout << "found match!" << endl;
			return (*it)->get_font();
		}
	}
	// If we make it here there was no match
	cout << "constructing new font" << endl;
	EMFTGLFontInstance* fi = new EMFTGLFontInstance(mode,file_name,face_size,depth,use_dl);
	font_instances.push_back(fi);
	return fi->get_font();
}

EMFTGL::EMFTGLFontInstance::~EMFTGLFontInstance(){
	if (font != 0) {
		delete font;
		font = 0;
	}
}

EMFTGL::EMFTGLFontInstance::EMFTGLFontInstance(EMFTGL::FontMode mode, const string& file_name, const unsigned int _face_size, const unsigned int _depth, const bool use_dl) :
		font_mode(mode), font_file_name(file_name), face_size(_face_size), depth(_depth),use_display_lists(use_dl), font(0)
{
	if (mode == EMFTGL::PIXMAP) {
		font = new FTGLPixmapFont(font_file_name.c_str());
	}
	else if (mode == EMFTGL::TEXTURE) {
		font = new FTGLTextureFont(font_file_name.c_str());
	}
	else if ( mode == EMFTGL::EXTRUDE ) {
		font = new FTGLExtrdFont(font_file_name.c_str());
		font->Depth(depth);
	}
	else {
		LOGERR("Error, unsupported mode ");
		return;
	}
	
	if (font->Error()) {
		delete font;
		LOGERR( string("Could not open font file " + font_file_name).c_str());
		font = 0;
	}
	else {
		font->UseDisplayList(use_display_lists);
		font->FaceSize(face_size);
	}
}

bool EMFTGL::EMFTGLFontInstance::params_match(EMFTGL::FontMode mode, const string& file_name, const unsigned int face_size, const unsigned int depth, const bool use_dl) {
	if ( this->font_mode == mode && this->font_file_name == file_name &&
			this->face_size == face_size && this->depth == depth &&	this->use_display_lists == use_dl ) return true;
	
	return false;
}



#endif //EMAN2_USING_FTGL