#include "renderer.h"

using EMAN::EMUtil;
using EMAN::Renderer;


EMUtil::EMDataType Renderer::rendered_dt(EMUtil::EMDataType dt) const {
	if(dt == EMUtil::EM_COMPRESSED) {
		if (renderbits <= 0)       dt = EMUtil::EM_FLOAT;
		else if (renderbits <= 8)  dt = EMUtil::EM_UCHAR;
		else if (renderbits <= 16) dt = EMUtil::EM_USHORT;
	}

	return dt;
}
