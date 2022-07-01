#include "renderer.h"

using EMAN::EMUtil;
using EMAN::Renderer;


EMUtil::EMDataType Renderer::rendered_dt(EMUtil::EMDataType dt,
                                         std::initializer_list<decltype(dt)> dts) const {
	if(dt == EMUtil::EM_COMPRESSED) {
		map<int, decltype(dt)> bits;
		for (const auto &l: dts) {
			int key = (l == EMUtil::EM_FLOAT ? 0 : EMDataTypeBits[l]);
			bits[key] = l;
		}

		dt = bits.lower_bound(renderbits)->second;
	}

	return dt;
}
