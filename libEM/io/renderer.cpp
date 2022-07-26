#include "renderer.h"

using EMAN::EMUtil;
using EMAN::Renderer;


EMUtil::EMDataType Renderer::rendered_dt(EMUtil::EMDataType dt,
                                         std::initializer_list<decltype(dt)> dts) {
	if(dt == EMUtil::EM_COMPRESSED) {
		map<int, decltype(dt)> bits;
		for (const auto &l: dts) {
			int key;
			switch (l) {
				case EMUtil::EM_FLOAT:
					key = 0;
					break;
				case EMUtil::EM_CHAR:
				case EMUtil::EM_SHORT:
					key = 2*EMDataTypeBits[l] + 1;
					break;
				default:
					key = 2*EMDataTypeBits[l];
			}
			bits[key] = l;
		}

		renderbits = std::min(renderbits, EMDataTypeBits[(int)rbegin(bits)->second]);

		auto it = bits.lower_bound(renderbits*2);
		dt = it->second;
	}

	return dt;
}
