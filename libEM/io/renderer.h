#ifndef EMAN_RENDERER_H
#define EMAN_RENDERER_H
#include <math.h>

namespace EMAN {

	template <EMUtil::EMDataType I>
	struct EM2Type {};

	template <>
	struct EM2Type<EMUtil::EM_CHAR> {
		using type = char;
	};

	template <>
	struct EM2Type<EMUtil::EM_UCHAR> {
		using type = unsigned char;
	};

	template <>
	struct EM2Type<EMUtil::EM_SHORT> {
		using type = short;
	};

	template <>
	struct EM2Type<EMUtil::EM_USHORT> {
		using type = unsigned short;
	};

	template <>
	struct EM2Type<EMUtil::EM_INT> {
		using type = int;
	};

	template <>
	struct EM2Type<EMUtil::EM_UINT> {
		using type = unsigned int;
	};

	template <>
	struct EM2Type<EMUtil::EM_FLOAT> {
		using type = float;
	};

	template <>
	struct EM2Type<EMUtil::EM_DOUBLE> {
		using type = double;
	};


	class Renderer {

	protected:
		Renderer() =default;

		int renderbits = 16;
		int renderlevel = 1;
		float rendermax = 0.0;
		float rendermin = 0.0;

		template<class T>
		auto getRenderedDataAndRendertrunc(float *data, size_t size);
	};

	template<class T>
	auto Renderer::getRenderedDataAndRendertrunc(float *data, size_t size) {
		if constexpr (!std::is_same<T, float>::value) {
			float RMIN;
			float RMAX = (1 << renderbits) - 1;

			if constexpr(std::is_unsigned<T>::value)
				RMIN = 0.0f;
			else
				RMIN = -(1 << (renderbits - 1));

			std::vector<T> rendered_data(size);
			size_t count = 0;

			for (size_t i = 0; i < size; ++i) {
				if (data[i] <= rendermin) {
					rendered_data[i] = T(RMIN);
					count++;
				}
				else if (data[i] >= rendermax) {
					rendered_data[i] = T(RMAX);
					count++;
				}
				else
					rendered_data[i] = (T)roundf(  (data[i] - rendermin)
					                             / (rendermax - rendermin)
					                             * (RMAX - RMIN) + RMIN);
			}
			return std::tuple(rendered_data, count);
		}
		else
			return std::tuple(std::vector<float>(data, data + size), (size_t)0);
	}

}


#endif //EMAN_RENDERER_H
