#ifndef EMAN_RENDERER_H
#define EMAN_RENDERER_H


namespace EMAN {

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
		float RMIN;
		float RMAX = (1 << renderbits) - 1;

		if constexpr(std::is_unsigned<T>::value)
			RMIN = 0.0f;
		else
			RMIN = -(1 << (renderbits - 1));

		auto rendered_data = new T[size];
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

}


#endif //EMAN_RENDERER_H
