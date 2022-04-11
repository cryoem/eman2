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
	};

}


#endif //EMAN_RENDERER_H
