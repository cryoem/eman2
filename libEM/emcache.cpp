/**
 * $Id$
 */
#include "emcache.h"
#include "imageio.h"
#include <assert.h>

#ifdef WIN32
#define access _access
#define F_OK 00
#else
#include <unistd.h>
#endif


using namespace EMAN;

GlobalCache *GlobalCache::global_cache = 0;

GlobalCache::GlobalCache()
{
	imageio_cache = new EMCache < ImageIO > (8);
}

GlobalCache::GlobalCache(const GlobalCache &)
{
}

GlobalCache::~GlobalCache()
{
	delete imageio_cache;
	imageio_cache = 0;
}


GlobalCache *GlobalCache::instance()
{
	if (!global_cache) {
		global_cache = new GlobalCache();
	}
	return global_cache;
}



ImageIO *GlobalCache::get_imageio(string filename, int rw_mode)
{
	ImageIO *io = imageio_cache->get(filename);
	if (io) {
		if (rw_mode != ImageIO::READ_ONLY) {
			imageio_cache->remove(filename);
			io = 0;
		}		
		else {
			if (access(filename.c_str(), F_OK) != 0) {
				printf("%s was removed\n", filename.c_str());
				imageio_cache->remove(filename);
				io = 0;
			}
		}
	}

	return io;
}


void GlobalCache::add_imageio(string filename, int rw_mode, ImageIO * io)
{
	if (io) {
		file_rw_dict[filename] = rw_mode;
		imageio_cache->add(filename, io);
	}
}

