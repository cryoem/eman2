/**
 * $Id$
 */
#include "emcache.h"
#include "imageio.h"
#include <assert.h>

using namespace EMAN;

GlobalCache *GlobalCache::global_cache = 0;

GlobalCache::GlobalCache()
{
    imageio_cache = new EMCache<ImageIO>(8);
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



ImageIO *GlobalCache::get_imageio(string filename)
{
    return imageio_cache->get(filename);
}


void GlobalCache::add_imageio(string filename, ImageIO * io)
{
    if (io) {
	imageio_cache->add(filename, io);
    }
}
