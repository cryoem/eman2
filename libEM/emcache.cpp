/**
 * $Id$
 */
 
/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
 
#include "emcache.h"
#include "imageio.h"
#include "util.h"

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
	if(imageio_cache)
	{
		delete imageio_cache;
		imageio_cache = 0;
	}
}


GlobalCache *GlobalCache::instance()
{
	if (!global_cache) {
		global_cache = new GlobalCache();
	}
	return global_cache;
}

ImageIO *GlobalCache::get_imageio(const string & filename, int rw_mode)
{
    ImageIO *io = imageio_cache->get(filename);
    if (io) {
        bool need_remove = false;
        
        int old_rw = file_rw_dict[filename];
        
        if (rw_mode == ImageIO::READ_ONLY) {
            if (old_rw == ImageIO::WRITE_ONLY) {
                need_remove = true;
            }
        }
        else if (rw_mode != old_rw) {
            need_remove = true;
        }
	else {
	    if (!Util::is_file_exist(filename)) {
		need_remove = true;
            }
        }
        if (need_remove) {
            imageio_cache->remove(filename);
            io = 0;
        }        
    }
    
    return io;
}


void GlobalCache::add_imageio(const string & filename, int rw_mode, ImageIO * io)
{
	if (io) {
		file_rw_dict[filename] = rw_mode;
		imageio_cache->add(filename, io);
	}
}

