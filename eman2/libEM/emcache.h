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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */
 

#ifdef IMAGEIO_CACHE    
#ifndef eman__emcache__h__
#define eman__emcache__h__ 1

#include <cstdlib>
#include <string>
#include <map>
using std::string;
using std::map;

namespace EMAN
{
	class ImageIO;

	/* GlobalCache is a Singleton class that handles cache across EMAN. */
	class GlobalCache
	{
	  public:
		static GlobalCache *instance();
		ImageIO *get_imageio(const string & filename, int rw);
        int contains(const string & filename);
		void add_imageio(const string & filename, int rw, int persist, ImageIO * io);
        void close_imageio(const string & filename);
        void delete_imageio(const string & filename);
        void clean();

	  private:
        pthread_mutex_t mutex;          
		static GlobalCache *global_cache;
		map < string, ImageIO* >file_imageio;
		map < string, int >file_rw;
        map < string, int >file_ref;
        map < string, int >file_time;
        map < string, int >file_persist;

		GlobalCache();
		~GlobalCache();

	};

}

#endif // Header check
#endif // IMAGEIO_CACHE
