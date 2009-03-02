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
 
#ifndef eman__emcache__h__
#define eman__emcache__h__ 1

#include <cstdlib>
#include <string>
#include <map>

using std::string;
using std::map;

namespace EMAN
{
	/** EMCache is a generic cache that can cache anything defined by 'T'.
     * Item in the cache is identified by its name. Use add() to add
     * an item to the cache, and get() to return the item. If the item
     * is not in the cache, a 0 pointer will be returned.
     */
	template < class T > class EMCache
	{
	  public:
		explicit EMCache(int cache_size)
		{
			item_cache = new T *[cache_size];
			name_cache = new string[cache_size];

			size = cache_size;
			nitems = 0;
		}

		~EMCache()
		{
			for (int i = 0; i < nitems; i++) {
				if( item_cache[i] )
				{
					delete item_cache[i];
					item_cache[i] = 0;
				}
			}

			if( item_cache )
			{
				delete[]item_cache;
				item_cache = 0;
			}

			if( name_cache )
			{
				delete[]name_cache;
				name_cache = 0;
			}
		}


		T *get(const string & itemname) const
		{
			T *result = 0;

			for (int i = 0; i < nitems; i++)
			{
				if (name_cache[i] == itemname) {
					result = item_cache[i];
					break;
				}
			}

			return result;
		}

		void add(const string & itemname, T * item)
		{
			if (!item) {
				return;
			}

			if (nitems < size) {
				item_cache[nitems] = item;
				name_cache[nitems] = itemname;
				nitems++;
			}
			else {
				int r = (int) (1.0 * size * rand() / (RAND_MAX + 1.0));
				if( item_cache[r] )
				{
					delete item_cache[r];
					item_cache[r] = 0;
				}

				item_cache[r] = item;
				name_cache[r] = itemname;
			}
		}

		void remove(const string & itemname)
		{
			int r = -1;
			for (int i = 0; i < nitems; i++) {
				if (name_cache[i] == itemname) {
					r = i;
					break;
				}
			}
			if (r >= 0) {
				if( item_cache[r] )
				{
					delete item_cache[r];
					item_cache[r] = 0;
				}
				name_cache[r] = "";
			}
		}

		int get_size() const
		{
			return size;
		}

	  private:
		T ** item_cache;
		string *name_cache;

		int size;
		int nitems;
	};
	
	class ImageIO;

	/** GlobalCache is a Singleton class that handles cache across EMAN.
     * It uses EMCache template in its internal implementation.
     */
	class GlobalCache
	{
	  public:
		static GlobalCache *instance();

		ImageIO *get_imageio(const string & filename, int rw_mode);
		void add_imageio(const string & filename, int rw_mode, ImageIO * io);

	  private:
		EMCache < ImageIO > *imageio_cache;
		static GlobalCache *global_cache;
		map < string, int >file_rw_dict;


		GlobalCache();
		GlobalCache(const GlobalCache & gc);
		 ~GlobalCache();

	};

}

#endif
