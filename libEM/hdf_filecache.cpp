/**
 * $Id$
 */

/*
 * Author: Grant Tang, 05/01/2012 (gtang@bcm.edu)
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

#ifdef HDFIO_CACHE

#include <list>
#include "hdf_filecache.h"
#include <boost/date_time.hpp>
#include <sys/time.h>
#include <sys/resource.h>

using namespace EMAN;
using std::list;

FileItem::FileItem() :
		_path(""), _imgio(0), _timestamp(0), _readonly(true)
{
}

FileItem::FileItem(const string& fpath, ImageIO* const fimageio, const time_t& ftimestamp, bool freadonly) :
		_path(fpath), _imgio(fimageio), _timestamp(ftimestamp), _readonly(freadonly)
{
}

FileItem::~FileItem()
{
	delete _imgio;
}

int FileItem::set_path(const string& fpath)
{
	_path = fpath;
	return 0;
}

string FileItem::get_path() const
{
	return _path;
}

int FileItem::set_imgio(ImageIO* const fimageio)
{
	_imgio = fimageio;
	return 0;
}

ImageIO * FileItem::get_imgio() const
{
	return _imgio;
}

int FileItem::set_timestamp(const time_t& ftimestamp)
{
	_timestamp = ftimestamp;
	return 0;
}

time_t FileItem::get_timestamp() const
{
	return _timestamp;
}

int FileItem::set_readonly(bool freadonly)
{
	_readonly = freadonly;
	return 0;
}

bool FileItem::get_readonly() const
{
	return _readonly;
}


HDFCache * HDFCache::_instance = 0;

HDFCache::HDFCache()
{
	_thread = boost::thread(&HDFCache::cache_routine, this);

	struct rlimit rl;
	getrlimit(RLIMIT_NOFILE, &rl);
	CACHESIZE = rl.rlim_cur/2;
//	std::cout << "CACHESIZE = " << CACHESIZE << std::endl;
}

HDFCache * HDFCache::instance()
{
	if(!_instance) {
		_instance = new HDFCache();
	}

	return _instance;
}

void HDFCache::cache_routine()
{
	boost::posix_time::seconds intervalTime(CHECK_INTERVAL);

	while(true) {
		try {
			force_clean();
			purge_file();
			boost::this_thread::sleep(intervalTime);
		}
		catch (std::exception& exc) {
			std::cerr << "Very Bad Thing Happen in HDFCache routine: " << exc.what() << std::endl;
		}
	}
}

FileItem * HDFCache::get_file(const string& filename)
{
	boost::mutex::scoped_lock l(m_mutex);

	if(file_pool.find(filename) != file_pool.end()) {
		FileItem * hdf_item = file_pool[filename];
		hdf_item->set_timestamp(time(0));
		return hdf_item;
	}
	else {
		return 0;
	}
}

int HDFCache::add_file(FileItem * newfile)
{
	if(file_pool.size() >= CACHESIZE) {
		force_clean();
	}

	boost::mutex::scoped_lock l(m_mutex);
	file_pool[newfile->get_path()] = newfile;
	file_pool2.push_back(newfile);

	return 0;
}

/**I don't put lock in this function to avoid race condition.
 * the lock is in upper level function purge_file() and force_clean()*/
int HDFCache::close_file(const string& filename)
{
	FileItem * hdfitem = file_pool[filename];

	ImageIO * hdfio = hdfitem->get_imgio();
	delete hdfio;
	file_pool[filename]->set_imgio(0);
	file_pool.erase(filename);

	vector<FileItem *>::iterator it = find(file_pool2.begin(), file_pool2.end(), hdfitem);
	file_pool2.erase(it);

	return 0;
}

int HDFCache::purge_file()
{
	boost::mutex::scoped_lock l(m_mutex);

	//close those files have access time over the threshold
	sort(file_pool2.begin(), file_pool2.end(), least_access());

	vector<FileItem *>::iterator it2 = lower_bound(file_pool2.begin(), file_pool2.end(), ACCESS_TIME_THRESHOLD, access_time_cmp());

	if(it2 != file_pool2.begin()) {
		list<FileItem *> expired_file;
		copy(file_pool2.begin(), it2, back_inserter(expired_file));
		list<FileItem *>::const_iterator lit;
		for(lit = expired_file.begin(); lit!=expired_file.end(); ++lit) {
			string filename = (*lit)->get_path();
			close_file(filename);
		}
	}

	return 0;
}

int HDFCache::force_clean()
{
	boost::mutex::scoped_lock l(m_mutex);

	//close 1/10 oldest files when cache size reach limit
	if(file_pool2.size() >= CACHESIZE) {
		sort(file_pool2.begin(), file_pool2.end(), least_access());
		for(size_t i=0; i<=CACHESIZE/10; ++i) {
			close_file(file_pool2.front()->get_path());
		}
	}

	return 0;
}

#endif	//HDFIO_CACHE
