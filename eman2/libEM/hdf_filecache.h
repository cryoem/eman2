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

#ifndef eman__hdf_filecache__h__
#define eman__hdf_filecache__h__ 1

#ifdef	HDFIO_CACHE

#include <string>
#include <cstdio>
#include <ctime>
#include <functional>

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>

#include "imageio.h"

using std::string;
using std::binary_function;

const static int CHECK_INTERVAL = 5;	//seconds for thread to check if the opened file need been closed
const static int ACCESS_TIME_THRESHOLD = 30;	//seconds since the last access to close HDF image file
const static int ACCESS_WRITE_TIME_THRESHOLD = 5;	//seconds since the last access to close HDF image file

namespace EMAN
{

/**A simple class wrap the ImageIO pointer. Although this file cache is used only for HDF5 file I/O,
 * I still use ImageIO as file handle for two reasons:
 * 1. This file cache may be used for other image format in the future.
 * 2. We have two kinds of HDF I/O object due to the internal format change.
 * */
class FileItem {
public:
	FileItem();
	FileItem(const string& fpath, ImageIO* const fimageio, const time_t& ftimestamp, bool freadonly);
	~FileItem();

	int set_path(const string& fpath);
	string get_path() const;

	int set_imgio(ImageIO* const fimageio);
	ImageIO * get_imgio() const;

	int set_timestamp(const time_t& ftimestamp);
	time_t get_timestamp() const;

	int set_readonly(bool freadonly);
	bool get_readonly() const;

private:
	string 		_path;		//full path name to the opened file
	ImageIO * 	_imgio;		//image I/O object to the opened file
	time_t		_timestamp;	//last access time
	bool		_readonly;	//Is the file is opened in read only mode
};

/**A file cache for HDF5 I/O object. A background thread will routinely check the file
 * cache to close those files have expired.*/
class HDFCache {
public:
	static HDFCache *instance();

	//file name is a full path name
	FileItem * get_file(const string& filename);

	int add_file(FileItem * newfile);

private:
	HDFCache();
	HDFCache(const HDFCache&);
	~HDFCache();

	static HDFCache * _instance;

	//Primary file handle container for fast accessing by file name
	typedef boost::unordered_map<string, FileItem *> MyHashtable;
	MyHashtable file_pool;

	//Secondary file handle container for fast sorting to close those files expired
	vector<FileItem *> file_pool2;

	boost::thread	_thread;
	boost::mutex 	m_mutex;

	size_t CACHESIZE;		//size limit for the HDF file cache, half of the file descriptor limit per process

	//These functions are running in the thread
	int purge_file();		//purge the file has last access time >= threshold
	int force_clean();		//close some files when the cache reach it's size limit
	void cache_routine();	//function for background thread
	int close_file(const string& filename);

	//functor for sorting file pointer by access time
	struct least_access : public binary_function<FileItem *, FileItem *, bool>{
		bool operator()(FileItem * x, FileItem * y) const {
			return x->get_timestamp() < y->get_timestamp();
		}
	};

	//functor for finding those files have access time over the threshold
	struct access_time_cmp : public binary_function<FileItem *, int, bool>{
		bool operator()(FileItem * f, int interval) const {
			return f->get_timestamp() < time(0)-ACCESS_TIME_THRESHOLD;
		}
	};
};

}

#endif	//HDFIO_CACHE

#endif	//eman__hdf_filecache__h__
