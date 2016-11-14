/**
 * $Id$
 */
 
/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu), Ian Rees 06/30/2013
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
// #define DEBUG_CACHE	1
 
#include "emcache.h"
#include "imageio.h"
#include "util.h"
#include <stdio.h>
#include <pthread.h>
using namespace EMAN;

void *thread_clean( void *ptr ) {
    while (1) {
        sleep(1);
        GlobalCache::instance()->clean();
    }
}
void thread_start() {
    pthread_t thread1;
    pthread_create(&thread1,NULL,&thread_clean,NULL);
}


// GlobalCache
GlobalCache *GlobalCache::global_cache = 0;

GlobalCache::GlobalCache()
{
    // pthread_mutex_init(&mutex, NULL);
    thread_start();
}

GlobalCache::~GlobalCache()
{
    clean();
    pthread_mutex_destroy(&mutex);
}

GlobalCache *GlobalCache::instance()
{
	if (!global_cache) {
		global_cache = new GlobalCache();
	}
	return global_cache;
}

ImageIO *GlobalCache::get_imageio(const string & filename, int rw)
{
    pthread_mutex_lock(&mutex);
#ifdef DEBUG_CACHE
    printf("get_imageio: filename: %s, rw: %d\n", filename.c_str(), rw);
#endif

    if (file_imageio.count(filename) == 0) {
        // printf("get_imageio: no cached imageio!\n");
        pthread_mutex_unlock(&mutex);
        return 0;        
    }     
    
    ImageIO *io = file_imageio[filename];    
#ifdef DEBUG_CACHE
    printf("get_imageio:      found cached; refs: %d, current rw: %d, request rw: %d\n", file_ref[filename], file_rw[filename], rw);
#endif
    time_t t = time(0);
    file_time[filename] = time(0);
    // printf("get_imageio:      updated time\n");
    if (file_rw[filename] == rw) {
        // current == request
        file_ref[filename]++;
        // printf("get_imageio:      ref += 1 to %d\n", file_ref[filename]);
    } else if (file_rw[filename] == 2 && rw == 1) {
        // current == read_write, request == read
        file_ref[filename]++;
        // printf("get_imageio:      ref += 1 to %d (passing write mode to read request)\n", file_ref[filename]);
    } else if (file_ref[filename] == 0) {
        // current != request, no open handles
        // printf("get_imageio:      re-opening in mode %d\n", rw);
        this->delete_imageio(filename);
        io = 0;        
    } else {
        // current != request, open handles
        // current is read, request is write -- we can't make this change with open handles. How to fail?
        file_ref[filename]++;
        printf("GlobalCache::get_imageio: %s: Cannot switch mode %d to %d with open references!\n", filename.c_str(), file_rw[filename], rw);
    }
    pthread_mutex_unlock(&mutex);
    return io;
}

void GlobalCache::delete_imageio(const string & filename) {
    // Lock acquired in caller
    if (file_imageio.count(filename) == 0) {
        // printf("delete_imageio: Not in cache\n");
    } else if (file_ref[filename] == 0) {
        ImageIO *io = file_imageio[filename];
        file_imageio.erase(filename);
        delete io;
    }
}

void GlobalCache::close_imageio(const string & filename) {
    pthread_mutex_lock(&mutex);
#ifdef DEBUG_CACHE
    printf("close_imageio: filename: %s\n", filename.c_str());
#endif
    if (file_imageio.count(filename) > 0) {
        ImageIO *io = file_imageio[filename];
        // Decrement the reference and update the access time.
        file_ref[filename]--;
        file_time[filename] = time(0);
        // printf("close_imageio:      ref -= 1 to %d\n", file_ref[filename]);
        if (file_ref[filename] == 0 && file_persist[filename] == 0) {
            // printf("close_imageio:      auto closing!\n");
            this->delete_imageio(filename);
        }
    } else {
        // printf("close_imageio:      No ImageIO to close!\n");
    }
    pthread_mutex_unlock(&mutex);
}

void GlobalCache::add_imageio(const string & filename, int rw, int persist, ImageIO * io)
{
    pthread_mutex_lock(&mutex);    
#ifdef DEBUG_CACHE
    printf("add_imageio: filename %s, rw %d, persist %d\n", filename.c_str(), rw, persist);
#endif
    if (file_imageio.count(filename) > 0) {
        // Todo: How to handle error if file is already in cache?
        // printf("File already in cache.")
    } else if (io && persist > 0) {
        // Todo: Make a class to hold all these values and the ImageIO.
        file_imageio[filename] = io;
        file_ref[filename] = 1;
		file_rw[filename] = rw;
        file_time[filename] = time(0);
        file_persist[filename] = persist;
	}
    pthread_mutex_unlock(&mutex);
}

int GlobalCache::contains(const string & filename) {
    pthread_mutex_lock(&mutex);
    int ret = 0;
    if (file_imageio.count(filename) > 0) {
        ret = 1;
    }
    pthread_mutex_unlock(&mutex);    
    return ret;
}

void GlobalCache::clean() {
    pthread_mutex_lock(&mutex);    
#ifdef DEBUG_CACHE
    printf("clean: %ld items in cache\n", file_imageio.size());
#endif
    
    // Find cached items to close
    vector<string> toclose;
    time_t now = time(0);
    int seconds;
    string filename;
    for (map<string, ImageIO*>::iterator iter=file_imageio.begin(); iter!=file_imageio.end(); ++iter) {
        filename = iter->first;
        seconds = int(difftime(now, file_time[filename]));
#ifdef DEBUG_CACHE
        printf("clean:       filename: %s, rw: %d, ref: %d, last access: %d, persist: %d\n", filename.c_str(), file_rw[filename], file_ref[filename], seconds, file_persist[filename]);
#endif
        if (seconds >= file_persist[filename]) {
            // printf("clean:      marking to close %s...\n", filename.c_str());
            toclose.push_back(filename);
        }
    }
    
    // Close the selected items
    for(vector<string>::iterator iter=toclose.begin(); iter!=toclose.end(); iter++) {
        filename = *iter;
        // printf("clean: closing %s...\n", filename.c_str());        
        this->delete_imageio(filename);
    }
    
    pthread_mutex_unlock(&mutex);
}

#endif // IMAGEIO_CACHE

