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

#include "log.h"
#include "util.h"
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>

#ifdef WIN32
#include <time.h>
#include <process.h>
#else
#include <unistd.h>
#endif

using namespace EMAN;
using std::string;

Log::Log()
{
	out = 0;
	log_level = ERROR_LOG;
#ifdef WIN32
	char c[2];
	c[0] = getenv("WINDIR")[0];
	c[1] = '\0';
	default_emandir = string(c) + string(":\\.eman");
#else
	default_emandir = string(getenv("HOME")) + "/.eman";
	mkdir(default_emandir.c_str(), 0xffff);
#endif
	default_emanlog = ".emanlog";
	location = "";
}

Log::Log(const Log &)
{
}

Log::~Log()
{
	if (out) {
		fclose(out);
		out = 0;
	}
}

Log *Log::instance = 0;

Log *Log::logger()
{
	if (!instance) {
		instance = new Log();
	}
	return instance;
}

void Log::loc(LogLevel level, const string & filename, int linenum, const string & func)
{
	if (log_level < level) {
		return;
	}

	location = Util::sbasename(filename) + ":" + Util::int2str(linenum);
	if (func != "") {
		location +=" " + func + "()";
	}
}

void Log::vlog(const char *format, LogLevel level, va_list arg)
{
	if (log_level < level) {
		return;
	}

	const char *key = "";

	switch (level) {
	case WARNING_LOG:
		key = "Warning: ";
		break;
	case ERROR_LOG:
		key = "Error: ";
		break;
	default:
		key = "";
	}

	FILE *file = stdout;
	if (out) {
		file = out;
	}

	fprintf(file, "%s", key);
	vfprintf(file, format, arg);
	if (location != "") {
		fprintf(file, " at %s", location.c_str());
	}
	fprintf(file, "\n");
}

void Log::variable(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	vlog(format, VARIABLE_LOG, arg);
	va_end(arg);
}

void Log::debug(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	vlog(format, DEBUG_LOG, arg);
	va_end(arg);
}

void Log::warn(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	vlog(format, WARNING_LOG, arg);
	va_end(arg);
}

void Log::error(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	vlog(format, ERROR_LOG, arg);
	va_end(arg);
}


void Log::set_level(int level)
{
	log_level = (LogLevel)level;
}

void Log::set_logfile(const char *filename)
{
	if (filename && !out) {
		out = fopen(filename, "wb");
	}
}

int Log::begin(int argc, char *argv[], int ppid)
{
	time_t tm = time(0);
	const char *pwd = getenv("PWD");
#ifdef	_WIN32
	int ref = _getpid();
#else
	int ref = getpid();
#endif

	string filename = Util::sbasename(argv[0]);

	char s[4048];
#ifndef WIN32
	sprintf(s, "%d\t%d\t%d\t%d\t%s", ref, (int)tm, 0, ppid ? ppid : getppid(), filename.c_str());
#else
	sprintf(s, "%d\t%d\t%d\t%d\t%s", ref, (int)tm, 0, ppid, filename.c_str());
#endif
	for (int i = 1; i < argc; i++) {
		sprintf(s + strlen(s), " %s", argv[i]);
	}
	sprintf(s + strlen(s), "\n");

	FILE *eman_file = fopen(default_emanlog.c_str(), "a");
	if (!eman_file) {
		return 0;
	}

	//Util::file_lock_wait(eman_file);
	fprintf(eman_file, "%s", s);
	fclose(eman_file);

	string dirlist = default_emandir + "./dirlist";
	FILE *in = fopen(dirlist.c_str(), "r");
	if (in) {
		char s[1024];
		int f = 0;
		while (fscanf(in, " %1023s", s) == 1) {
			if (strcmp(s, pwd) == 0) {
				f = 1;
				break;
			}
		}

		fclose(in);
		if (!f) {
			in = 0;
		}
	}

	if (!in) {
		FILE *dirout = fopen(dirlist.c_str(), "a");
		if (dirout) {
			fprintf(dirout, "%s\n", pwd);
			fclose(dirout);
		}
	}

	return ref;
}


void Log::end(int ref, const string& file, const string& text)
{
	FILE *out = fopen(".emanlog", "a");

	if (out) {
		time_t tm = time(0);
		//Util::file_lock_wait(out);
		fprintf(out, "%d\t%ld\t%s\t%s\n", ref, tm, file.c_str(), text.c_str());
		fclose(out);
	}
}
