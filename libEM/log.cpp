/**
 * $Id$
 */
#include "log.h"
#include "util.h"
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef WIN32
#include <windows.h>
#include <io.h>
#include <time.h>
#include <process.h>
#include <direct.h>
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

	location = Util::sbasename(filename) + ":" + Util::int2str(linenum) + " " + func;
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


void Log::set_level(LogLevel level)
{
	log_level = level;
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
	int ref = getpid();

	string filename = Util::sbasename(argv[0]);
	
	char s[4048];
#ifndef WIN32
	sprintf(s, "%d\t%d\t%d\t%d\t%s", ref, tm, 0, ppid ? ppid : getppid(), filename.c_str());
#else
	sprintf(s, "%d\t%d\t%d\t%d\t%s", ref, tm, 0, ppid, filename.c_str());
#endif
	for (int i = 1; i < argc; i++) {
		sprintf(s + strlen(s), " %s", argv[i]);
	}
	sprintf(s + strlen(s), "\n");

	FILE *eman_file = fopen(default_emanlog.c_str(), "a");
	if (!eman_file) {
		return 0;
	}

	Util::file_lock_wait(eman_file);
	fprintf(eman_file, "%s", s);
	Util::file_unlock(eman_file);
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


void Log::end(int ref, char *file, char *text)
{
	FILE *out = fopen(".emanlog", "a");

	if (out) {
		time_t tm = time(0);
		Util::file_lock_wait(out);
		fprintf(out, "%d\t%d\t%s\t%s\n", ref, tm, file, text);
		Util::file_unlock(out);
		fclose(out);
	}
}
