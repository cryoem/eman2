/**
 * $Id$
 */
#ifndef eman__log__h__
#define eman__log__h__ 1

#include <stdio.h>
#include <stdarg.h>
#include <string>

using std::string;

namespace EMAN
{
    /** todo:
     * redesign begin() and end() methods
     */
    
    /** 
     * class Log defines a way to output logging information.
     * 1) The logs can either go to standard output (default), or go to a
     *    user-given file. 
     * 2) 4 verbose log levels are defined. by default, ERROR_LOG is used.
     * 3) Typical usage:
     *    Log::logger()->set_level(Log::WARNING_LEVEL);
     *    Log::logger()->error("cannot open file");
     */
    class Log
    {
    public:
	enum LogLevel {
	    ERROR_LOG,
	    WARNING_LOG,
	    NORMAL_LOG,
	    VARIABLE_LOG
	};

    public:
	static Log *logger();

	int begin(int argc, char *argv[], int ppid);
	void end(int ref, char *file = "-", char *text = "");

	
	void variable(const char *format, ...);
	void log(const char *format, ...);
	void warn(const char *format, ...);
	void error(const char *format, ...);

	void set_level(LogLevel level);
	void set_logfile(const char *filename);

    private:
	Log();
	Log(const Log &);
	~Log();

	void vlog(const char *format, LogLevel level, va_list arg);
	
	static Log *instance;
	FILE *out;
	LogLevel log_level;
	string default_emandir;
	string default_emanlog;
    };
}

#endif
