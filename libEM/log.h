#ifndef __log__h__
#define __log__h__

#include <stdio.h>
#include <stdarg.h>

namespace EMAN {
    class Log {
    public:
	enum LogLevel { ERROR_LOG, WARNING_LOG, NORMAL_LOG, VARIABLE_LOG };
	
    public:
	static Log* logger();

	void variable(const char* format, ...);
	void log(const char* format, ...);
	void warn(const char* format, ...);
	void error(const char* format, ...);
	
	void set_level(LogLevel level);
	void set_logfile(const char* filename);
	
    private:
	Log();
	Log(const Log&);
	~Log();
	
	void vlog(const char* format, LogLevel level, va_list arg);
	
	static Log* instance;
	FILE* out;
	LogLevel log_level;
    };
}

#endif
