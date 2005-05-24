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
#if 1
#ifndef __func__
#define __func__ ""
#endif
#endif
	
#define ENTERFUNC  LOGDEBUG("Enter ")
#define EXITFUNC   LOGDEBUG("Exit ")
	
#define LOGERR Log::logger()->loc(Log::ERROR_LOG, __FILE__, __LINE__, __func__); Log::logger()->error

#define LOGWARN Log::logger()->loc(Log::WARNING_LOG, __FILE__, __LINE__, __func__); Log::logger()->warn

#define LOGDEBUG Log::logger()->loc(Log::DEBUG_LOG, __FILE__, __LINE__, __func__); Log::logger()->debug

#define LOGVAR Log::logger()->loc(Log::VARIABLE_LOG, __FILE__, __LINE__, __func__); Log::logger()->variable
	

	/** Log defines a way to output logging information.
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
		enum LogLevel
		{
			ERROR_LOG,        // error message
			WARNING_LOG,      // warning message
			DEBUG_LOG,        // debug message, usually at function level.
			VARIABLE_LOG      // very-detailed-level debug message
		};

	  public:
		static Log *logger();

		/** begin() and start() are used by command-line programs*/
		int begin(int argc, char *argv[], int ppid);
		void end(int ref, char *file = "-", char *text = "");

		/** log an error message. log level = ERROR_LOG.
		 * Its args are the same as printf().*/
		void error(const char *format, ...);

		/** log a warning message. log level = WARNING_LOG.
		 * Its args are the same as printf().*/
		void warn(const char *format, ...);

		/** log a debug message. log level = DEBUG_LOG.
		 * Its args are the same as printf().*/
		void debug(const char *format, ...);

		/** log a very-detailed-level debug message. log level = VARIABLE_LOG.
		 * Its args are the same as printf().*/
		void variable(const char *format, ...);
		
		void set_level(int level);

		/** set log output file. If this function is not called,
		 * output is standart output.
		 */
		void set_logfile(const char *filename);

		void loc(LogLevel level, const string & file, int linenum, const string & func);
		
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
		string location;
	};
}

#endif
