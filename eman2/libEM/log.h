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

#ifndef eman__log__h__
#define eman__log__h__ 1

#include <cstdarg>
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
		void end(int ref, const string& file = "-", const string& text ="");

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
