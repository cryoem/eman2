#ifndef eman__exception_h__
#define eman__exception_h__ 1

#include <string>

using std::string;

namespace EMAN {
    
    class Exception : public std::exception {
    public:
		string filename;
		int line;
		string desc;
		string objname;
		
		Exception(const string& file = "", int line_num = 0,
				  const string& desc1 = "", const string& objname1 = "")
			: filename(file), line(line_num), desc(desc1), objname(objname1) {}
	
		virtual ~Exception() throw() {}

		virtual const char *what() const throw() ;
    };


	class _NotExistingObjectException : public Exception
	{
	public:
		_NotExistingObjectException(const string& objname, const string& file = "unknown",
									int line = 0, const string& desc = "none")
			: Exception(file, line, desc, objname) {}
	};
	
#define NotExistingObjectException(objname, desc)  _NotExistingObjectException(objname, __FILE__, __LINE__, desc)

	
	class _ImageFormatException : public Exception
	{
	public:
		_ImageFormatException(const string& desc, const string& file = "unknown",
							 int line = 0)
			: Exception(file, line, desc) {}
	};

#define ImageFormatException(desc) _ImageFormatException(desc, __FILE__, __LINE__)

	class _ImageReadException : public Exception
	{
	public:
		_ImageReadException(const string& imagename, const string& file = "unknown",
						   int line = 0, const string& desc = "")
			: Exception(file, line, desc, imagename) {}
	};
#define ImageReadException(imagename, desc) _ImageReadException(imagename, __FILE__, __LINE__, desc)
	
	
	class _ImageWriteException : public Exception
	{
	public:
		_ImageWriteException(const string& imagename, const string& file = "unknown",
							int line = 0, const string& desc = "")
			: Exception(file, line, desc, imagename) {}
	};
#define ImageWriteException(imagename, desc) _ImageWriteException(imagename, __FILE__, __LINE__, desc)

	class _NullPointerException : public Exception
	{
	public:
		_NullPointerException(const string& file = "unknown",
							 int line = 0, const string& desc = "")
			: Exception(file, line, desc) {}
	};
#define NullPointerException(desc) _NullPointerException(__FILE__, __LINE__, desc)
	
	
	class ShrinkFactorException : public Exception
	{
	public:
		ShrinkFactorException(int sf, const string& file = "unknown",
							  int line = 0, const string& desc = "")
			: Exception(file, line, desc)
		{
			char s[32];
			sprintf(s, "%d", sf);
			objname = string(s);
		}
	};

	
	class InvalidModeException : public Exception
	{
	public:
		InvalidModeException(int low_mode, int high_mode, int wrong_mode,
							 const string& file = "unknown",
							 int line = 0, const string& desc = "")
			: Exception(file, line, desc)
		{
			char s[128];
			sprintf(s, "mode %d out of range [%d,%d]",
					wrong_mode, low_mode, high_mode);
			objname = string(s);
		}
	};
					

}
    
#endif
