#ifndef eman__exception_h__
#define eman__exception_h__ 1

#include <string>

using std::string;

namespace EMAN {
    
    class Exception : public std::exception {
    public:
		string filename;
		int linenum;
		string desc;
		string objname;
		
		Exception(const string& file = "", int line = 0,
				  const string& desc1 = "", const string& objname1 = "")
			: filename(file), linenum(line), desc(desc1), objname(objname1) {}
	
		virtual ~Exception() throw() {}

		virtual const char *what() const throw();

		virtual const char *name() const { return "Exception"; }
    };


	class _NotExistingObjectException : public Exception
	{
	public:
		_NotExistingObjectException(const string& objname1, const string& file = "unknown",
									int line = 0, const string& desc1 = "none")
			: Exception(file, line, desc1, objname1) {}
		
		const char *name() const { return "NotExistingObjectException"; }
		
	};
	
#define NotExistingObjectException(objname, desc)  _NotExistingObjectException(objname, __FILE__, __LINE__, desc)

	
	class _ImageFormatException : public Exception
	{
	public:
		_ImageFormatException(const string& desc1, const string& file = "unknown",
							  int line = 0)
			: Exception(file, line, desc1) {}
		
		const char *name() const { return "ImageFormatException"; }
		
	};

#define ImageFormatException(desc) _ImageFormatException(desc, __FILE__, __LINE__)

	
	class _ImageDimensionException : public Exception
	{
	public:
		_ImageDimensionException(const string& desc1, const string& file = "unknown",
								 int line = 0)
			: Exception(file, line, desc1) {}
		
		const char *name() const { return "ImageDimensionException"; }
		
	};

#define ImageDimensionException(desc) _ImageDimensionException(desc, __FILE__, __LINE__)

	class _FileAccessException : public Exception
	{
	public:
		_FileAccessException(const string& filename1, const string& file = "unknown",
							 int line = 0, const string& desc1 = "")
			: Exception(file, line, desc1, filename1)
		{
			desc = "cannot access file '" + filename1 + "'";
		}
		
		const char *name() const { return "FileAccessException"; }
		
	};
#define FileAccessException(filename) _FileAccessException(filename, __FILE__, __LINE__)
	
	
	class _ImageReadException : public Exception
	{
	public:
		_ImageReadException(const string& imagename, const string& file = "unknown",
							int line = 0, const string& desc1 = "")
			: Exception(file, line, desc1, imagename) {}
		
		const char *name() const { return "ImageReadException"; }
		
	};
#define ImageReadException(filename, desc) _ImageReadException(filename, __FILE__, __LINE__, desc)
	
	
	class _ImageWriteException : public Exception
	{
	public:
		_ImageWriteException(const string& imagename, const string& file = "unknown",
							 int line = 0, const string& desc1 = "")
			: Exception(file, line, desc1, imagename) {}
		
		const char *name() const { return "ImageWriteException"; }
		
	};
#define ImageWriteException(imagename, desc) _ImageWriteException(imagename, __FILE__, __LINE__, desc)

	class _NullPointerException : public Exception
	{
	public:
		_NullPointerException(const string& file = "unknown",
							  int line = 0, const string& desc1 = "")
			: Exception(file, line, desc1) {}
		
		const char *name() const { return "NullPointerException"; }
		
	};
#define NullPointerException(desc) _NullPointerException(__FILE__, __LINE__, desc)
	

	class _TypeException : public Exception
	{
	public:
		_TypeException(const string & desc1, const string & type,
					   const string & file = "unknown", int line = 0)
			: Exception(file, line, desc1, type) {}
		
		const char *name() const { return "TypeException"; }
		
	};
#define TypeException(desc, type) _TypeException(desc, type, __FILE__, __LINE__)

	
	class _InvalidValueException : public Exception
	{
	public:
		_InvalidValueException(int val, const string& file = "unknown",
							   int line = 0, const string& desc1 = "")
			: Exception(file, line, desc1)
		{
			char s[32];
			sprintf(s, "%d", val);
			objname = string(s);
		}		
		const char *name() const { return "InvalidValueException"; }
		
	};
#define InvalidValueException(val, desc) _InvalidValueException(val, __FILE__, __LINE__, desc)


	
	class _OutofRangeException : public Exception
	{
	public:
		_OutofRangeException(int low, int high, int input,
							 const string& file = "unknown", 
							 int line = 0, const string & desc1 = "",
							 const string& objname1 = "")
			: Exception(file, line, desc1, objname1)
		{

			char s[128];
			sprintf(s, "%d out of range [%d,%d]", input, low, high);
			desc = string(s);
		}
		
		const char *name() const { return "OutofRangeException"; }
		
	};
#define OutofRangeException(low, high, input, objname) \
 _OutofRangeException(low, high, input,  __FILE__, __LINE__, objname)

}
    
#endif
