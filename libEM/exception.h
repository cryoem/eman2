#ifndef eman__exception_h__
#define eman__exception_h__ 1

#include <string>

using std::string;

/** Exception class
 *
 * Exception class is a subclass of std::exception; All EMAN2
 * exception classes are subclass of Exception class.
 *
 * A XYZ exception class is defined in the following way:
 *   1) The class is named _XYZException.
 *   2) The class has a function to return its name "XYZException".
 *   3) A macro called "XYZException" is defined to simplify the usage
 *      of _XYZException class. 
 */
namespace EMAN {

	/** Exception class is the parent class of all EMAN2 exceptions.
	 */
    class Exception : public std::exception {
    public:
		/** Contructor.
		 * @params file The name of the file where an exception is
		 * thrown.
		 * @params line The line number in the file where the
		 * exception is thrown.
		 * @params desc_str Exception description string.
		 * @params objname_str Exception involved object name.
		 */ 
		Exception(const string& file = "", int line = 0,
				  const string& desc_str = "", const string& objname_str = "")
			: filename(file), linenum(line), desc(desc_str), objname(objname_str) {}
		
		virtual ~Exception() throw() {}

		/** The exception information.
		 * @return The exception information including exception
		 * location (filename, line number,function name) and
		 * description.
		 */
		virtual const char *what() const throw();

		/** The name of this exception class.
		 * @return The name of this exception class.
		 */
		virtual const char *name() const { return "Exception"; }
	protected:
		string filename;
		int linenum;
		string desc;
		string objname;
    };


	/** Used when an object type, like an EMObject type, doesn't
	 * exist.
	 * Parameters:
	 *   1. objname  The name of the not-existing object.
	 *   2. desc The description of the situation.
	 */
	class _NotExistingObjectException : public Exception
	{
	public:
		_NotExistingObjectException(const string& objname_str,
									const string& file = "unknown",
									int line = 0,
									const string& desc_str = "none")
			: Exception(file, line, desc_str, objname_str) {}
		
		const char *name() const { return "NotExistingObjectException"; }
		
	};
	
#define NotExistingObjectException(objname, desc)  \
 _NotExistingObjectException(objname, __FILE__, __LINE__, desc)

	/** Used when an image is not in the expected format.
	 */
	class _ImageFormatException : public Exception
	{
	public:
		_ImageFormatException(const string& desc_str,
							  const string& file = "unknown",
							  int line = 0)
			: Exception(file, line, desc_str) {}
		
		const char *name() const { return "ImageFormatException"; }
		
	};

#define ImageFormatException(desc) _ImageFormatException(desc, __FILE__, __LINE__)

	/** Used when an image is not in the expected dimension. For
	 * example, a 2D image is given when a 3D image is expected.
	 *  Parameters:
	 *    1. desc  The description of the situation.
	 */
	class _ImageDimensionException : public Exception
	{
	public:
		_ImageDimensionException(const string& desc_str, const string& file = "unknown",
								 int line = 0)
			: Exception(file, line, desc_str) {}
		
		const char *name() const { return "ImageDimensionException"; }
		
	};

#define ImageDimensionException(desc) _ImageDimensionException(desc, __FILE__, __LINE__)

	/** Used when a file access error occurs. For example, when you
	 * try to open a non-existing file or directory.
	 * Parameters:
	 *   1. filename  The name of the file with access error.
	 */
	class _FileAccessException : public Exception
	{
	public:
		_FileAccessException(const string& filename_str, const string& file = "unknown",
							 int line = 0, const string& desc_str = "")
			: Exception(file, line, desc_str, filename_str)
		{
			desc = "cannot access file '" + filename_str + "'";
		}
		
		const char *name() const { return "FileAccessException"; }
		
	};
#define FileAccessException(filename) _FileAccessException(filename, __FILE__, __LINE__)
	
	/** Used when an error occurs at image reading time.
	 * Parameters:
	 *   1. imagename The name of image with writing error.
	 *   2. desc The description of the situation.
	 */
	class _ImageReadException : public Exception
	{
	public:
		_ImageReadException(const string& imagename, const string& file = "unknown",
							int line = 0, const string& desc_str = "")
			: Exception(file, line, desc_str, imagename) {}
		
		const char *name() const { return "ImageReadException"; }
		
	};
#define ImageReadException(filename, desc) \
 _ImageReadException(filename, __FILE__, __LINE__, desc)
	
	/** Used when an error occurs at image writing time.
	 * Parameters:
	 *   1. imagename The name of image with writing error.
	 *   2. desc The description of the situation.
	 */
	class _ImageWriteException : public Exception
	{
	public:
		_ImageWriteException(const string& imagename, const string& file = "unknown",
							 int line = 0, const string& desc_str = "")
			: Exception(file, line, desc_str, imagename) {}
		
		const char *name() const { return "ImageWriteException"; }
		
	};
#define ImageWriteException(imagename, desc) \
 _ImageWriteException(imagename, __FILE__, __LINE__, desc)

	/** Used when a NULL is given to a pointer that should not be NULL.
	 * Parameter:
	 *   1. desc: The description string.
	 */
	class _NullPointerException : public Exception
	{
	public:
		_NullPointerException(const string& file = "unknown",
							  int line = 0, const string& desc_str = "")
			: Exception(file, line, desc_str) {}
		
		const char *name() const { return "NullPointerException"; }
		
	};
#define NullPointerException(desc) _NullPointerException(__FILE__, __LINE__, desc)
	
	/** Used when a type cast error occurs. For example, when casting
	 * an EMData* type to a float type.
	 *
	 * Parameters:
	 *  1. desc  Description of the situation.
	 *  2. type  The name of type causing trouble.
	 */
	class _TypeException : public Exception
	{
	public:
		_TypeException(const string & desc_str, const string & type,
					   const string & file = "unknown", int line = 0)
			: Exception(file, line, desc_str, type) {}
		
		const char *name() const { return "TypeException"; }
		
	};
#define TypeException(desc, type) _TypeException(desc, type, __FILE__, __LINE__)

	/** Used when an invalid integer value is given.
	 * Parameters:
	 *  1. val  The invalid integer value.
	 *  2. desc Description of the situation.
	 */
	class _InvalidValueException : public Exception
	{
	public:
		_InvalidValueException(int val, const string& file = "unknown",
							   int line = 0, const string& desc_str = "")
			: Exception(file, line, desc_str)
		{
			char s[32];
			sprintf(s, "%d", val);
			objname = string(s);
		}		
		const char *name() const { return "InvalidValueException"; }
		
	};
#define InvalidValueException(val, desc) \
 _InvalidValueException(val, __FILE__, __LINE__, desc)


	/** Used when the given value is out of range.
	 * parameters:
	 *  1. low    The lower limit of the valid range.
	 *  2. high   The uppper limit of the valid range.
	 *  3. input  The given, out-of-range value.
	 *  4. objname The name of the variable holding the value.
	 */
	class _OutofRangeException : public Exception
	{
	public:
		_OutofRangeException(int low, int high, int input,
							 const string& file = "unknown", 
							 int line = 0, const string & desc_str = "",
							 const string& objname_str = "")
			: Exception(file, line, desc_str, objname_str)
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
