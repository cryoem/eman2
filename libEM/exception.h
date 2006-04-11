#ifndef eman__exception_h__
#define eman__exception_h__ 1

#include <string>

using std::string;

/** E2Exception class
 *
 * E2Exception class is a subclass of std::exception; All EMAN2
 * exception classes are subclass of E2Exception class.
 *
 * A XYZ Exception class is defined in the following way:
 *   0) It will extend E2Exception class.
 *   1) The class is named _XYZException.
 *   2) The class has a function to return its name "XYZException".
 *   3) A macro called "XYZException" is defined to simplify the usage
 *      of _XYZException class. So that filename, function name, and
 *      line number can be handled automatically.
 *
 * How to use XYZException:
 *
 *   1) To throw exception, use "throw XYZException(...)";
 *   2) To catch exception, use "catch (_XYZException & e) ...".
 */

namespace EMAN {

	/** E2Exception class is the parent class of all EMAN2 E2Exceptions.
	 */
    class E2Exception : public std::exception {
    public:
		/** Contructor.
		 * @param file The name of the file where an E2Exception is
		 * thrown.
		 * @param line The line number in the file where the
		 * E2Exception is thrown.
		 * @param desc_str E2Exception description string.
		 * @param objname_str E2Exception involved object name.
		 */ 
		E2Exception(const string& file = "", int line = 0,
				  const string& desc_str = "", const string& objname_str = "")
			: filename(file), linenum(line), desc(desc_str), objname(objname_str) {}
		
		virtual ~E2Exception() throw() {}

		/** The E2Exception information.
		 * @return The E2Exception information including E2Exception
		 * location (filename, line number,function name) and
		 * description.
		 */
		virtual const char *what() const throw();

		/** The name of this E2Exception class.
		 * @return The name of this E2Exception class.
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
	class _NotExistingObjectException : public E2Exception
	{
	public:
		_NotExistingObjectException(const string& objname_str,
									const string& file = "unknown",
									int line = 0,
									const string& desc_str = "none")
			: E2Exception(file, line, desc_str, objname_str) {}
		
		const char *name() const { return "NotExistingObjectException"; }
		
	};	
#define NotExistingObjectException(objname, desc)  \
 _NotExistingObjectException(objname, __FILE__, __LINE__, desc)


	/** Used when an image is not in the expected format.
	 */
	class _ImageFormatException : public E2Exception
	{
	public:
		_ImageFormatException(const string& desc_str,
							  const string& file = "unknown",
							  int line = 0)
			: E2Exception(file, line, desc_str) {}
		
		const char *name() const { return "ImageFormatException"; }
		
	};
#define ImageFormatException(desc) _ImageFormatException(desc, __FILE__, __LINE__)


	/** Used when an image is not in the expected dimension. For
	 * example, a 2D image is given when a 3D image is expected.
	 *  Parameters:
	 *    1. desc  The description of the situation.
	 */
	class _ImageDimensionException : public E2Exception
	{
	public:
		_ImageDimensionException(const string& desc_str, const string& file = "unknown",
								 int line = 0)
			: E2Exception(file, line, desc_str) {}
		
		const char *name() const { return "ImageDimensionException"; }
		
	};
#define ImageDimensionException(desc) _ImageDimensionException(desc, __FILE__, __LINE__)


	/** Used when a file access error occurs. For example, when you
	 * try to open a non-existing file or directory.
	 * Parameters:
	 *   1. filename  The name of the file with access error.
	 */
	class _FileAccessException : public E2Exception
	{
	public:
		_FileAccessException(const string& filename_str, const string& file = "unknown",
							 int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str, filename_str)
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
	class _ImageReadException : public E2Exception
	{
	public:
		_ImageReadException(const string& imagename, const string& file = "unknown",
							int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str, imagename) {}
		
		const char *name() const { return "ImageReadException"; }
		
	};
#define ImageReadException(filename, desc) \
 _ImageReadException(filename, __FILE__, __LINE__, desc)

	
	/** Used when an error occurs at image writing time.
	 * Parameters:
	 *   1. imagename The name of image with writing error.
	 *   2. desc The description of the situation.
	 */
	class _ImageWriteException : public E2Exception
	{
	public:
		_ImageWriteException(const string& imagename, const string& file = "unknown",
							 int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str, imagename) {}
		
		const char *name() const { return "ImageWriteException"; }
		
	};
#define ImageWriteException(imagename, desc) \
 _ImageWriteException(imagename, __FILE__, __LINE__, desc)


	/** Used when a NULL is given to a pointer that should not be NULL.
	 * Parameter:
	 *   1. desc: The description string.
	 */
	class _NullPointerException : public E2Exception
	{
	public:
		_NullPointerException(const string& file = "unknown",
							  int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str) {}
		
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
	class _TypeException : public E2Exception
	{
	public:
		_TypeException(const string & desc_str, const string & type,
					   const string & file = "unknown", int line = 0)
			: E2Exception(file, line, desc_str, type) {}
		
		const char *name() const { return "TypeException"; }
		
	};
#define TypeException(desc, type) _TypeException(desc, type, __FILE__, __LINE__)


	/** Used when an invalid integer value is given.
	 * Parameters:
	 *  1. val  The invalid integer value.
	 *  2. desc Description of the situation.
	 */
	class _InvalidValueException : public E2Exception
	{
	public:
		_InvalidValueException(int val, const string& file = "unknown",
							   int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str)
		{
			char s[32];
			sprintf(s, "%d", val);
			objname = string(s);
		}		
		_InvalidValueException(float val, const string& file = "unknown",
							   int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str)
		{
			char s[32];
			sprintf(s, "%g", val);
			objname = string(s);
		}		
		const char *name() const { return "InvalidValueException"; }
		
	};
#define InvalidValueException(val, desc) \
 _InvalidValueException(val, __FILE__, __LINE__, desc)


	/** Used when an invalid (format) string is given.
	 * Parameters:
	 *  1. str  The invalid integer value.
	 *  2. desc Description of the situation.
	 */
	class _InvalidStringException : public E2Exception
	{
	public:
		_InvalidStringException(const string& str, const string& file = "unknown",
							   int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str)
		{
			objname = str;
		}		
		const char *name() const { return "InvalidStringException"; }
		
	};
#define InvalidStringException(str, desc) \
 _InvalidStringException(str, __FILE__, __LINE__, desc)


	/** Used when the given value is out of range.
	 * parameters:
	 *  1. low    The lower limit of the valid range.
	 *  2. high   The uppper limit of the valid range.
	 *  3. input  The given, out-of-range value.
	 *  4. objname The name of the variable holding the value.
	 */	 
	class _OutofRangeException : public E2Exception
	{
	public:
		_OutofRangeException(int low, int high, int input,
							 const string& file = "unknown", 
							 int line = 0, const string & desc_str = "",
							 const string& objname_str = "")
			: E2Exception(file, line, desc_str, objname_str)
		{

			char s[128];
			sprintf(s, "%d out of range [%d,%d]", input, low, high);
			desc = string(s);
		}
		
		const char *name() const { return "OutofRangeException"; }
		
	};
#define OutofRangeException(low, high, input, objname) \
 _OutofRangeException(low, high, input,  __FILE__, __LINE__, objname)
 
 
 	class _InvalidCallException : public E2Exception
 	{
 	public:
 		_InvalidCallException(const string& file = "unknown",
							  int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str) {}
		
		const char *name() const { return "Invalid function call for this type"; }
		
	};
#define InvalidCallException(desc) _InvalidCallException(__FILE__, __LINE__, desc)
 		
	
	/***/
	class _InvalidParameterException : public E2Exception
	{
	public:
		_InvalidParameterException(const string& file = "unknown",
									int line = 0, const string& desc_str = "")
			: E2Exception(file, line, desc_str) {}
		
		const char *name() const { return "Invalid Parameter"; }
	};
#define InvalidParameterException(desc) _InvalidParameterException(__FILE__, __LINE__, desc)

}
    
    
#endif
