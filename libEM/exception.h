#ifndef eman__exception_h__
#define eman__exception_h__ 1

#include <string>

using std::string;

namespace EMAN {
    
    class Exception : public std::exception {
    public:
		Exception(const string& file = "unknown", int line_num = 0,
				  const string& desc1 = "none") : filename(file), line(line_num),
												 desc(desc1) {}
	
		virtual ~Exception() throw() {}

		virtual void set_desc(const string& desc1)
		{
			desc = desc1;
		}
	
		virtual const char* get_file() const
		{
			return filename.c_str();
		}

		virtual const char* get_desc() const
		{
			return desc.c_str();
		}
		
		virtual int get_line_num() const
		{
			return line;
		}
		
	
    private:
		string filename;
		int line;
		string desc;
    };
}
    
#endif
