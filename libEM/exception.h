#ifndef __exception_h__
#define __exception_h__

#include <string>

using std::string;

namespace EMAN {
    
    class Exception : public std::exception {
    public:
	Exception(const string& file, int line_num = 0, const string& desc = "none");
	
	virtual ~Exception() throw() {}

	virtual void set_desc(const string& desc);
	
	virtual const char* get_file() const;
	virtual const char* get_desc() const;
	virtual int get_line_num() const;
	virtual const char* what() const throw();
	
    private:
	string filename;
	int line;
	string desc;
    };

    
#endif
