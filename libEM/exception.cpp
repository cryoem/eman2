#include "exception.h"
#include "util.h"


using namespace EMAN;

const char* Exception::what() const throw()
{
	string err1 = "";
	if (objname != "") {
		err1 = "error with '" + objname + "': ";
	}

	string msg = filename + ":" + Util::int2str(line) + ": " + err1 + desc;
	return msg.c_str();
}
