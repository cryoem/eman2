#include "exception.h"
#include "util.h"


using namespace EMAN;

const char* Exception::what() const throw()
{
	string err1 = "";
	if (objname != "") {
		err1 = "error with '" + objname + "': ";
	}

	string msg = "Exception at " + filename + ":" + Util::int2str(linenum);
	msg += ": " + err1 + "'" + desc + "' caught";
	return msg.c_str();
}
