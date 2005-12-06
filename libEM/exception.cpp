#include "util.h"


using namespace EMAN;

const char* E2Exception::what() const throw()
{
	string err1 = "";
	if (objname != "") {
		err1 = "error with '" + objname + "': ";
	}

	string msg = string(name()) + " at " + filename + ":" + Util::int2str(linenum);
	msg += ": " + err1 + "'" + desc + "' caught\n";
	return msg.c_str();
}
