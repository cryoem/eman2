#include "exception.h"


using namespace EMAN;

const char* Exception::what() const throw()
{
	string err1 = "";
	if (objname != "") {
		err1 = "error with '" + objname + "': ";
	}
	string msg = filename + ":";
	char sline[32];
	sprintf(sline, "%d", line);
	msg += string(sline) + ": " + err1 + desc;
	return msg.c_str();
}
