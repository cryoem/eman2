#include "exception.h"
#include "log.h"

using namespace EMAN;

void Exception::dump() const
{
	string err1 = "";
	if (objname != "") {
		err1 = "error with '" + objname + "': ";
	}
	
	Log::logger()->error("%s:%d: %s%s\n",
						 filename.c_str(), line, err1.c_str(), desc.c_str());
	
}
