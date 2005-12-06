#ifndef __em_assert_h_
#define __em_assert_h_

/** Define Assert() function that is effective only when -DDEBUG is used.
 */

#ifdef DEBUG
#include <cassert>
#define Assert(s) assert(s)
#else
#define Assert(s) ((void) (0))

#endif

#endif
