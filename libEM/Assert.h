#ifndef __em_assert_h_
#define __em_assert_h_

#ifdef DEBUG
#include <assert.h>
#define Assert(s) assert(s)
#else
#define Assert(s) ((void) (0))

#endif

#endif
