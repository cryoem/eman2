#ifndef eman__portable_fileio_h__
#define eman__portable_fileio_h__ 1

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1

#include <stdio.h>
#include <sys/types.h>

typedef off_t largefile_off_t;

inline int portable_fseek(FILE * fp, off_t offset, int whence)
{
    return fseeko(fp, offset, whence);
}

inline off_t portable_ftell(FILE * fp)
{
    return ftello(fp);
}

#endif
