#ifndef __portable_fileio_h__
#define __portable_fileio_h__

#include <stdio.h>
#include <sys/types.h>


/* Our very own off_t-like type, 64-bit if possible */
#if !defined(HAVE_LARGEFILE_SUPPORT)
typedef off_t largefile_off_t;
#elif SIZEOF_OFF_T >= 8
typedef off_t largefile_off_t;
#elif SIZEOF_FPOS_T >= 8
typedef fpos_t largefile_off_t;
#else
#error "Large file support, but neither off_t nor fpos_t is large enough."
#endif


/* a portable fseek() function
   return 0 on success, non-zero on failure (with errno set) */
inline int portable_fseek(FILE *fp, largefile_off_t offset, int whence)
{
#if !defined(HAVE_LARGEFILE_SUPPORT)
    return fseek(fp, offset, whence);
#elif defined(HAVE_FSEEKO) && SIZEOF_OFF_T >= 8
    return fseeko(fp, offset, whence);
#elif defined(HAVE_FSEEK64)
    return fseek64(fp, offset, whence);
#elif defined(__BEOS__)
    return _fseek(fp, offset, whence);
#elif SIZEOF_FPOS_T >= 8
    /* lacking a 64-bit capable fseek(), use a 64-bit capable fsetpos()
       and fgetpos() to implement fseek()*/
    fpos_t pos;
    switch (whence) {
    case SEEK_END:
#ifdef _WIN32
	fflush(fp);
	if (_lseeki64(fileno(fp), 0, 2) == -1)
	    return -1;
#else
	if (fseek(fp, 0, SEEK_END) != 0)
	    return -1;
#endif
	/* fall through */
    case SEEK_CUR:
	if (fgetpos(fp, &pos) != 0)
	    return -1;
	offset += pos;
	break;
        /* case SEEK_SET: break; */
    }
    return fsetpos(fp, &offset);
#else
#error "Large file support, but no way to fseek."
#endif
}

/* a portable ftell() function
   Return -1 on failure with errno set appropriately, current file
   position on success */
inline largefile_off_t portable_ftell(FILE* fp)
{
#if !defined(HAVE_LARGEFILE_SUPPORT)
    return ftell(fp);
#elif defined(HAVE_FTELLO) && SIZEOF_OFF_T >= 8
    return ftello(fp);
#elif defined(HAVE_FTELL64)
    return ftell64(fp);
#elif SIZEOF_FPOS_T >= 8
    fpos_t pos;
    if (fgetpos(fp, &pos) != 0)
	return -1;
    return pos;
#else
#error "Large file support, but no way to ftell."
#endif
}



#endif
