#ifndef _ALLOC_H
#define _ALLOC_H

#include <stdlib.h>
#include <stdio.h>

#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CC_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define CC_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define CC_OURABS(a) (((a) >= 0) ? (a) : -(a))

#define CCcheck_val_2(val,msg){                                             \
if((val)){                                                              \
fflush(stdout);                                                     \
fprintf(stderr, "%s at %s, line %d\n",(msg),__FILE__,__LINE__ );    \
goto CLEAN;                                                         \
}                                                                       \
}

#define CCcheck_val(val,msg){                                               \
if((val)){                                                              \
fflush(stdout);                                                     \
fprintf(stderr, "%s at %s, line %d\n",(msg),__FILE__,__LINE__ );    \
}                                                                       \
}


#define CCcheck_NULL(item,msg){                                           \
if(!(item)){                                                            \
fflush(stdout);                                                     \
fprintf(stderr, "%s at %s, line %d\n",(msg),__FILE__,__LINE__ );    \
}                                                                       \
}

#define CCcheck_NULL_2(item,msg){                                           \
if(!(item)){                                                            \
fflush(stdout);                                                     \
fprintf(stderr, "%s at %s, line %d\n",(msg),__FILE__,__LINE__ );    \
val = 1;                                                             \
goto CLEAN;                                                         \
}                                                                       \
}

#define CC_SAFE_MALLOC(nnum,type)   (type *) CCutil_allocrus (((size_t) (nnum)) * sizeof (type))           \

#define CC_FREE(object,type) {                                             \
CCutil_freerus ((void *) (object));                                    \
object = (type *) NULL;                                                \
}

#define CC_IFFREE(object,type) {                                           \
if ((object)) CC_FREE ((object),type);                                 \
}

void
*CCutil_allocrus (size_t size),
*CCutil_reallocrus (void *ptr, size_t size),
CCutil_freerus (void *p);

/****************************************************************************/
/*                                                                          */
/*                             copy.c                                       */
/*                                                                          */
/****************************************************************************/

void fill_int(int *dst,int n,int v);
void fill_dbl(double *dst,int n,double v);
void fill_float(float *dst,int n,float v);
void fill_char( char *dst, int n, char v );
void acopy_int(const int *src,int *dst,int n);
void acopy_dbl(const double *src,double *dst,int n);

#endif


