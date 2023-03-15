#ifndef _RND_BYTES_H_
#define _RND_BYTES_H_

#include "stdlib.h"

static inline
void randombytes( unsigned char * v , unsigned len ) { for(unsigned i=0;i<len;i++) v[i]=rand()&0xff; }

#endif
