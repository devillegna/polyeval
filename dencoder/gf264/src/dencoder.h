
#ifndef _DENCODE_H_
#define _DENCODE_H_


#include <stdint.h>


#ifdef  __cplusplus
extern  "C" {
#endif


// sizeof output = 2 x sizeof input // assert(n_u64>=4)
void bc_encode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 );

// sizeof output = sizeof input // assert(n_u64>=8)
void ibc_decode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 );

// sizeof output = 2 x sizeof input // assert(n_u64>=4)
void encode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 );

// sizeof output = sizeof input // assert(n_u64>=8)
void decode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 );




#ifdef  __cplusplus
}
#endif


#endif
