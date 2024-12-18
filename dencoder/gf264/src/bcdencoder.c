
#include "dencoder.h"

#include "bc_1.h"
// n_byte >= 32
//void bc_1( void * poly , unsigned n_byte );
//void ibc_1( void * poly , unsigned n_byte );

#include "string.h"


// sizeof output = 2 x sizeof input // assert(n_u64>=4)
void bc_encode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
{
    uint64_t temp[n_u64];
    memcpy( temp , fx , n_u64*8 );
    bc_1(temp, n_u64*8);
    encode_64(rfx,temp,n_u64);
}

// sizeof output = sizeof input // assert(n_u64>=8)
void ibc_decode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
{
    decode_64(rfx,fx,n_u64);
    ibc_1(rfx, n_u64*8);
}




