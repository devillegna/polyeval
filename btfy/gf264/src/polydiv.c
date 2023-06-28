
#include "polydiv.h"


static inline void xor_down( uint64_t *poly_st_term , unsigned blk_size, int n_blk ) {
    uint64_t *dest = poly_st_term;
    for(int i=0;i<n_blk;i++) {
        dest -= blk_size;
        for(int j=(blk_size-1); j>=0;j--) dest[j] ^= poly_st_term[j];
    }
}



void polydiv( uint64_t *poly , int polylen , unsigned si )
{
    unsigned blk_size = 1<<si;
    unsigned len_blk = polylen>>si;
    for(unsigned i=2;i<len_blk;i+=2) {
        int log_n_blk = __builtin_ctz(i);
        xor_down( poly + i*blk_size , blk_size , (1<<log_n_blk)-1 );
    }
}



void ipolydiv( uint64_t *poly , int polylen , unsigned si )
{
    unsigned blk_size = 1<<si;
    unsigned len_blk = polylen>>si;
    len_blk = (len_blk&1)? len_blk : len_blk-1;
    for(unsigned i=len_blk-1;i>0;i-=2) {
        int log_n_blk = __builtin_ctz(i);
        xor_down( poly + i*blk_size , blk_size , (1<<log_n_blk)-1 );
    }
}



