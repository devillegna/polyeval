#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<stdint.h>
#include<assert.h>

#include "btfy.h"
#include "gf264.h"
#include "dencoder.h"

static void polymul_ref( uint64_t * c , const uint64_t * a , const uint64_t * b , unsigned len )
{
    for(unsigned i=0;i<len*2;i++) c[i] = 0;
    for(unsigned i=0;i<len;i++){
        for(unsigned j=0;j<len;j++){
            uint64_t ai = a[i];
            uint64_t bj = b[j];
            uint64_t c0=0;
            uint64_t c1=0;
            if( ai&1 ) { c0 ^= bj; }
            for(int k=1;k<64;k++){
                if ((ai>>k)&1) {
                    c0 ^= (bj << k);
                    c1 ^= (bj >> (64-k));
                }
            }
            c[i+j]   ^= c0;
            c[i+j+1] ^= c1;
        }
    }
}

static inline
unsigned log_floor(unsigned a)
{
    for(unsigned i=0;i<32;i++){
        if(0==a) return i-1;
        a &=  ~(1<<i);
    }
    return 32;
}

static int polymul_fafft( uint64_t *c, const uint64_t * a, const uint64_t * b , unsigned n_u64 )
{
    if (0==n_u64) return 0;
    uint64_t *a1 = malloc(n_u64*8*2*2);
    if (NULL==a1) return -1;
    uint64_t *b1 = &a1[2*n_u64];

    unsigned loglen = log_floor(n_u64);
    // input transform
    bc_encode_64(a1,a,n_u64);
    btfy_64(a1,loglen+1,1ULL<<(32+loglen+1));
    bc_encode_64(b1,b,n_u64);
    btfy_64(b1,loglen+1,1ULL<<(32+loglen+1));

    for(unsigned i=0;i<2*n_u64;i++) { a1[i] = gf264_mul(a1[i],b1[i]); }

    //int r = fafft_output_transform( c , a1 , 2*n_u64 );
    ibtfy_64(a1,loglen+1,1ULL<<(32+loglen+1));
    ibc_decode_64(c,a1,n_u64*2);

    free(a1);
    return 0;
}


#define MAXLEN   (128)


int main(void)
{
    uint64_t a[MAXLEN*2] = {0};
    uint64_t b[MAXLEN*2] = {0};
    uint64_t c0[MAXLEN*2];
    uint64_t c1[MAXLEN*2];
    int fail = 0;
    for(int j = 0; j < 1000; j++)
    {
        if (0==j) {
            a[0] = 0x41a7;
            b[0] = 0x60b7acd9;
        } else if (1==j) {
            a[1] = 0x1ffffffff;
            b[1] = 0x1ffffffff;
        } else if (2==j) {
            a[MAXLEN-1] = 0x1ffffffff;
            b[MAXLEN-1] = 0x1ffffffff;
        } else {
            int qq = (j < MAXLEN )? j : MAXLEN;
            for(int i = 0; i < qq; i++) {
                a[i] = ((uint64_t)rand())^(((uint64_t)rand())<<32);
                b[i] = ((uint64_t)rand())^(((uint64_t)rand())<<32);
            }
        }
        //for(int i=0;i<LEN2048*2;i++) {c0[i]=0;}
        //polymul(a,b,c0,LEN2048);
        polymul_ref(c0,a,b,MAXLEN);
        polymul_fafft(c1,a,b,MAXLEN);

        for(int i = 0; i < MAXLEN*2; i++)
        {
            if(c0[i]!=c1[i])
            {
                printf("test FAIL [%d,%d]: ",j,i);
                printf("%llx vs %llx \n", c0[i] , c1[i]);
                if (0==i) {
                    printf("a0: %llx , b0: %llx\n", a[0], b[0]);
                }
                fail = 1;
            }
            if( fail ) { printf("\n"); break; }
        }
        if( fail ) break;
    }
    if(!fail) {
        printf("test fafft_polymul( [%d]xu64 ) pass [%d].\n",MAXLEN, 1000);
    }
    return 0;
}
