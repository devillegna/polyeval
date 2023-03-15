
#include "bc_64.h"

static inline
void bc_64_256(uint64_t *poly, unsigned n_256 ) {
    for(unsigned i=0;i<n_256;i++){
        uint64_t *p = poly+i*4;
        p[2] ^= p[3];
        p[1] ^= p[2];
    }
}


static inline
void ibc_64_256(uint64_t *poly, unsigned n_256 ) {
    for(unsigned i=0;i<n_256;i++){
        uint64_t *p = poly+i*4;
        p[1] ^= p[2];
        p[2] ^= p[3];
    }
}





//////////////////////////////////

static inline
void _div_2terms_poly( uint64_t * poly , unsigned deg_h, unsigned deg_l ) {
    for(unsigned i=(deg_h*2-1);i>=deg_h;i--) poly[i-deg_h+deg_l] ^= poly[i];
}

// s8 = x^4 - x
// div   s8^i , s8^(i-1) , ... , s8^1
static
void repr_s2_64( uint64_t * poly , unsigned n_64 )
{
  unsigned log_n = __builtin_ctz( n_64 );
  for( int i=log_n-2-1;i>=0;i--) {
    unsigned s_h = 1<<(i+2);
    unsigned s_l = 1<<i;
    unsigned len = 1<<(i+3);
    // div X^{4*s_ht} - X^{s_tt}
    for(unsigned j=0;j<n_64;j+=len) { _div_2terms_poly( poly+j , s_h , s_l ); }
  }
}


static inline
void _idiv_2terms_poly( uint64_t * poly , unsigned deg_h, unsigned deg_l ) {
    for(unsigned i=deg_h; i<deg_h*2;i++) poly[i-deg_h+deg_l] ^= poly[i];
}

// s8 = x^4 - x
// div   s8^i , s8^(i-1) , ... , s8^1
static
void irepr_s2_64( uint64_t * poly , unsigned n_64 )
{
  unsigned log_n = __builtin_ctz( n_64 );
  for( unsigned i=0;i<log_n-2;i++) {
    unsigned s_h = 1<<(i+2);
    unsigned s_l = 1<<i;
    unsigned len = 1<<(i+3);
    // div X^{4*s_ht} - X^{s_tt}
    for(unsigned j=0;j<n_64;j+=len) { _idiv_2terms_poly( poly+j , s_h , s_l ); }
  }
}


/////////////////////////////////////////

#include "bc_256.h"


void bc_64( uint64_t * poly , unsigned n_64 )
{
  if(2>=n_64) return;
  if(4>=n_64) { bc_64_256(poly,1); return; }
  repr_s2_64( poly , n_64 );
  bc_64_256( poly , n_64>>2 );
  bc_256( poly , n_64>>2 );
}

void ibc_64( uint64_t * poly , unsigned n_64 )
{
  if(2>=n_64) return;
  if(4>=n_64) { ibc_64_256(poly,1); return; }
  ibc_256( poly , n_64>>2 );
  ibc_64_256( poly , n_64>>2 );
  irepr_s2_64( poly , n_64 );
}

