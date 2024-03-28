
#include "bc_32.h"


//////////////////////////////////

static inline
void _div_2terms_poly( uint32_t * poly , unsigned deg_h, unsigned deg_l ) {
    for(unsigned i=(deg_h*2-1);i>=deg_h;i--) poly[i-deg_h+deg_l] ^= poly[i];
}

// s1 = x^2 - x
// div   s1^i , s1^(i-1) , ... , s1^1
static
void repr_s1_32( uint32_t * poly , unsigned n_32 )
{
  unsigned log_n = __builtin_ctz( n_32 );
  for( int i=log_n-1-1;i>=0;i--) {
    unsigned s_h = 1<<(i+1);
    unsigned s_l = 1<<i;
    unsigned len = 1<<(i+2);
    // div X^{2*s_ht} - X^{s_tt}
    for(unsigned j=0;j<n_32;j+=len) { _div_2terms_poly( poly+j , s_h , s_l ); }
  }
}


static inline
void _idiv_2terms_poly( uint32_t * poly , unsigned deg_h, unsigned deg_l ) {
    for(unsigned i=deg_h; i<deg_h*2;i++) poly[i-deg_h+deg_l] ^= poly[i];
}

// s1 = x^2 - x
// div   s1^i , s1^(i-1) , ... , s1^1
static
void irepr_s1_32( uint32_t * poly , unsigned n_32 )
{
  unsigned log_n = __builtin_ctz( n_32 );
  for( unsigned i=0;i<log_n-1;i++) {
    unsigned s_h = 1<<(i+1);
    unsigned s_l = 1<<i;
    unsigned len = 1<<(i+2);
    // div X^{2*s_ht} - X^{s_tt}
    for(unsigned j=0;j<n_32;j+=len) { _idiv_2terms_poly( poly+j , s_h , s_l ); }
  }
}


/////////////////////////////////////////

#include "bc_64.h"


void bc_32( uint32_t * poly , unsigned n_32 )
{
  if(2>=n_32) return;
  repr_s1_32( poly , n_32 );
  bc_64( (uint64_t*)poly , n_32>>1 );  // XXX: supress warning. might have alignment issue.
}

void ibc_32( uint32_t * poly , unsigned n_32 )
{
  if(2>=n_32) return;
  ibc_64( (uint64_t*)poly , n_32>>1 );  // XXX: supress warning. might have alignment issue.
  irepr_s1_32( poly , n_32 );
}

