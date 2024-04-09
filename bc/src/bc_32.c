
#include "bc_32.h"

#include "bc_64.h"


//////////////////////////////////

inline void bc_32_128( uint32_t * p )
{
    p[2] ^= p[3];
    p[1] ^= p[2];
}

inline void bc_32_256( uint32_t * p )
{
  // div x^4 - x
  p[4] ^= p[7];
  p[3] ^= p[6];
  p[2] ^= p[5];
  p[1] ^= p[4];
  // div x^2 - x
  p[6] ^= p[7];
  p[5] ^= p[6];
  p[2] ^= p[3];
  p[1] ^= p[2];
}

inline void bc_32_512( uint32_t * poly , unsigned n_512 )
{
  for(unsigned i=0;i<n_512;i++){
      uint32_t *p = poly+i*16;
      // repr s2 = x^4 - x
      // div x^8 - x^2
      p[9] ^= p[15];
      p[8] ^= p[14];
      p[7] ^= p[13];
      p[6] ^= p[12];
      p[5] ^= p[11];
      p[4] ^= p[10];
      p[3] ^= p[9];
      p[2] ^= p[8];
      // div x^4 - x
      p[4+8] ^= p[7+8];
      p[3+8] ^= p[6+8];
      p[2+8] ^= p[5+8];
      p[1+8] ^= p[4+8];
      p[4] ^= p[7];
      p[3] ^= p[6];
      p[2] ^= p[5];
      p[1] ^= p[4];
      // div (s2)^2 - (s2)  -> div x^8 - x^4
      p[11] ^= p[15];
      p[10] ^= p[14];
      p[9] ^= p[13];
      p[8] ^= p[12];
      p[7] ^= p[11];
      p[6] ^= p[10];
      p[5] ^= p[9];
      p[4] ^= p[8];
      // div x^2 - x
      p[6+8] ^= p[7+8];
      p[5+8] ^= p[6+8];
      p[2+8] ^= p[3+8];
      p[1+8] ^= p[2+8];
      p[6] ^= p[7];
      p[5] ^= p[6];
      p[2] ^= p[3];
      p[1] ^= p[2];
  }
}

static inline
void _div_2terms_poly( uint32_t * poly , unsigned deg_h, unsigned deg_l ) {
    for(unsigned i=(deg_h*2-1);i>=deg_h;i--) poly[i-deg_h+deg_l] ^= poly[i];
}

// s4 = x^16 - x
// div   s4^i , s4^(i-1) , ... , s4^1
static
void repr_s4_32( uint32_t * poly , unsigned n_32 )
{
  unsigned log_n = __builtin_ctz( n_32 );
  for( int i=log_n-4-1;i>=0;i--) {
    unsigned s_h = 1<<(i+4);
    unsigned s_l = 1<<i;
    unsigned len = 1<<(i+5);
    // div X^{16*s_ht} - X^{s_tt}
    for(unsigned j=0;j<n_32;j+=len) { _div_2terms_poly( poly+j , s_h , s_l ); }
  }
}

#include "bc_512.h"

void bc_32( uint32_t * poly , unsigned n_32 )
{
  if(2>=n_32) return;
  if(4>=n_32) { bc_32_128(poly); return; }
  if(8>=n_32) { bc_32_256(poly); return; }
  repr_s4_32( poly , n_32 );
  bc_32_512( poly , n_32>>4 );
  bc_512( poly , n_32>>4 );
}

////////////////////////////////////////////////////////////

inline void ibc_32_128( uint32_t * p )
{
    p[1] ^= p[2];
    p[2] ^= p[3];
}

inline void ibc_32_256( uint32_t * p )
{
  // div x^2 - x
  p[5] ^= p[6];
  p[6] ^= p[7];
  p[1] ^= p[2];
  p[2] ^= p[3];
  // div x^4 - x
  p[1] ^= p[4];
  p[2] ^= p[5];
  p[3] ^= p[6];
  p[4] ^= p[7];
}

inline void ibc_32_512( uint32_t * poly , unsigned n_512 )
{
  for(unsigned i=0;i<n_512;i++){
      uint32_t *p = poly+i*16;
      // idiv x^2 - x
      p[5+8] ^= p[6+8];
      p[6+8] ^= p[7+8];
      p[1+8] ^= p[2+8];
      p[2+8] ^= p[3+8];
      p[5] ^= p[6];
      p[6] ^= p[7];
      p[1] ^= p[2];
      p[2] ^= p[3];

      // idiv (s2)^2 - (s2)  -> div x^8 - x^4
      p[7] ^= p[11];
      p[6] ^= p[10];
      p[5] ^= p[9];
      p[4] ^= p[8];
      p[11] ^= p[15];
      p[10] ^= p[14];
      p[9] ^= p[13];
      p[8] ^= p[12];

      // irepr s2 = x^4 - x
      // idiv x^4 - x
      p[1+8] ^= p[4+8];
      p[2+8] ^= p[5+8];
      p[3+8] ^= p[6+8];
      p[4+8] ^= p[7+8];
      p[1] ^= p[4];
      p[2] ^= p[5];
      p[3] ^= p[6];
      p[4] ^= p[7];

      // idiv x^8 - x^2
      p[2] ^= p[8];
      p[3] ^= p[9];
      p[4] ^= p[10];
      p[5] ^= p[11];
      p[6] ^= p[12];
      p[7] ^= p[13];
      p[8] ^= p[14];
      p[9] ^= p[15];
  }
}

static inline
void _idiv_2terms_poly( uint32_t * poly , unsigned deg_h, unsigned deg_l ) {
    for(unsigned i=deg_h; i<deg_h*2;i++) poly[i-deg_h+deg_l] ^= poly[i];
}

// s4 = x^16 - x
// div   s4^i , s4^(i-1) , ... , s4^1
static
void irepr_s4_32( uint32_t * poly , unsigned n_32 )
{
  unsigned log_n = __builtin_ctz( n_32 );
  for(unsigned i=0;i<log_n-4;i++) {
    unsigned s_h = 1<<(i+4);
    unsigned s_l = 1<<i;
    unsigned len = 1<<(i+5);
    // div X^{16*s_ht} - X^{s_tt}
    for(unsigned j=0;j<n_32;j+=len) { _idiv_2terms_poly( poly+j , s_h , s_l ); }
  }
}

/////////////////////////////////////////


void ibc_32( uint32_t * poly , unsigned n_32 )
{
  if(2>=n_32) return;
  if(4>=n_32) { ibc_32_128(poly); return; }
  if(8>=n_32) { ibc_32_256(poly); return; }
  ibc_512( poly , n_32>>4 );
  ibc_32_512( poly , n_32>>4 );
  irepr_s4_32( poly , n_32 );
}

