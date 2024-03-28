
#ifndef _GF232_H_
#define _GF232_H_


#include <stdint.h>


/// GF(2^32) := x^32 + x^7 + x^3 + x^2 + 1   // 0x8d
static inline
uint32_t _gf232_mul( uint32_t a , uint32_t b )
{
  uint32_t r = a&(-(b&1));
  for(int i=1;i<32;i++) {
    uint32_t ar = ((a>>31)&1)*0x8d;
    a = (a<<1)^ar;  // a = a gf* 2;
    b >>=1;
    r ^= a&(-(b&1));
  }
  return r;
}

uint32_t gf232_mul( uint32_t a , uint32_t b );

uint32_t gf232_inv( uint32_t a );

#endif

