
#ifndef _CANTOR_TO_GF232_H_
#define _CANTOR_TO_GF232_H_


#include <stdint.h>

// 1024*32
#define SIZE_TBL_CANTOR2X  (1<<15)

#ifdef  __cplusplus
extern  "C" {
#endif

extern const uint32_t cantor_to_gf232_2x[SIZE_TBL_CANTOR2X];

extern const uint32_t cantor_basis32[32];

static inline
uint32_t cantor_to_gf232( uint32_t cantor_idx )
{
  uint64_t r = 0;
  while( cantor_idx ) {
//#if 4==sizeof(unsigned)
//    r ^= cantor_basis32[ __builtin_ctz(cantor_idx) ];
//#else
    r ^= cantor_basis32[ __builtin_ctzl(cantor_idx) ];
//#endif
    cantor_idx &= (cantor_idx-1);
  }
  return r;
}

uint32_t index_to_gf232( uint32_t fft_position_index );


#ifdef  __cplusplus
}
#endif


#endif
