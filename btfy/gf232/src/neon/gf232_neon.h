
#ifndef _GF232_NEON_H_
#define _GF232_NEON_H_


#include <stdint.h>

#include <arm_neon.h>


/// GF(2^32) := x^32 + x^7 + x^3 + x^2 + 1   // 0x8d

#if 0
static inline
uint64x2_t _gf232_mulx2_neon( uint64x2_t a , uint64x2_t b , uint64x2_t mask_0x1b)
{
  //uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);

  uint64x2_t x0 = vreinterpretq_u64_p128( vmull_p64(vget_low_u64(a),vget_low_u64(b)) );
  uint64x2_t y0 = vreinterpretq_u64_p128( vmull_high_p64(a,b) );
  // reduce
  uint64x2_t xr0 = vreinterpretq_u64_p128( vmull_high_p64(x0,mask_0x1b) );  // X,XXXXXXXX XXXXXXXX  <-- 3bits in the high 64bits
  uint64x2_t yr0 = vreinterpretq_u64_p128( vmull_high_p64(y0,mask_0x1b) );  // X,XXXXXXXX XXXXXXXX  <-- 3bits in the high 64bits

  uint64x2_t x1 = x0 ^ xr0;
  uint64x2_t y1 = y0 ^ yr0;

  uint64x2_t yrxr = vzip2q_u64(xr0,yr0);
  uint64x2_t yx = vzip1q_u64(x1,y1);

  return yx ^ vmulq_p8(yrxr,mask_0x1b);
}
#endif

static inline
uint32_t _gf232_mulx1_neon( uint32_t a , uint32_t b )
{
  poly64x1_t mask_0x8d = vcreate_p64(0x8d);

  uint64x2_t ab0 = vreinterpretq_u64_p128( vmull_p64(vcreate_p64(a),vcreate_p64(b)) );
  // reduce
  uint64x2_t ab1r = vreinterpretq_u64_p128( vmull_p64(vget_low_p64(vshrq_n_u64(ab0,32)),mask_0x8d ));

  uint64x2_t r0 = ab0 ^ ab1r ^ vmull_p8(vget_low_p8(vshrq_n_u64(ab1r,32)),mask_0x8d);

  return vdups_laneq_u32(r0,0);
}



#endif


