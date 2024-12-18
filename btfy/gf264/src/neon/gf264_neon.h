
#ifndef _GF264_NEON_H_
#define _GF264_NEON_H_


#include <stdint.h>

#include <arm_neon.h>


/// X^64 + X^4 + X^3 + X + 1
/// 0x1b
//static const uint64_t _gf264_reducer[2] __attribute__((aligned(32)))  = {0x1bULL,0x1bULL};


static inline
uint64x2_t _gf264_mul_neon( uint64x2_t a , uint64x2_t b , uint64x2_t mask_0x1b)
{
  //uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);

  uint64x2_t x0 = vreinterpretq_u64_p128( vmull_p64((poly64_t)vdupd_laneq_u64(a,0),(poly64_t)vdupd_laneq_u64(b,0)) );
  uint64x2_t y0 = vreinterpretq_u64_p128( vmull_high_p64(a,b) );
  // reduce
  uint64x2_t xr0 = vreinterpretq_u64_p128( vmull_high_p64(x0,mask_0x1b) );  // X,XXXXXXXX XXXXXXXX  <-- 3bits in the high 64bits
  uint64x2_t yr0 = vreinterpretq_u64_p128( vmull_high_p64(y0,mask_0x1b) );  // X,XXXXXXXX XXXXXXXX  <-- 3bits in the high 64bits

  uint64x2_t x1 = x0 ^ xr0;
  uint64x2_t y1 = y0 ^ yr0;

  uint64x2_t yrxr = vzip2q_u64(xr0,yr0);
  uint64x2_t yx = vzip1q_u64(x1,y1);

  return yx ^ vreinterpretq_u64_p8(vmulq_p8(yrxr,mask_0x1b));
}

static inline
uint64_t _gf264_mulx1_neon( uint64_t a , uint64_t b )
{
  uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);

  uint64x2_t ab0 = vreinterpretq_u64_p128( vmull_p64((poly64_t)a,(poly64_t)b) );
  // reduce
  uint64x2_t ab1r = vreinterpretq_u64_p128( vmull_high_p64(ab0,mask_0x1b) );  // X,XXXXXXXX XXXXXXXX  <-- 3bits in the high 64bits

  uint64x2_t r0 = ab0 ^ ab1r ^ vreinterpretq_u64_p8(vmulq_p8(vtrn2q_u64(ab1r,ab1r),mask_0x1b));  // only low 64 bits is meaningful

  return vdupd_laneq_u64(r0,0);
}



#endif


