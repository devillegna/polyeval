
#ifndef _GF232_NEON_H_
#define _GF232_NEON_H_


#include <stdint.h>

#include <arm_neon.h>


/// GF(2^32) := x^32 + x^7 + x^3 + x^2 + 1   // 0x8d

static inline
uint32x4_t _gf232_reducex4_neon( uint32x4_t ab0123_low , uint32x4_t ab0123_high )
{
  poly8x8_t _0x8d  = vdup_n_p8(0x8d);
  poly8x16_t h13 = vuzp2q_u8( ab0123_high , ab0123_high );
  uint32x4_t h13x0x8d = vreinterpretq_u32_p16( vmull_p8( vget_low_p8(h13) , _0x8d ) );
  ab0123_high ^= vshrq_n_u32( h13x0x8d , 24 );
  ab0123_low  ^= vshlq_n_u32( h13x0x8d , 8 );
  poly8x16_t h02 = vuzp1q_u8( ab0123_high , ab0123_high );
  ab0123_low ^= vreinterpretq_u32_p16( vmull_p8( vget_low_p8(h02) , _0x8d ));
  return ab0123_low;
}

static inline
uint32x4_t _gf232_mul4x4_neon( uint32x4_t a_x4 , uint32x4_t b_x4 )
{
  uint64x2_t a0_a1 = vmovl_u32( vget_low_u32(a_x4) );
  uint64x2_t a2_a3 = vmovl_high_u32( a_x4 );
  uint64x2_t b0_b1 = vmovl_u32( vget_low_u32(b_x4) );
  uint64x2_t b2_b3 = vmovl_high_u32( b_x4 );
  uint64x2_t ab0   = vreinterpretq_u64_p128( vmull_p64((poly64_t)vdupd_laneq_u64(a0_a1,0),(poly64_t)vdupd_laneq_u64(b0_b1,0)) );
  uint64x2_t ab1   = vreinterpretq_u64_p128( vmull_high_p64(a0_a1,b0_b1) );
  uint64x2_t ab2   = vreinterpretq_u64_p128( vmull_p64((poly64_t)vdupd_laneq_u64(a2_a3,0),(poly64_t)vdupd_laneq_u64(b2_b3,0)) );
  uint64x2_t ab3   = vreinterpretq_u64_p128( vmull_high_p64(a2_a3,b2_b3) );
  uint32x4_t ab01  = vzip1q_u64(ab0,ab1);
  uint32x4_t ab23  = vzip1q_u64(ab2,ab3);
  uint32x4_t ab0123_low  = vuzp1q_u32(ab01,ab23);
  uint32x4_t ab0123_high = vuzp2q_u32(ab01,ab23);
  return _gf232_reducex4_neon( ab0123_low , ab0123_high );
}

static inline
uint32x4_t _gf232_mul4x1_neon( uint32x4_t a_x4 , uint32_t b )
{
  uint64x2_t a0_a1 = vmovl_u32( vget_low_u32(a_x4) );
  uint64x2_t a2_a3 = vmovl_high_u32( a_x4 );
  uint64x2_t b0_b1 = vdupq_n_u64(b);
  uint64x2_t b2_b3 = b0_b1;
  uint64x2_t ab0   = vreinterpretq_u64_p128( vmull_p64((poly64_t)vdupd_laneq_u64(a0_a1,0),(poly64_t)vdupd_laneq_u64(b0_b1,0)) );
  uint64x2_t ab1   = vreinterpretq_u64_p128( vmull_high_p64(a0_a1,b0_b1) );
  uint64x2_t ab2   = vreinterpretq_u64_p128( vmull_p64((poly64_t)vdupd_laneq_u64(a2_a3,0),(poly64_t)vdupd_laneq_u64(b2_b3,0)) );
  uint64x2_t ab3   = vreinterpretq_u64_p128( vmull_high_p64(a2_a3,b2_b3) );
  uint32x4_t ab01  = vzip1q_u64(ab0,ab1);
  uint32x4_t ab23  = vzip1q_u64(ab2,ab3);
  uint32x4_t ab0123_low  = vuzp1q_u32(ab01,ab23);
  uint32x4_t ab0123_high = vuzp2q_u32(ab01,ab23);
  return _gf232_reducex4_neon( ab0123_low , ab0123_high );
}

//static inline uint64x2_t _gf232_mul2x2_neon( uint64x2_t a_x2 , uint64x2_t b_x2 ) { return a_x2^b_x2; }

static inline
uint32_t _gf232_mulx1_neon( uint32_t a , uint32_t b )
{
  uint32x4_t ab0 = vreinterpretq_u32_p128( vmull_p64((poly64_t)(a),(poly64_t)(b)) );
  // reduce
  uint32x4_t ab1r = vreinterpretq_u32_p128( vmull_p64((poly64_t)vgetq_lane_u32(ab0,1),(poly64_t)0x8d ));

  uint32x4_t r0 = ab0 ^ ab1r ^ vreinterpretq_u32_p16(vmull_p8( vreinterpret_p8_u32(vdup_laneq_u32(ab1r,1)),vreinterpret_p8_p64(vcreate_p64(0x8d))));

  return vgetq_lane_u32(r0,0);
}



#endif


