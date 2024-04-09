
#ifndef _GF232_AESNI_H_
#define _GF232_AESNI_H_


#include <stdint.h>

#include <immintrin.h>
#include <wmmintrin.h>

// GF(2^32) := x^32 + x^7 + x^3 + x^2 + 1   // 0x8d


static inline
__m128i _gf232_reducex4_sse( __m128i cl , __m128i ch )
{
    __m128i rd0 = ch ^ _mm_slli_epi32(ch,2) ^ _mm_slli_epi32(ch,3) ^ _mm_slli_epi32(ch,7);
    __m128i rd1 = _mm_srli_epi32(ch, 30 ) ^ _mm_srli_epi32(ch,29) ^ _mm_srli_epi32(ch,25);  // 6-bit data
	__m128i rd2 = rd1 ^ _mm_slli_epi32(rd1,2) ^ _mm_slli_epi32(rd1,3) ^ _mm_slli_epi32(rd1,7);
	return cl ^ rd0 ^ rd2;
}


static inline
__m128i gf232_mul4x1_sse( __m128i a , uint32_t _b )
{
	__m128i b    = _mm_set1_epi64x(_b);
	__m128i zero = _mm_setzero_si128();
	__m128i a01 = _mm_unpacklo_epi32( a , zero );
	__m128i a23 = _mm_unpackhi_epi32( a , zero );

	__m128i c0 = _mm_clmulepi64_si128( a01 , b , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a01 , b , 1 );
	__m128i c2 = _mm_clmulepi64_si128( a23 , b , 0 );
	__m128i c3 = _mm_clmulepi64_si128( a23 , b , 1 );
	__m128i c02 = _mm_unpacklo_epi64( c0 , c2 );
	__m128i c13 = _mm_unpacklo_epi64( c1 , c3 );
	__m128i c0123_lo = _mm_blend_epi32( c02 , _mm_slli_si128(c13,4) , 0xa );
	__m128i c0123_hi = _mm_blend_epi32( _mm_srli_si128(c02,4) , c13 , 0xa );

    return _gf232_reducex4_sse( c0123_lo , c0123_hi );
}


static inline
__m256i _gf232_reducex8_avx2( __m256i cl , __m256i ch )
{
    __m256i rd0 = ch ^ _mm256_slli_epi32(ch,2) ^ _mm256_slli_epi32(ch,3) ^ _mm256_slli_epi32(ch,7);
    __m256i rd1 = _mm256_srli_epi32(ch, 30 ) ^ _mm256_srli_epi32(ch,29) ^ _mm256_srli_epi32(ch,25);
	__m256i rd2 = rd1 ^ _mm256_slli_epi32(rd1,2) ^ _mm256_slli_epi32(rd1,3) ^ _mm256_slli_epi32(rd1,7);
	return cl ^ rd0 ^ rd2;
}


static inline
__m256i gf232_mul8x1_avx2( __m256i a , uint32_t _b )
{
	__m128i b    = _mm_set1_epi64x(_b);
	__m128i zero = _mm_setzero_si128();
	__m128i a0123 = _mm256_extracti128_si256( a , 0 );
	__m128i a4567 = _mm256_extracti128_si256( a , 1 );
	__m128i a01 = _mm_unpacklo_epi32( a0123 , zero );
	__m128i a23 = _mm_unpackhi_epi32( a0123 , zero );
	__m128i a45 = _mm_unpacklo_epi32( a4567 , zero );
	__m128i a67 = _mm_unpackhi_epi32( a4567 , zero );

	__m128i c0 = _mm_clmulepi64_si128( a01 , b , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a01 , b , 1 );
	__m128i c2 = _mm_clmulepi64_si128( a23 , b , 0 );
	__m128i c3 = _mm_clmulepi64_si128( a23 , b , 1 );
	__m128i c4 = _mm_clmulepi64_si128( a45 , b , 0 );
	__m128i c5 = _mm_clmulepi64_si128( a45 , b , 1 );
	__m128i c6 = _mm_clmulepi64_si128( a67 , b , 0 );
	__m128i c7 = _mm_clmulepi64_si128( a67 , b , 1 );

	__m256i c0_4 = _mm256_set_m128i( c4 , c0 );
	__m256i c1_5 = _mm256_set_m128i( c5 , c1 );
	__m256i c2_6 = _mm256_set_m128i( c6 , c2 );
	__m256i c3_7 = _mm256_set_m128i( c7 , c3 );

	__m256i c02 = _mm256_unpacklo_epi64( c0_4 , c2_6 );
	__m256i c13 = _mm256_unpacklo_epi64( c1_5 , c3_7 );
	__m256i c0123_lo = _mm256_blend_epi32( c02 , _mm256_slli_epi64(c13,32) , 0xaa );
	__m256i c0123_hi = _mm256_blend_epi32( _mm256_srli_epi64(c02,32) , c13 , 0xaa );

    return _gf232_reducex8_avx2( c0123_lo , c0123_hi );
}

static inline
__m256i gf232_mul8x2_avx2( __m256i a , uint32_t _b0 , uint32_t _b1 )
{
	__m128i b0    = _mm_set1_epi64x(_b0);
	__m128i b1    = _mm_set1_epi64x(_b1);
	__m128i zero = _mm_setzero_si128();
	__m128i a0123 = _mm256_extracti128_si256( a , 0 );
	__m128i a4567 = _mm256_extracti128_si256( a , 1 );
	__m128i a01 = _mm_unpacklo_epi32( a0123 , zero );
	__m128i a23 = _mm_unpackhi_epi32( a0123 , zero );
	__m128i a45 = _mm_unpacklo_epi32( a4567 , zero );
	__m128i a67 = _mm_unpackhi_epi32( a4567 , zero );

	__m128i c0 = _mm_clmulepi64_si128( a01 , b0 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a01 , b0 , 1 );
	__m128i c2 = _mm_clmulepi64_si128( a23 , b0 , 0 );
	__m128i c3 = _mm_clmulepi64_si128( a23 , b0 , 1 );
	__m128i c4 = _mm_clmulepi64_si128( a45 , b1 , 0 );
	__m128i c5 = _mm_clmulepi64_si128( a45 , b1 , 1 );
	__m128i c6 = _mm_clmulepi64_si128( a67 , b1 , 0 );
	__m128i c7 = _mm_clmulepi64_si128( a67 , b1 , 1 );

	__m256i c0_4 = _mm256_set_m128i( c4 , c0 );
	__m256i c1_5 = _mm256_set_m128i( c5 , c1 );
	__m256i c2_6 = _mm256_set_m128i( c6 , c2 );
	__m256i c3_7 = _mm256_set_m128i( c7 , c3 );

	__m256i c02 = _mm256_unpacklo_epi64( c0_4 , c2_6 );
	__m256i c13 = _mm256_unpacklo_epi64( c1_5 , c3_7 );
	__m256i c0123_lo = _mm256_blend_epi32( c02 , _mm256_slli_epi64(c13,32) , 0xaa );
	__m256i c0123_hi = _mm256_blend_epi32( _mm256_srli_epi64(c02,32) , c13 , 0xaa );

    return _gf232_reducex8_avx2( c0123_lo , c0123_hi );
}

static inline
__m256i gf232_mul8x4_avx2( __m256i a , __m128i b )
{
	__m128i zero = _mm_setzero_si128();
	__m128i b01   = _mm_unpacklo_epi32( b , zero );
	__m128i b23   = _mm_unpackhi_epi32( b , zero );
	__m128i a0123 = _mm256_extracti128_si256( a , 0 );
	__m128i a4567 = _mm256_extracti128_si256( a , 1 );
	__m128i a01 = _mm_unpacklo_epi32( a0123 , zero );
	__m128i a23 = _mm_unpackhi_epi32( a0123 , zero );
	__m128i a45 = _mm_unpacklo_epi32( a4567 , zero );
	__m128i a67 = _mm_unpackhi_epi32( a4567 , zero );

	__m128i c0 = _mm_clmulepi64_si128( a01 , b01 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a01 , b01 , 1 );
	__m128i c2 = _mm_clmulepi64_si128( a23 , b01 , 0x10 );
	__m128i c3 = _mm_clmulepi64_si128( a23 , b01 , 0x11 );
	__m128i c4 = _mm_clmulepi64_si128( a45 , b23 , 0 );
	__m128i c5 = _mm_clmulepi64_si128( a45 , b23 , 1 );
	__m128i c6 = _mm_clmulepi64_si128( a67 , b23 , 0x10 );
	__m128i c7 = _mm_clmulepi64_si128( a67 , b23 , 0x11 );

	__m256i c0_4 = _mm256_set_m128i( c4 , c0 );
	__m256i c1_5 = _mm256_set_m128i( c5 , c1 );
	__m256i c2_6 = _mm256_set_m128i( c6 , c2 );
	__m256i c3_7 = _mm256_set_m128i( c7 , c3 );

	__m256i c02 = _mm256_unpacklo_epi64( c0_4 , c2_6 );
	__m256i c13 = _mm256_unpacklo_epi64( c1_5 , c3_7 );
	__m256i c0123_lo = _mm256_blend_epi32( c02 , _mm256_slli_epi64(c13,32) , 0xaa );
	__m256i c0123_hi = _mm256_blend_epi32( _mm256_srli_epi64(c02,32) , c13 , 0xaa );

    return _gf232_reducex8_avx2( c0123_lo , c0123_hi );
}

static inline
__m256i gf232_mul8x8_avx2( __m256i a , __m256i b )
{
	__m128i zero = _mm_setzero_si128();
	__m128i b0123 = _mm256_extracti128_si256( b , 0 );
	__m128i b4567 = _mm256_extracti128_si256( b , 1 );
	__m128i b01   = _mm_unpacklo_epi32( b0123 , zero );
	__m128i b23   = _mm_unpackhi_epi32( b0123 , zero );
	__m128i b45   = _mm_unpacklo_epi32( b4567 , zero );
	__m128i b67   = _mm_unpackhi_epi32( b4567 , zero );
	__m128i a0123 = _mm256_extracti128_si256( a , 0 );
	__m128i a4567 = _mm256_extracti128_si256( a , 1 );
	__m128i a01 = _mm_unpacklo_epi32( a0123 , zero );
	__m128i a23 = _mm_unpackhi_epi32( a0123 , zero );
	__m128i a45 = _mm_unpacklo_epi32( a4567 , zero );
	__m128i a67 = _mm_unpackhi_epi32( a4567 , zero );

	__m128i c0 = _mm_clmulepi64_si128( a01 , b01 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a01 , b01 , 0x11 );
	__m128i c2 = _mm_clmulepi64_si128( a23 , b23 , 0 );
	__m128i c3 = _mm_clmulepi64_si128( a23 , b23 , 0x11 );
	__m128i c4 = _mm_clmulepi64_si128( a45 , b45 , 0 );
	__m128i c5 = _mm_clmulepi64_si128( a45 , b45 , 0x11 );
	__m128i c6 = _mm_clmulepi64_si128( a67 , b67 , 0 );
	__m128i c7 = _mm_clmulepi64_si128( a67 , b67 , 0x11 );

	__m256i c0_4 = _mm256_set_m128i( c4 , c0 );
	__m256i c1_5 = _mm256_set_m128i( c5 , c1 );
	__m256i c2_6 = _mm256_set_m128i( c6 , c2 );
	__m256i c3_7 = _mm256_set_m128i( c7 , c3 );

	__m256i c02 = _mm256_unpacklo_epi64( c0_4 , c2_6 );
	__m256i c13 = _mm256_unpacklo_epi64( c1_5 , c3_7 );
	__m256i c0123_lo = _mm256_blend_epi32( c02 , _mm256_slli_epi64(c13,32) , 0xaa );
	__m256i c0123_hi = _mm256_blend_epi32( _mm256_srli_epi64(c02,32) , c13 , 0xaa );

    return _gf232_reducex8_avx2( c0123_lo , c0123_hi );
}


static inline
uint32_t gf232_mul_u32_sse( uint32_t a , uint32_t b )
{
	uint32_t tmp[4] __attribute__((aligned(32))) = {0};
	tmp[0] = a;
	__m128i a0 = _mm_load_si128( (__m128i const *)tmp );
	tmp[0] = b;
	__m128i b0 = _mm_load_si128( (__m128i const *)tmp );
	__m128i c0 = _mm_clmulepi64_si128( a0 , b0 , 0 );

	__m128i reducer = _mm_set1_epi64x( 0x8d );
	__m128i r0 = _mm_clmulepi64_si128( _mm_srli_si128(c0,4) , reducer , 0 );
	__m128i r1 = _mm_clmulepi64_si128( _mm_srli_si128(r0,4) , reducer , 0 );

	return _mm_extract_epi32( c0^r0^r1 , 0 );
}



#endif
