
#ifndef _GF232_AESNI_H_
#define _GF232_AESNI_H_


#include <stdint.h>



#if 1

#include <immintrin.h>

// GF(2^32) := x^32 + x^7 + x^3 + x^2 + 1   // 0x8d

static const uint64_t _gf232_reducer[4] __attribute__((aligned(32)))  = {0xffffffffULL,0x8dULL,0x8dULL,0x8dULL};

static inline
__m128i _gf232_reduce_sse( __m128i x0 )
{
	__m128i reducer = _mm_load_si128( (__m128i const*)_gf232_reducer );
	//__m128i *reducer = (__m128i *)_gf2ext128_reducer;
	//__m128i r0 = _mm_clmulepi64_si128( x0 , reducer , 0x11 );


	__m128i r0 = (x0&reducer)^_mm_clmulepi64_si128( _mm_slli_si128(x0,4) , reducer , 0x11 );
	__m128i r1 = (r0&reducer)^_mm_clmulepi64_si128( _mm_slli_si128(r0,4) , reducer , 0x11 );
	return r1;
}


static inline
__m128i gf232_mulx1_sse( __m128i a0 , __m128i b0 )
{
	__m128i c0 = _mm_clmulepi64_si128( a0 , b0 , 0 );
	__m128i c3 = _gf232_reduce_sse( c0 );
	return c3;
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

	__m128i reducer = _mm_load_si128( (__m128i const*)(&_gf232_reducer[2]) );
	__m128i r0 = _mm_clmulepi64_si128( _mm_srli_si128(c0,4) , reducer , 0 );
	__m128i r1 = _mm_clmulepi64_si128( _mm_srli_si128(r0,4) , reducer , 0 );

	return _mm_extract_epi32( c0^r0^r1 , 0 );
}


#else

#include <immintrin.h>


/// X^64 + X^4 + X^3 + X + 1
/// 0x1b
static const uint64_t _gf2ext64_reducer[4] __attribute__((aligned(32)))  = {0x415A776C2D361B00ULL,0x1bULL,0x415A776C2D361B00ULL,0x1bULL};

static inline
__m128i _gf2ext64_reduce_sse( __m128i x0 )
{
	__m128i reducer = _mm_load_si128( (__m128i const*)_gf2ext64_reducer );
	//__m128i *reducer = (__m128i *)_gf2ext128_reducer;
	__m128i r0 = _mm_clmulepi64_si128( x0 , reducer , 0x11 );
	__m128i r1 = _mm_shuffle_epi8( reducer , r0 );
	__m128i r2 = x0 ^ r0;
	r1 ^= _mm_slli_si128( r2 , 8 );
	return _mm_srli_si128( r1 , 8 );
}


static inline
__m128i _gf2ext64_reduce_x2_sse( __m128i x0 , __m128i y0 )
{
	__m128i reducer = _mm_load_si128( (__m128i const*)_gf2ext64_reducer );
	//__m128i *reducer = (__m128i *)_gf2ext128_reducer;
	__m128i r0 = _mm_clmulepi64_si128( x0 , reducer , 0x11 );
	__m128i s0 = _mm_clmulepi64_si128( y0 , reducer , 0x11 );
	__m128i r2 = x0 ^ r0;
	__m128i s2 = y0 ^ s0;
	__m128i pr = _mm_unpacklo_epi64( r2 , s2 );

	__m128i rr = _mm_unpackhi_epi64( r0 , s0 );
	__m128i rr2 = _mm_shuffle_epi8( reducer , rr );

	return pr^rr2;
}

static inline
__m256i _gf2ext64_reduce_x4_avx2( __m128i w0 , __m128i x0 , __m128i y0 , __m128i z0 )
{
	__m256i reducer2 = _mm256_load_si256( (__m256i const*)_gf2ext64_reducer );
	//__m128i *reducer = (__m128i *)_gf2ext128_reducer;
	__m128i reducer = _mm256_castsi256_si128( reducer2 );
	__m128i r0 = _mm_clmulepi64_si128( w0 , reducer , 0x11 );
	__m128i s0 = _mm_clmulepi64_si128( x0 , reducer , 0x11 );
	__m128i t0 = _mm_clmulepi64_si128( y0 , reducer , 0x11 );
	__m128i u0 = _mm_clmulepi64_si128( z0 , reducer , 0x11 );
	__m128i r2 = w0 ^ r0;
	__m128i s2 = x0 ^ s0;
	__m128i t2 = y0 ^ t0;
	__m128i u2 = z0 ^ u0;
	__m256i pr1 = _mm256_castsi128_si256( _mm_unpacklo_epi64( r2 , s2 ) );
	__m256i pr2 = _mm256_inserti128_si256( pr1 , _mm_unpacklo_epi64( t2 , u2 ) , 1 );

	__m256i rr1 = _mm256_castsi128_si256( _mm_unpackhi_epi64( r0 , s0 ) );
	__m256i rr2 = _mm256_inserti128_si256( rr1 , _mm_unpackhi_epi64( t0 , u0 ) , 1 );
	return pr2 ^ _mm256_shuffle_epi8( reducer2 , rr2 );
}



static inline
__m128i _gf2ext64_mul_sse( __m128i a0 , __m128i b0 )
{
	__m128i c0 = _mm_clmulepi64_si128( a0 , b0 , 0 );
	__m128i c3 = _gf2ext64_reduce_sse( c0 );
	return c3;
}

static inline
__m128i _gf2ext64_mul_hi_sse( __m128i a0 , __m128i b0 )
{
	__m128i c0 = _mm_clmulepi64_si128( a0 , b0 , 1 );
	__m128i c3 = _gf2ext64_reduce_sse( c0 );
	return c3;
}

static inline
__m128i _gf2ext64_mul_2x1_sse( __m128i a0a1 , __m128i b0 )
{
	__m128i c0 = _mm_clmulepi64_si128( a0a1 , b0 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a0a1 , b0 , 1 );
	__m128i c3 = _gf2ext64_reduce_x2_sse( c0 , c1 );
	return c3;
}

static inline
__m128i _gf2ext64_mul_2x2_sse( __m128i a0a1 , __m128i b0b1 )
{
	__m128i c0 = _mm_clmulepi64_si128( a0a1 , b0b1 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( a0a1 , b0b1 , 0x11 );
	__m128i c3 = _gf2ext64_reduce_x2_sse( c0 , c1 );
	return c3;
}

static inline
__m256i _gf2ext64_mul_4x1_avx2( __m256i a , __m128i b0 )
{
	__m128i al = _mm256_castsi256_si128( a );
	__m128i c0 = _mm_clmulepi64_si128( al , b0 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( al , b0 , 1 );
	__m128i ah = _mm256_extracti128_si256( a , 1 );
	__m128i c2 = _mm_clmulepi64_si128( ah , b0 , 0 );
	__m128i c3 = _mm_clmulepi64_si128( ah , b0 , 1 );

	return _gf2ext64_reduce_x4_avx2( c0 , c1 , c2 , c3 );
}

static inline
__m256i _gf2ext64_mul_2x1_x2_avx2( __m256i a , __m128i b01 )
{
	__m128i al = _mm256_castsi256_si128( a );
	__m128i c0 = _mm_clmulepi64_si128( al , b01 , 0 );
	__m128i c1 = _mm_clmulepi64_si128( al , b01 , 1 );
	__m128i ah = _mm256_extracti128_si256( a , 1 );
	__m128i c2 = _mm_clmulepi64_si128( ah , b01 , 0x10 );
	__m128i c3 = _mm_clmulepi64_si128( ah , b01 , 0x11 );

	return _gf2ext64_reduce_x4_avx2( c0 , c1 , c2 , c3 );
}

static inline
__m256i _gf2ext64_mul_4x4_avx2( __m256i a , __m256i b )
{
	__m128i al = _mm256_castsi256_si128( a );
	__m128i bl = _mm256_castsi256_si128( b );
	__m128i c0 = _mm_clmulepi64_si128( al , bl , 0 );
	__m128i c1 = _mm_clmulepi64_si128( al , bl , 0x11 );
	__m128i ah = _mm256_extracti128_si256( a , 1 );
	__m128i bh = _mm256_extracti128_si256( b , 1 );
	__m128i c2 = _mm_clmulepi64_si128( ah , bh , 0 );
	__m128i c3 = _mm_clmulepi64_si128( ah , bh , 0x11 );

	return _gf2ext64_reduce_x4_avx2( c0 , c1 , c2 , c3 );
}



static inline
uint64_t gf2ext64_mul_u64( uint64_t a , uint64_t b )
{
	uint64_t tmp[4] __attribute__((aligned(32)));
	tmp[0] = a;
	tmp[1] = 0;
	__m128i a0 = _mm_load_si128( (__m128i const *)tmp );
	tmp[0] = b;
	__m128i b0 = _mm_load_si128( (__m128i const *)tmp );
	__m128i c0 = _gf2ext64_mul_sse( a0 , b0 );

	_mm_store_si128((__m128i*) tmp , c0 );
	return tmp[0];
}



static inline
void gf2ext64_mul_sse( uint8_t * c , const uint8_t * a , const uint8_t * b )
{
	uint64_t tmp[4] __attribute__((aligned(32)));
	for(unsigned i=0;i<8;i++) ((uint8_t*)tmp)[i] = a[i];
	tmp[1] = 0;
	__m128i a0 = _mm_load_si128( (__m128i const *)tmp );
	for(unsigned i=0;i<8;i++) ((uint8_t*)tmp)[i] = b[i];
	__m128i b0 = _mm_load_si128( (__m128i const *)tmp );
	__m128i c0 = _gf2ext64_mul_sse( a0 , b0 );

	_mm_store_si128((__m128i*) tmp , c0 );
	for(unsigned i=0;i<8;i++) c[i] = ((uint8_t*)tmp)[i];
#if 0
	__m128i a0 = _mm_load_si128( (__m128i const *)a );
	__m128i b0 = _mm_load_si128( (__m128i const *)b );
	__m128i c0 = _gf2ext64_mul_sse( a0 , b0 );

	_mm_store_si128((__m128i*) c , c0 );
#endif
}


static inline
void gf2ext64_mul_2x2_sse( uint8_t * c , const uint8_t * a , const uint8_t * b )
{
	__m128i a0a1 = _mm_load_si128( (__m128i const *)a );
	__m128i b0b1 = _mm_load_si128( (__m128i const *)b );
	__m128i c0c1 = _gf2ext64_mul_2x2_sse( a0a1 , b0b1 );

	_mm_store_si128((__m128i*) c , c0c1 );
}

static inline
void gf2ext64_mul_4x4_avx2( uint8_t * c , const uint8_t * a , const uint8_t * b )
{
	__m256i as = _mm256_load_si256( (__m256i const *)a );
	__m256i bs = _mm256_load_si256( (__m256i const *)b );
	__m256i cs = _gf2ext64_mul_4x4_avx2( as , bs );

	_mm256_store_si256((__m256i*) c , cs );
}



#endif

#endif
