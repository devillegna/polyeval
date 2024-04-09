

#include <stdint.h>
#include <string.h>

#include "btfy32.h"
#include "cantor_to_gf232.h"

#include "gf232_aesni.h"

static inline unsigned min(unsigned a,unsigned b) { return (a<b)?a:b; }


static inline
void butterfly_32_x8( uint32_t * poly , unsigned unit ,  uint32_t a )
{
	unsigned unit_2= unit/2;  // ASSERT  (uint_2 & 7) == 0
	for(unsigned i=0;i<unit_2;i+=8) {
		__m256i p0 = _mm256_loadu_si256( (__m256i*)(poly+i) );
		__m256i p1 = _mm256_loadu_si256( (__m256i*)(poly+unit_2+i) );
		p0 ^= gf232_mul8x1_avx2( p1 , a );
		p1 ^= p0;
		_mm256_storeu_si256( (__m256i*)(poly+i) , p0 );
		_mm256_storeu_si256( (__m256i*)(poly+unit_2+i) , p1 );
	}
}



static inline
void btfy_s0_x1_32( uint32_t * poly , uint32_t a )
{
	poly[0] ^= gf232_mul_u32_sse( poly[1] , a );
	poly[1] ^= poly[0];
}

static inline
void btfy_s1s0_x1_32( uint32_t * poly , uint32_t idx_offset )
{
	uint32_t a_s1 = cantor_to_gf232( idx_offset>>1 );
	uint32_t a0_s0 = cantor_to_gf232( idx_offset );
	uint32_t a2_s0 = a0_s0 ^ cantor_basis32[1];

	uint32_t p0 = poly[0];
	uint32_t p1 = poly[1];
	uint32_t p2 = poly[2];
	uint32_t p3 = poly[3];

	p0 ^= gf232_mul_u32_sse( p2 , a_s1 );
	p1 ^= gf232_mul_u32_sse( p3 , a_s1 );
	p2 ^= p0;
	p3 ^= p1;
	p0 ^= gf232_mul_u32_sse( p1 , a0_s0 );
	p1 ^= p0;
	p2 ^= gf232_mul_u32_sse( p3 , a2_s0 );
	p3 ^= p2;

	poly[0] = p0;
	poly[1] = p1;
	poly[2] = p2;
	poly[3] = p3;
}

static inline
void btfy_s2s1s0_x1_32( uint32_t * poly , uint32_t idx_offset )
{
	uint32_t a_s2   = cantor_to_gf232( idx_offset>>2 );
	uint32_t a_s1_0 = cantor_to_gf232( idx_offset>>1 );
	uint32_t a_s1_1 = a_s1_0 ^ cantor_basis32[1];
	uint32_t a_s0_0 = cantor_to_gf232( idx_offset );
	uint32_t a_s0_1 = a_s0_0 ^ cantor_basis32[1];
	uint32_t a_s0_2 = a_s0_0 ^ cantor_basis32[2];
	uint32_t a_s0_3 = a_s0_1 ^ cantor_basis32[2];

	uint32_t p0 = poly[0];
	uint32_t p1 = poly[1];
	uint32_t p2 = poly[2];
	uint32_t p3 = poly[3];
	uint32_t p4 = poly[4];
	uint32_t p5 = poly[5];
	uint32_t p6 = poly[6];
	uint32_t p7 = poly[7];

	p0 ^= gf232_mul_u32_sse( p4 , a_s2 );
	p1 ^= gf232_mul_u32_sse( p5 , a_s2 );
	p2 ^= gf232_mul_u32_sse( p6 , a_s2 );
	p3 ^= gf232_mul_u32_sse( p7 , a_s2 );
	p4 ^= p0;
	p5 ^= p1;
	p6 ^= p2;
	p7 ^= p3;

	p0 ^= gf232_mul_u32_sse( p2 , a_s1_0 );
	p1 ^= gf232_mul_u32_sse( p3 , a_s1_0 );
	p2 ^= p0;
	p3 ^= p1;
	p4 ^= gf232_mul_u32_sse( p6 , a_s1_1 );
	p5 ^= gf232_mul_u32_sse( p7 , a_s1_1 );
	p6 ^= p4;
	p7 ^= p5;

	p0 ^= gf232_mul_u32_sse( p1 , a_s0_0 );
	p1 ^= p0;
	p2 ^= gf232_mul_u32_sse( p3 , a_s0_1 );
	p3 ^= p2;
	p4 ^= gf232_mul_u32_sse( p5 , a_s0_2 );
	p5 ^= p4;
	p6 ^= gf232_mul_u32_sse( p7 , a_s0_3 );
	p7 ^= p6;

	poly[0] = p0;
	poly[1] = p1;
	poly[2] = p2;
	poly[3] = p3;
	poly[4] = p4;
	poly[5] = p5;
	poly[6] = p6;
	poly[7] = p7;
}

///////////////////////////////////////////


static inline
void btfy_s3s2s1s0_x1_32( uint32_t * poly ,  uint32_t a_s3, uint32_t a_s2, uint32_t a_s1, uint32_t a_s0 , uint32_t idx_offset )
{
	uint32_t a_s3_0 = a_s3 ^ cantor_to_gf232_2x[ idx_offset>>3 ];
	uint32_t a_s2_0 = a_s2 ^ cantor_to_gf232_2x[ (idx_offset>>2) ];
	uint32_t a_s2_1 = a_s2 ^ cantor_to_gf232_2x[ (idx_offset>>2)+1 ];

	__m128i s1_a = _mm_set1_epi32( a_s1 ) ^ _mm_load_si128( (__m128i*)(cantor_to_gf232_2x + (idx_offset>>1)) );
	__m256i s0_a = _mm256_set1_epi32( a_s0 ) ^ _mm256_load_si256( (__m256i*)(cantor_to_gf232_2x + idx_offset) );


	__m256i p0 = _mm256_loadu_si256( (__m256i*)(poly) );
	__m256i p1 = _mm256_loadu_si256( (__m256i*)(poly+8) );
	p0 ^= gf232_mul8x1_avx2( p1 , a_s3_0 );
	p1 ^= p0;

	__m256i p2,p3;
	p2 = _mm256_permute2x128_si256( p0 , p1 , 0x20 );
	p3 = _mm256_permute2x128_si256( p0 , p1 , 0x31 );
	p2 ^= gf232_mul8x2_avx2( p3 , a_s2_0 , a_s2_1 );
	p3 ^= p2;

	p0 = _mm256_shuffle_pd( p2 , p3 , 0x0 );
	p1 = _mm256_shuffle_pd( p2 , p3 , 0xf );
	p0 ^= gf232_mul8x4_avx2( p1 , s1_a );
	p1 ^= p0;

	p2 = _mm256_blend_epi32( p0 , _mm256_slli_epi64(p1,32) , 0xaa );
	p3 = _mm256_blend_epi32( _mm256_srli_epi64(p0,32) , p1 , 0xaa );
	p2 ^= gf232_mul8x8_avx2( p3 , s0_a );
	p3 ^= p2;

	p0 = _mm256_blend_epi32( p2 , _mm256_slli_epi64(p3,32) , 0xaa );
	p1 = _mm256_blend_epi32( _mm256_srli_epi64(p2,32) , p3 , 0xaa );

	p2 = _mm256_shuffle_pd( p0 , p1 , 0x0 );
	p3 = _mm256_shuffle_pd( p0 , p1 , 0xf );

	p0 = _mm256_permute2x128_si256( p2 , p3 , 0x20 );
	p1 = _mm256_permute2x128_si256( p2 , p3 , 0x31 );

	_mm256_storeu_si256( (__m256i*)(poly) , p0 );
	_mm256_storeu_si256( (__m256i*)(poly+8) , p1 );

}


/////////////////////////////////////////////////////////


void btfy_32( uint32_t * poly , unsigned log_n , uint32_t idx_offset )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { btfy_s0_x1_32( poly , cantor_to_gf232(idx_offset) ); return; }
	if( 2 == log_n ) { btfy_s1s0_x1_32( poly , idx_offset ); return; }
	if( 3 == log_n ) { btfy_s2s1s0_x1_32( poly , idx_offset ); return; }

    unsigned log_unit=log_n;
	for(;log_unit>4; log_unit--) {
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));

		unsigned i = log_unit-1;
		uint32_t extra_a = cantor_to_gf232(idx_offset>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf232(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint32_t a = cantor_to_gf232_2x[k]^extra_a;
				butterfly_32_x8( poly + (j+k)*unit , unit , a );
			}
		}
	}
	if(4==log_unit){
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));
		uint32_t s3_idx = cantor_to_gf232(idx_offset>>3);
		uint32_t s2_idx = cantor_to_gf232(idx_offset>>2);
		uint32_t s1_idx = cantor_to_gf232(idx_offset>>1);
		uint32_t s0_idx = cantor_to_gf232(idx_offset);

		uint32_t last_s0_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X/8,num);  // use 8 element per iteration.

		for(unsigned j=0;j<num;j+=tab_size) {
			uint32_t s0_diff = (j<<4)^last_s0_idx;
			last_s0_idx = (j<<4);
			s3_idx ^= cantor_to_gf232(s0_diff>>4);
			s2_idx ^= cantor_to_gf232(s0_diff>>3);
			s1_idx ^= cantor_to_gf232(s0_diff>>2);
			s0_idx ^= cantor_to_gf232(s0_diff>>1);
			for(unsigned k=0;k<tab_size;k++) {
				btfy_s3s2s1s0_x1_32( poly+(j+k)*unit , s3_idx, s2_idx, s1_idx, s0_idx , k*8 );  // use 8 element per iteration
			}
		}
	}
}





/////////////////////////////////////////////////////////


static inline
void i_butterfly_32_x8( uint32_t * poly , unsigned unit ,  uint32_t a )
{
	unsigned unit_2= unit/2;  // ASSERT  (uint_2 & 7) == 0
	for(unsigned i=0;i<unit_2;i+=8) {
		__m256i p0 = _mm256_loadu_si256( (__m256i*)(poly+i) );
		__m256i p1 = _mm256_loadu_si256( (__m256i*)(poly+unit_2+i) );
		p1 ^= p0;
		p0 ^= gf232_mul8x1_avx2( p1 , a );
		_mm256_storeu_si256( (__m256i*)(poly+i) , p0 );
		_mm256_storeu_si256( (__m256i*)(poly+unit_2+i) , p1 );
	}
}

///////////////////////////////////////////


static inline
void i_btfy_s0_x1_32( uint32_t * poly , uint32_t a )
{
	poly[1] ^= poly[0];
	poly[0] ^= gf232_mul_u32_sse( poly[1] , a );
}

static inline
void i_btfy_s1s0_x1_32( uint32_t * poly , uint32_t idx_offset )
{
	uint32_t a_s1 = cantor_to_gf232( idx_offset>>1 );
	uint32_t a0_s0 = cantor_to_gf232( idx_offset );
	uint32_t a2_s0 = a0_s0 ^ cantor_basis32[1];

	uint32_t p0 = poly[0];
	uint32_t p1 = poly[1];
	uint32_t p2 = poly[2];
	uint32_t p3 = poly[3];

	p1 ^= p0;
	p0 ^= gf232_mul_u32_sse( p1 , a0_s0 );
	p3 ^= p2;
	p2 ^= gf232_mul_u32_sse( p3 , a2_s0 );

	p2 ^= p0;
	p3 ^= p1;
	p0 ^= gf232_mul_u32_sse( p2 , a_s1 );
	p1 ^= gf232_mul_u32_sse( p3 , a_s1 );

	poly[0] = p0;
	poly[1] = p1;
	poly[2] = p2;
	poly[3] = p3;
}

static inline
void i_btfy_s2s1s0_x1_32( uint32_t * poly , uint32_t idx_offset )
{
	uint32_t a_s2   = cantor_to_gf232( idx_offset>>2 );
	uint32_t a_s1_0 = cantor_to_gf232( idx_offset>>1 );
	uint32_t a_s1_1 = a_s1_0 ^ cantor_basis32[1];
	uint32_t a_s0_0 = cantor_to_gf232( idx_offset );
	uint32_t a_s0_1 = a_s0_0 ^ cantor_basis32[1];
	uint32_t a_s0_2 = a_s0_0 ^ cantor_basis32[2];
	uint32_t a_s0_3 = a_s0_1 ^ cantor_basis32[2];

	uint32_t p0 = poly[0];
	uint32_t p1 = poly[1];
	uint32_t p2 = poly[2];
	uint32_t p3 = poly[3];
	uint32_t p4 = poly[4];
	uint32_t p5 = poly[5];
	uint32_t p6 = poly[6];
	uint32_t p7 = poly[7];

	p1 ^= p0;
	p3 ^= p2;
	p5 ^= p4;
	p7 ^= p6;
	p0 ^= gf232_mul_u32_sse( p1 , a_s0_0 );
	p2 ^= gf232_mul_u32_sse( p3 , a_s0_1 );
	p4 ^= gf232_mul_u32_sse( p5 , a_s0_2 );
	p6 ^= gf232_mul_u32_sse( p7 , a_s0_3 );

	p2 ^= p0;
	p3 ^= p1;
	p6 ^= p4;
	p7 ^= p5;
	p0 ^= gf232_mul_u32_sse( p2 , a_s1_0 );
	p1 ^= gf232_mul_u32_sse( p3 , a_s1_0 );
	p4 ^= gf232_mul_u32_sse( p6 , a_s1_1 );
	p5 ^= gf232_mul_u32_sse( p7 , a_s1_1 );

	p4 ^= p0;
	p5 ^= p1;
	p6 ^= p2;
	p7 ^= p3;
	p0 ^= gf232_mul_u32_sse( p4 , a_s2 );
	p1 ^= gf232_mul_u32_sse( p5 , a_s2 );
	p2 ^= gf232_mul_u32_sse( p6 , a_s2 );
	p3 ^= gf232_mul_u32_sse( p7 , a_s2 );

	poly[0] = p0;
	poly[1] = p1;
	poly[2] = p2;
	poly[3] = p3;
	poly[4] = p4;
	poly[5] = p5;
	poly[6] = p6;
	poly[7] = p7;
}

///////////////////////////////////////////


static inline
void i_btfy_s3s2s1s0_x1_32( uint32_t * poly ,  uint32_t a_s3, uint32_t a_s2, uint32_t a_s1, uint32_t a_s0 , uint32_t idx_offset )
{
	uint32_t a_s3_0 = a_s3 ^ cantor_to_gf232_2x[ idx_offset>>3 ];
	uint32_t a_s2_0 = a_s2 ^ cantor_to_gf232_2x[ (idx_offset>>2) ];
	uint32_t a_s2_1 = a_s2 ^ cantor_to_gf232_2x[ (idx_offset>>2)+1 ];

	__m128i s1_a = _mm_set1_epi32( a_s1 ) ^ _mm_load_si128( (__m128i*)(cantor_to_gf232_2x + (idx_offset>>1)) );
	__m256i s0_a = _mm256_set1_epi32( a_s0 ) ^ _mm256_load_si256( (__m256i*)(cantor_to_gf232_2x + idx_offset) );

	__m256i p2 = _mm256_loadu_si256( (__m256i*)(poly) );
	__m256i p3 = _mm256_loadu_si256( (__m256i*)(poly+8) );
	__m256i p0,p1;

	p0 = _mm256_permute2x128_si256( p2 , p3 , 0x20 );
	p1 = _mm256_permute2x128_si256( p2 , p3 , 0x31 );

	p2 = _mm256_shuffle_pd( p0 , p1 , 0x0 );
	p3 = _mm256_shuffle_pd( p0 , p1 , 0xf );

	p0 = _mm256_blend_epi32( p2 , _mm256_slli_epi64(p3,32) , 0xaa );
	p1 = _mm256_blend_epi32( _mm256_srli_epi64(p2,32) , p3 , 0xaa );

	p1 ^= p0;
	p0 ^= gf232_mul8x8_avx2( p1 , s0_a );

	p2 = _mm256_blend_epi32( p0 , _mm256_slli_epi64(p1,32) , 0xaa );
	p3 = _mm256_blend_epi32( _mm256_srli_epi64(p0,32) , p1 , 0xaa );

	p3 ^= p2;
	p2 ^= gf232_mul8x4_avx2( p3 , s1_a );

	p0 = _mm256_shuffle_pd( p2 , p3 , 0x0 );
	p1 = _mm256_shuffle_pd( p2 , p3 , 0xf );

	p1 ^= p0;
	p0 ^= gf232_mul8x2_avx2( p1 , a_s2_0 , a_s2_1 );

	p2 = _mm256_permute2x128_si256( p0 , p1 , 0x20 );
	p3 = _mm256_permute2x128_si256( p0 , p1 , 0x31 );

	p3 ^= p2;
	p2 ^= gf232_mul8x1_avx2( p3 , a_s3_0 );

	_mm256_storeu_si256( (__m256i*)(poly) , p2 );
	_mm256_storeu_si256( (__m256i*)(poly+8) , p3 );

}


///////////////////////////////////////////



void ibtfy_32( uint32_t * poly , unsigned log_n , uint32_t idx_offset )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { i_btfy_s0_x1_32( poly , cantor_to_gf232(idx_offset) ); return; }
	if( 2 == log_n ) { i_btfy_s1s0_x1_32( poly , idx_offset ); return; }
	if( 3 == log_n ) { i_btfy_s2s1s0_x1_32( poly , idx_offset ); return; }
	if( 4 == log_n ) { i_btfy_s3s2s1s0_x1_32( poly, cantor_to_gf232(idx_offset>>3), cantor_to_gf232(idx_offset>>2), cantor_to_gf232(idx_offset>>1), cantor_to_gf232(idx_offset) , 0 ); return; }

    unsigned log_unit=4;
	do { //if(log_unit<log_n)
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));
		uint32_t s3_idx = cantor_to_gf232(idx_offset>>3);
		uint32_t s2_idx = cantor_to_gf232(idx_offset>>2);
		uint32_t s1_idx = cantor_to_gf232(idx_offset>>1);
		uint32_t s0_idx = cantor_to_gf232(idx_offset);

		uint32_t last_s0_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X/8,num);  // use 8 element per iteration.

		for(unsigned j=0;j<num;j+=tab_size) {
			uint32_t s0_diff = (j<<4)^last_s0_idx;
			last_s0_idx = (j<<4);
			s3_idx ^= cantor_to_gf232(s0_diff>>4);
			s2_idx ^= cantor_to_gf232(s0_diff>>3);
			s1_idx ^= cantor_to_gf232(s0_diff>>2);
			s0_idx ^= cantor_to_gf232(s0_diff>>1);
			for(unsigned k=0;k<tab_size;k++) {
				i_btfy_s3s2s1s0_x1_32( poly+(j+k)*unit , s3_idx, s2_idx, s1_idx, s0_idx , k*8 );  // use 8 element per iteration
			}
		}
	} while(0);
	for(log_unit=5; log_unit<=log_n; log_unit++) {
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));

		unsigned i = log_unit-1;
		uint64_t extra_a = cantor_to_gf232(idx_offset>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf232(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint64_t a = cantor_to_gf232_2x[k]^extra_a;
				i_butterfly_32_x8( poly + (j+k)*unit , unit , a );
			}
		}
	}
}



