

#include <stdint.h>
#include <string.h>

#include "btfy32.h"
#include "cantor_to_gf232.h"

#include "gf232_aesni.h"




#if 1



static inline
void butterfly_32( uint32_t * poly , unsigned unit ,  uint32_t a )
{
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i++) {
		poly[i] ^= gf232_mul_u32_sse( poly[unit_2+i] , a );
		poly[unit_2+i] ^= poly[i];
	}
}


static inline
void i_butterfly_32( uint32_t * poly , unsigned unit , uint32_t a )
{
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i++) {
		poly[unit_2+i] ^= poly[i];
		poly[i] ^= gf232_mul_u32_sse( poly[unit_2+i] , a );
	}
}

/////////////////////////////////////////////////////////

static inline unsigned min(unsigned a,unsigned b) { return (a<b)?a:b; }

void btfy_32( uint32_t * poly , unsigned log_n , uint32_t scalar_a )
{
	if( 0 == log_n ) return;

	for(unsigned log_unit=log_n; log_unit>0; log_unit--) {
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));

		unsigned i = log_unit-1;
		uint32_t extra_a = cantor_to_gf232(scalar_a>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf232(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint64_t a = cantor_to_gf232_2x[k]^extra_a;
				butterfly_32( poly + (j+k)*unit , unit , a );
			}
		}
	}
}

void ibtfy_32( uint32_t * poly , unsigned log_n , uint32_t scalar_a )
{
	if( 0 == log_n ) return;

	for(unsigned log_unit=1; log_unit<=log_n; log_unit++) {
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));

		unsigned i = log_unit-1;
		uint64_t extra_a = cantor_to_gf232(scalar_a>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf232(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint64_t a = cantor_to_gf232_2x[k]^extra_a;
				i_butterfly_32( poly + (j+k)*unit , unit , a );
			}
		}
	}
}




#else

//
// This macro increase performance but reorder the outputs
//
//#define _REORDER_OUTPUT_


/////////////////////////////////////////////////////////////




static inline
void butterfly_64_avx2( __m256i * poly , unsigned unit_256 ,  __m128i a ) /// unit_256 >= 2
{
	unsigned unit_2= unit_256/2;
	for(unsigned i=0;i<unit_2;i++) {
		__m256i poly_i = _mm256_loadu_si256( poly + i );
		__m256i poly_unit2_i = _mm256_loadu_si256( poly + unit_2 + i );
		_mm_prefetch( &poly[i+1] , _MM_HINT_T0 );
		_mm_prefetch( &poly[unit_2+i+1] , _MM_HINT_T0 );

		poly_i ^= _gf2ext64_mul_4x1_avx2( poly_unit2_i , a );
		poly_unit2_i ^= poly_i;

		_mm256_storeu_si256( poly+i , poly_i );
		_mm256_storeu_si256( poly+unit_2+i , poly_unit2_i );
		// slower:
		//_mm256_stream_si256( poly+i , poly_i );
		//_mm256_stream_si256( poly+unit_2+i , poly_unit2_i );
	}
}




//////////////////////////////////


static inline
void i_butterfly_64_avx2( __m256i * poly , unsigned unit_256 , __m128i a ) /// unit_256 >= 2
{
	unsigned unit_2= unit_256/2;
	for(unsigned i=0;i<unit_2;i++) {
		__m256i poly_i = _mm256_loadu_si256( poly + i );
		__m256i poly_unit2_i = _mm256_loadu_si256( poly + unit_2 + i );
		_mm_prefetch( &poly[i+1] , _MM_HINT_T0 );
		_mm_prefetch( &poly[unit_2+i+1] , _MM_HINT_T0 );

		poly_unit2_i ^= poly_i;
		poly_i ^= _gf2ext64_mul_4x1_avx2( poly_unit2_i , a );

		_mm256_storeu_si256( poly+i , poly_i );
		_mm256_storeu_si256( poly+unit_2+i , poly_unit2_i );
	}
}



///////////////////////////////////////////


static inline
void btfy_s1_x1( uint64_t * poly4 , uint64_t extra_a ) {
	__m128i p01 = _mm_loadu_si128( (__m128i*)poly4 );
	__m128i p23 = _mm_loadu_si128( (__m128i*)(poly4+2) );
	__m128i eaea = _mm_set1_epi64x( extra_a );

	p01 ^= _gf2ext64_mul_2x1_sse( p23 , eaea );
	p23 ^= p01;

	_mm_storeu_si128( (__m128i*)poly4 , p01 );
	_mm_storeu_si128( (__m128i*)(poly4+2) , p23 );
}

static inline
void i_btfy_s1_x1( uint64_t * poly4 , uint64_t extra_a ) {
	__m128i p01 = _mm_loadu_si128( (__m128i*)poly4 );
	__m128i p23 = _mm_loadu_si128( (__m128i*)(poly4+2) );
	__m128i eaea = _mm_set1_epi64x( extra_a );

	p23 ^= p01;
	p01 ^= _gf2ext64_mul_2x1_sse( p23 , eaea );

	_mm_storeu_si128( (__m128i*)poly4 , p01 );
	_mm_storeu_si128( (__m128i*)(poly4+2) , p23 );
}

static inline
void btfy_s0_x1( uint64_t * poly2 , uint64_t extra_a ) {
	__m128i p01 = _mm_loadu_si128( (__m128i*)poly2 );
	__m128i eaea = _mm_set1_epi64x( extra_a );

	p01 ^= _gf2ext64_mul_hi_sse( p01 , eaea );
	p01 ^= _mm_slli_si128( p01 , 8 );

	_mm_storeu_si128( (__m128i*)poly2 , p01 );
}

static inline
void i_btfy_s0_x1( uint64_t * poly2 , uint64_t extra_a ) {
	__m128i p01 = _mm_loadu_si128( (__m128i*)poly2 );
	__m128i eaea = _mm_set1_epi64x( extra_a );

	p01 ^= _mm_slli_si128( p01 , 8 );
	p01 ^= _gf2ext64_mul_hi_sse( p01 , eaea );

	_mm_storeu_si128( (__m128i*)poly2 , p01 );
}


///////////////////////////////////////////

static inline unsigned min(unsigned a,unsigned b) { return (a<b)?a:b; }


static
void btfy_s1( __m256i* ptr , unsigned n_m256, uint64_t extra_a )
{
	unsigned num = n_m256;   // one 256-bit data comsumes 1 elements in the table
	unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

	unsigned last_idx = 0;
	for(unsigned j=0;j<num;j+=tab_size) {
		unsigned diff = j^last_idx;
		last_idx = j;
		extra_a ^= cantor_to_gf264(diff<<1);
		__m128i s1_ea = _mm_set1_epi64x(extra_a);

		for(unsigned k=0;k<tab_size;k+=2) {
			__m128i s1_a = _mm_load_si128( (__m128i*) &cantor_to_gf264_2x[k] )^s1_ea;

			__m256i t0123 = _mm256_loadu_si256(ptr);
			__m256i t4567 = _mm256_loadu_si256(ptr+1);

			__m256i t0145 = _mm256_permute2x128_si256( t0123 , t4567 , 0x20 );
			__m256i t2367 = _mm256_permute2x128_si256( t0123 , t4567 , 0x31 );

			t0145 ^= _gf2ext64_mul_2x1_x2_avx2( t2367 , s1_a );
			t2367 ^= t0145;

			_mm256_storeu_si256( ptr , t0145 );
			_mm256_storeu_si256( ptr+1 , t2367 );
			ptr += 2;
		}
	}
}


static
void btfy_s0( __m256i* ptr , unsigned n_m256, uint64_t extra_a )
{
	unsigned num = n_m256*2;   // one 256-bit data comsumes 2 elements in the table
	unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

	unsigned last_idx = 0;
	for(unsigned j=0;j<num;j+=tab_size) {
		unsigned diff = j^last_idx;
		last_idx = j;
		extra_a ^= cantor_to_gf264(diff<<1);
		__m256i s0_ea = _mm256_broadcastq_epi64( _mm_set1_epi64x(extra_a) );

		for(unsigned k=0;k<tab_size;k+=4) {
			__m256i s0_a = _mm256_load_si256( (__m256i*) &cantor_to_gf264_2x[k] )^s0_ea;

			__m256i t0145 = _mm256_loadu_si256(ptr);
			__m256i t2367 = _mm256_loadu_si256(ptr+1);

			__m256i t0246 = _mm256_unpacklo_epi64( t0145 , t2367 );
			__m256i t1357 = _mm256_unpackhi_epi64( t0145 , t2367 );

			t0246 ^= _gf2ext64_mul_4x4_avx2( t1357 , s0_a );
			t1357 ^= t0246;
#if defined(_REORDER_OUTPUT_)
			_mm256_storeu_si256( ptr , t0246 );
			_mm256_storeu_si256( ptr+1 , t1357 );
#else
			t0145 = _mm256_unpacklo_epi64( t0246 , t1357 );
			t2367 = _mm256_unpackhi_epi64( t0246 , t1357 );

			__m256i t0123 = _mm256_permute2x128_si256( t0145 , t2367 , 0x20 );
			__m256i t4567 = _mm256_permute2x128_si256( t0145 , t2367 , 0x31 );

			_mm256_storeu_si256( ptr , t0123 );
			_mm256_storeu_si256( ptr+1 , t4567 );
#endif
			ptr += 2;
		}
	}
}


/////////////////////////////////////////////////////


static
void i_btfy_s1( __m256i* ptr , unsigned n_m256, uint64_t extra_a )
{
	unsigned num = n_m256;   // one 256-bit data comsumes 1 elements in the table
	unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

	unsigned last_idx = 0;
	for(unsigned j=0;j<num;j+=tab_size) {
		unsigned diff = j^last_idx;
		last_idx = j;
		extra_a ^= cantor_to_gf264(diff<<1);
		__m128i s1_ea = _mm_set1_epi64x(extra_a);

		for(unsigned k=0;k<tab_size;k+=2) {
			__m128i s1_a = _mm_load_si128( (__m128i*) &cantor_to_gf264_2x[k] )^s1_ea;

			__m256i t0145 = _mm256_loadu_si256(ptr);
			__m256i t2367 = _mm256_loadu_si256(ptr+1);

			t2367 ^= t0145;
			t0145 ^= _gf2ext64_mul_2x1_x2_avx2( t2367 , s1_a );

			__m256i t0123 = _mm256_permute2x128_si256( t0145 , t2367 , 0x20 );
			__m256i t4567 = _mm256_permute2x128_si256( t0145 , t2367 , 0x31 );

			_mm256_storeu_si256( ptr , t0123 );
			_mm256_storeu_si256( ptr+1 , t4567 );
			ptr += 2;
		}
	}
}

static
void i_btfy_s0( __m256i* ptr , unsigned n_m256, uint64_t extra_a )
{
	unsigned num = n_m256*2;   // one 256-bit data comsumes 2 elements in the table
	unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

	unsigned last_idx = 0;
	for(unsigned j=0;j<num;j+=tab_size) {
		unsigned diff = j^last_idx;
		last_idx = j;
		extra_a ^= cantor_to_gf264(diff<<1);
		__m256i s0_ea = _mm256_broadcastq_epi64( _mm_set1_epi64x(extra_a) );

		for(unsigned k=0;k<tab_size;k+=4) {
			__m256i s0_a = _mm256_load_si256( (__m256i*) &cantor_to_gf264_2x[k] )^s0_ea;

#if defined(_REORDER_OUTPUT_)
			__m256i t0246 = _mm256_loadu_si256(ptr);
			__m256i t1357 = _mm256_loadu_si256(ptr+1);
			__m256i t0145, t2367;
#else
			__m256i t0123 = _mm256_loadu_si256(ptr);
			__m256i t4567 = _mm256_loadu_si256(ptr+1);

			__m256i t0145 = _mm256_permute2x128_si256( t0123 , t4567 , 0x20 );
			__m256i t2367 = _mm256_permute2x128_si256( t0123 , t4567 , 0x31 );

			__m256i t0246 = _mm256_unpacklo_epi64( t0145 , t2367 );
			__m256i t1357 = _mm256_unpackhi_epi64( t0145 , t2367 );
#endif

			t1357 ^= t0246;
			t0246 ^= _gf2ext64_mul_4x4_avx2( t1357 , s0_a );

			t0145 = _mm256_unpacklo_epi64( t0246 , t1357 );
			t2367 = _mm256_unpackhi_epi64( t0246 , t1357 );

			_mm256_storeu_si256( ptr , t0145 );
			_mm256_storeu_si256( ptr+1 , t2367 );
			ptr += 2;
		}
	}
}



/////////////////////////////////////////////////////////





void btfy_64( uint64_t * poly , unsigned log_n , uint64_t scalar_a )
{
	if( 0 == log_n ) { return; }
	if( 1 == log_n ) { btfy_s0_x1( poly , cantor_to_gf264(scalar_a) ); return; }
	if( 2 == log_n ) {
		btfy_s1_x1( poly , cantor_to_gf264(scalar_a>>1) );
		uint64_t extra_a = cantor_to_gf264(scalar_a);
		btfy_s0_x1( poly   , extra_a );
		btfy_s0_x1( poly+2 , extra_a^cantor_basis[1] );
		return;
	}

	// now log_n > 2

	unsigned n_terms = 1<<log_n;
	unsigned i=log_n-1;

	__m128i zero = _mm_setzero_si128();
	for( ; i>1; i-- ) {
		unsigned unit = 1<<(i+1); /// u >= 8
		unsigned num  = 1<<(log_n-i-1);

		uint64_t extra_a = cantor_to_gf264(scalar_a>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {  // process constants outside of the cantor_to_gf264_2x[]
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf264(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint64_t a = cantor_to_gf264_2x[k]^extra_a;
				__m128i xmm_a = _mm_insert_epi64( zero , a , 0 );

				butterfly_64_avx2( (__m256i*)(poly + (j+k)*unit) , unit>>2 , xmm_a );
			}
		}
	}

	btfy_s1( (__m256i*)poly , n_terms/4 , cantor_to_gf264(scalar_a>>1) );
	btfy_s0( (__m256i*)poly , n_terms/4 , cantor_to_gf264(scalar_a) );
}

void ibtfy_64( uint64_t * poly , unsigned log_n , uint64_t scalar_a )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { i_btfy_s0_x1( poly , cantor_to_gf264(scalar_a) ); return; }
	if( 2 == log_n ) {
		uint64_t extra_a = cantor_to_gf264(scalar_a);
		i_btfy_s0_x1( poly   , extra_a );
		i_btfy_s0_x1( poly+2 , extra_a^cantor_basis[1] );
		i_btfy_s1_x1( poly , cantor_to_gf264(scalar_a>>1) );
		return;
	}

	// now, log_n > 2

	unsigned n_terms = 1<<log_n;
	unsigned i=0;

	__m128i zero = _mm_setzero_si128();
	i_btfy_s0( (__m256i*)poly , n_terms/4 , cantor_to_gf264(scalar_a) );
	i_btfy_s1( (__m256i*)poly , n_terms/4 , cantor_to_gf264(scalar_a>>1) );
	i = 2;

	for( ; i<log_n; i++ ) {
		unsigned unit = 1<<(i+1);
		unsigned num  = 1<<(log_n-i-1);

		uint64_t extra_a = cantor_to_gf264(scalar_a>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {  // process constants outside of the cantor_to_gf264_2x[]
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf264(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint64_t a = cantor_to_gf264_2x[k]^extra_a;
				__m128i xmm_a = _mm_insert_epi64( zero , a , 0 );

				i_butterfly_64_avx2( (__m256i*)(poly + (j+k)*unit) , unit>>2 , xmm_a );
			}
		}
	}
}

#endif
