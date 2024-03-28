

#include <stdint.h>
#include <string.h>


#include "cantor_to_gf232.h"
#include "btfy32.h"
#include "gf232_neon.h"




#if 1


static inline
void butterfly_32( uint32_t * poly , unsigned unit ,  uint32_t a )
{
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i++) {
		poly[i] ^= _gf232_mulx1_neon( poly[unit_2+i] , a );
		poly[unit_2+i] ^= poly[i];
	}
}


static inline
void i_butterfly_32( uint32_t * poly , unsigned unit , uint32_t a )
{
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i++) {
		poly[unit_2+i] ^= poly[i];
		poly[i] ^= _gf232_mulx1_neon( poly[unit_2+i] , a );
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

//////////////////////////////////


static inline
void butterfly_64( uint64_t * poly , unsigned unit ,  uint64_t a )
{
	uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);
	uint64x2_t aa = vdupq_n_u64(a);
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i+=2) {
		uint64x2_t p0 = vld1q_u64( poly+i );
		uint64x2_t p1 = vld1q_u64( poly+unit_2+i );
		p0 ^= _gf264_mul_neon( aa , p1 , mask_0x1b );
		p1 ^= p0;
		vst1q_u64( poly+i , p0 );
		vst1q_u64( poly+unit_2+i , p1 );
	}
}


static inline
void i_butterfly_64( uint64_t * poly , unsigned unit , uint64_t a )
{
	uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);
	uint64x2_t aa = vdupq_n_u64(a);
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i+=2) {
		uint64x2_t p0 = vld1q_u64( poly+i );
		uint64x2_t p1 = vld1q_u64( poly+unit_2+i );
		p1 ^= p0;
		p0 ^= _gf264_mul_neon( aa , p1 , mask_0x1b );
		vst1q_u64( poly+i , p0 );
		vst1q_u64( poly+unit_2+i , p1 );
	}
}

/////////////////////////////////////////////////////////

static inline
void btfy_s0( uint64_t * poly , uint64_t scalar_a )
{
	uint64_t extra_a = cantor_to_gf264(scalar_a);
	uint64_t p0 = poly[0];
	uint64_t p1 = poly[1];
	p0 ^= _gf264_mulx1_neon( p1 , extra_a );
	p1 ^= p0;
	poly[0] = p0;
	poly[1] = p1;
}

static inline
void i_btfy_s0( uint64_t * poly , uint64_t scalar_a )
{
	uint64_t extra_a = cantor_to_gf264(scalar_a);
	uint64_t p0 = poly[0];
	uint64_t p1 = poly[1];
	p1 ^= p0;
	p0 ^= _gf264_mulx1_neon( p1 , extra_a );
	poly[0] = p0;
	poly[1] = p1;
}

///////////////////////////////

static inline unsigned min(unsigned a,unsigned b) { return (a<b)?a:b; }

static
void btfy_stage_s0_2x( uint64_t * poly , unsigned n_terms , uint64_t scalar_a )
{
	uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);
	uint64x2_t extra_a = vdupq_n_u64( cantor_to_gf264(scalar_a) );
	unsigned num = n_terms/2;
	unsigned last_idx = 0;
	unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

	for(unsigned j=0;j<num;j+=tab_size) {
		unsigned diff = j^last_idx;
		last_idx = j;
		extra_a ^= vdupq_n_u64(cantor_to_gf264(diff<<1));

		for(unsigned k=0;k<tab_size;k+=2) {
			uint64x2x2_t p = vld2q_u64( poly+(j+k)*2 );
			uint64x2_t a = vld1q_u64( &cantor_to_gf264_2x[k] ) ^ extra_a;
			p.val[0] ^= _gf264_mul_neon( p.val[1] , a , mask_0x1b );
			p.val[1] ^= p.val[0];
			vst2q_u64( poly+(j+k)*2 , p );
		}
	}
}

static
void i_btfy_stage_s0_2x( uint64_t * poly , unsigned n_terms , uint64_t scalar_a )
{
	uint64x2_t mask_0x1b = vdupq_n_u64(0x1b);
	uint64x2_t extra_a = vdupq_n_u64( cantor_to_gf264(scalar_a) );
	unsigned num = n_terms/2;
	unsigned last_idx = 0;
	unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

	for(unsigned j=0;j<num;j+=tab_size) {
		unsigned diff = j^last_idx;
		last_idx = j;
		extra_a ^= vdupq_n_u64(cantor_to_gf264(diff<<1));

		for(unsigned k=0;k<tab_size;k+=2) {
			uint64x2x2_t p = vld2q_u64( poly+(j+k)*2 );
			uint64x2_t a = vld1q_u64( &cantor_to_gf264_2x[k] ) ^ extra_a;
			p.val[1] ^= p.val[0];
			p.val[0] ^= _gf264_mul_neon( p.val[1] , a , mask_0x1b );
			vst2q_u64( poly+(j+k)*2 , p );
		}
	}
}

//////////////////////////////////////////////////////


void btfy_64( uint64_t * poly , unsigned log_n , uint64_t scalar_a )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { btfy_s0(poly,scalar_a); return; }

    unsigned n_terms = 1<<log_n;

	for(unsigned i=log_n-1; i>0; i--) {
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
				butterfly_64( poly + (j+k)*unit , unit , a );
			}
		}
	}
	btfy_stage_s0_2x( poly , n_terms , scalar_a );
}

void ibtfy_64( uint64_t * poly , unsigned log_n , uint64_t scalar_a )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { i_btfy_s0(poly,scalar_a); return; }

    unsigned n_terms = 1<<log_n;

	i_btfy_stage_s0_2x( poly , n_terms , scalar_a );

	for(unsigned i=1; i<log_n; i++) {
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
				i_butterfly_64( poly + (j+k)*unit , unit , a );
			}
		}
	}
}


#endif
