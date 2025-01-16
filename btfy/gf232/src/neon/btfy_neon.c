

#include <stdint.h>
#include <string.h>


#include "cantor_to_gf232.h"
#include "btfy32.h"
#include "gf232_neon.h"

static inline unsigned min(unsigned a,unsigned b) { return (a<b)?a:b; }



static inline
void butterfly_32_x4( uint32_t * poly , unsigned unit ,  uint32_t a )
{
	unsigned unit_2= unit/2;  // ASSERT  (uint_2 & 3) == 0
	for(unsigned i=0;i<unit_2;i+=4) {
		uint32x4_t p0 = vld1q_u32( poly+i );
		uint32x4_t p1 = vld1q_u32( poly+unit_2+i );
		p0 ^= _gf232_mul4x1_neon( p1 , a );
		p1 ^= p0;
		vst1q_u32( poly+i , p0 );
		vst1q_u32( poly+unit_2+i , p1 );
	}
}


static inline
void btfy_s0_x1_32( uint32_t * poly , uint32_t a )
{
	poly[0] ^= _gf232_mulx1_neon( poly[1] , a );
	poly[1] ^= poly[0];
}

static inline
void btfy_s1s0_x1_32( uint32_t * poly , uint32_t idx_offset )
{
	uint32_t s1_a = cantor_to_gf232( idx_offset>>1 );
	uint32_t s0_a_0 = cantor_to_gf232( idx_offset );
	uint32_t s0_a_1 = s0_a_0 ^ cantor_basis32[1];

	uint32_t p0 = poly[0];
	uint32_t p1 = poly[1];
	uint32_t p2 = poly[2];
	uint32_t p3 = poly[3];

	p0 ^= _gf232_mulx1_neon( p2 , s1_a );
	p1 ^= _gf232_mulx1_neon( p3 , s1_a );
	p2 ^= p0;
	p3 ^= p1;
	p0 ^= _gf232_mulx1_neon( p1 , s0_a_0 );
	p1 ^= p0;
	p2 ^= _gf232_mulx1_neon( p3 , s0_a_1 );
	p3 ^= p2;

	poly[0] = p0;
	poly[1] = p1;
	poly[2] = p2;
	poly[3] = p3;
}

static inline
void btfy_s2s1s0_x1_32( uint32_t * poly , uint32_t a_s2, uint32_t a_s1, uint32_t a_s0 , uint32_t idx_offset )
{
	uint32_t s2_a_0 = a_s2 ^ cantor_to_gf232_2x[ (idx_offset>>2) ];
	uint32x4_t s1_a = vdupq_n_u32( a_s1 ) ^ vcombine_u32( vdup_n_u32(cantor_to_gf232_2x[ (idx_offset>>1) ]),vdup_n_u32(cantor_to_gf232_2x[ (idx_offset>>1)+1 ]) );
	uint32x4_t s0_a = vdupq_n_u32( a_s0 ) ^ vld1q_u32( cantor_to_gf232_2x + idx_offset );

	uint32x4_t p0 = vld1q_u32( poly );
	uint32x4_t p1 = vld1q_u32( poly+4 );

	p0 ^= _gf232_mul4x1_neon( p1 , s2_a_0 );
	p1 ^= p0;

	uint32x4_t p2 = vreinterpretq_u32_u64(vuzp1q_u64( vreinterpretq_u64_u32(p0) , vreinterpretq_u64_u32(p1) ));
	uint32x4_t p3 = vreinterpretq_u32_u64(vuzp2q_u64( vreinterpretq_u64_u32(p0) , vreinterpretq_u64_u32(p1) ));
	p2 ^= _gf232_mul4x4_neon( p3 , s1_a );
	p3 ^= p2;
	p0 = vreinterpretq_u32_u64(vuzp1q_u64( vreinterpretq_u64_u32(p2) , vreinterpretq_u64_u32(p3) ));
	p1 = vreinterpretq_u32_u64(vuzp2q_u64( vreinterpretq_u64_u32(p2) , vreinterpretq_u64_u32(p3) ));

	p2 = vuzp1q_u32( p0 , p1 );
	p3 = vuzp2q_u32( p0 , p1 );
	p2 ^= _gf232_mul4x4_neon( p3 , s0_a );
	p3 ^= p2;
	p0 = vzip1q_u32( p2 , p3 );
	p1 = vzip2q_u32( p2 , p3 );

	vst1q_u32( poly , p0 );
	vst1q_u32( poly+4 , p1 );
}




/////////////////////////////////////////////////////////


void btfy_32( uint32_t * poly , unsigned log_n , uint32_t idx_offset )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { btfy_s0_x1_32( poly , cantor_to_gf232(idx_offset) ); return; }
	if( 2 == log_n ) { btfy_s1s0_x1_32( poly , idx_offset ); return; }
	if( 3 == log_n ) { btfy_s2s1s0_x1_32( poly , cantor_to_gf232(idx_offset>>2) , cantor_to_gf232(idx_offset>>1) , cantor_to_gf232(idx_offset) , 0 ); return; }

    unsigned log_unit=log_n;
	for(;log_unit>3; log_unit--) {
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
				butterfly_32_x4( poly + (j+k)*unit , unit , a );
			}
		}
	}
	if(3==log_unit){
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));
		uint32_t s2_idx = cantor_to_gf232(idx_offset>>2);
		uint32_t s1_idx = cantor_to_gf232(idx_offset>>1);
		uint32_t s0_idx = cantor_to_gf232(idx_offset);

		uint32_t last_s0_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X/4,num);  // use 4 element per iteration.

		for(unsigned j=0;j<num;j+=tab_size) {
			uint32_t s0_diff = (j<<3)^last_s0_idx;
			last_s0_idx = (j<<3);
			s2_idx ^= cantor_to_gf232(s0_diff>>3);
			s1_idx ^= cantor_to_gf232(s0_diff>>2);
			s0_idx ^= cantor_to_gf232(s0_diff>>1);
			for(unsigned k=0;k<tab_size;k++) {
				btfy_s2s1s0_x1_32( poly+(j+k)*unit , s2_idx, s1_idx, s0_idx , k*4 );  // use 4 element per iteration
			}
		}
	}
}



/////////////////////////////////////////////////////////

static inline
void i_butterfly_32_x4( uint32_t * poly , unsigned unit ,  uint32_t a )
{
	unsigned unit_2= unit/2;  // ASSERT  (uint_2 & 3) == 0
	for(unsigned i=0;i<unit_2;i+=4) {
		uint32x4_t p0 = vld1q_u32( poly+i );
		uint32x4_t p1 = vld1q_u32( poly+unit_2+i );
		p1 ^= p0;
		p0 ^= _gf232_mul4x1_neon( p1 , a );
		vst1q_u32( poly+i , p0 );
		vst1q_u32( poly+unit_2+i , p1 );
	}
}

static inline
void i_btfy_s0_x1_32( uint32_t * poly , uint32_t a )
{
	poly[1] ^= poly[0];
	poly[0] ^= _gf232_mulx1_neon( poly[1] , a );
}

static inline
void i_btfy_s1s0_x1_32( uint32_t * poly , uint32_t idx_offset )
{
	uint32_t s1_a = cantor_to_gf232( idx_offset>>1 );
	uint32_t s0_a_0 = cantor_to_gf232( idx_offset );
	uint32_t s0_a_1 = s0_a_0 ^ cantor_basis32[1];

	uint32_t p0 = poly[0];
	uint32_t p1 = poly[1];
	uint32_t p2 = poly[2];
	uint32_t p3 = poly[3];

	p1 ^= p0;
	p0 ^= _gf232_mulx1_neon( p1 , s0_a_0 );
	p3 ^= p2;
	p2 ^= _gf232_mulx1_neon( p3 , s0_a_1 );

	p2 ^= p0;
	p3 ^= p1;
	p0 ^= _gf232_mulx1_neon( p2 , s1_a );
	p1 ^= _gf232_mulx1_neon( p3 , s1_a );

	poly[0] = p0;
	poly[1] = p1;
	poly[2] = p2;
	poly[3] = p3;
}

static inline
void i_btfy_s2s1s0_x1_32( uint32_t * poly , uint32_t a_s2, uint32_t a_s1, uint32_t a_s0 , uint32_t idx_offset )
{
	uint32_t s2_a_0 = a_s2 ^ cantor_to_gf232_2x[ (idx_offset>>2) ];
	uint32x4_t s1_a = vdupq_n_u32( a_s1 ) ^ vcombine_u32( vdup_n_u32(cantor_to_gf232_2x[ (idx_offset>>1) ]),vdup_n_u32(cantor_to_gf232_2x[ (idx_offset>>1)+1 ]) );
	uint32x4_t s0_a = vdupq_n_u32( a_s0 ) ^ vld1q_u32( cantor_to_gf232_2x + idx_offset );

	uint32x4_t p0,p1,p2,p3;
	p0 = vld1q_u32( poly );
	p1 = vld1q_u32( poly+4 );

	p2 = vuzp1q_u32( p0 , p1 );
	p3 = vuzp2q_u32( p0 , p1 );
	p3 ^= p2;
	p2 ^= _gf232_mul4x4_neon( p3 , s0_a );
	p0 = vzip1q_u32( p2 , p3 );
	p1 = vzip2q_u32( p2 , p3 );

	p2 = vreinterpretq_u32_u64(vuzp1q_u64( vreinterpretq_u64_u32(p0) , vreinterpretq_u64_u32(p1) ));
	p3 = vreinterpretq_u32_u64(vuzp2q_u64( vreinterpretq_u64_u32(p0) , vreinterpretq_u64_u32(p1) ));
	p3 ^= p2;
	p2 ^= _gf232_mul4x4_neon( p3 , s1_a );
	p0 = vreinterpretq_u32_u64(vuzp1q_u64( vreinterpretq_u64_u32(p2) , vreinterpretq_u64_u32(p3) ));
	p1 = vreinterpretq_u32_u64(vuzp2q_u64( vreinterpretq_u64_u32(p2) , vreinterpretq_u64_u32(p3) ));

	p1 ^= p0;
	p0 ^= _gf232_mul4x1_neon( p1 , s2_a_0 );

	vst1q_u32( poly , p0 );
	vst1q_u32( poly+4 , p1 );
}


/////////////////////////////////////////////////////////


void ibtfy_32( uint32_t * poly , unsigned log_n , uint32_t idx_offset )
{
	if( 0 == log_n ) return;
	if( 1 == log_n ) { i_btfy_s0_x1_32( poly , cantor_to_gf232(idx_offset) ); return; }
	if( 2 == log_n ) { i_btfy_s1s0_x1_32( poly , idx_offset ); return; }
	if( 3 == log_n ) { i_btfy_s2s1s0_x1_32( poly , cantor_to_gf232(idx_offset>>2) , cantor_to_gf232(idx_offset>>1) , cantor_to_gf232(idx_offset) , 0 ); return; }

    unsigned log_unit=3;
	do { //if(log_unit<log_n)
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));
		uint32_t s2_idx = cantor_to_gf232(idx_offset>>2);
		uint32_t s1_idx = cantor_to_gf232(idx_offset>>1);
		uint32_t s0_idx = cantor_to_gf232(idx_offset);

		uint32_t last_s0_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X/4,num);  // use 4 element per iteration.

		for(unsigned j=0;j<num;j+=tab_size) {
			uint32_t s0_diff = (j<<3)^last_s0_idx;
			last_s0_idx = (j<<3);
			s2_idx ^= cantor_to_gf232(s0_diff>>3);
			s1_idx ^= cantor_to_gf232(s0_diff>>2);
			s0_idx ^= cantor_to_gf232(s0_diff>>1);
			for(unsigned k=0;k<tab_size;k++) {
				i_btfy_s2s1s0_x1_32( poly+(j+k)*unit , s2_idx, s1_idx, s0_idx , k*4 );  // use 4 element per iteration
			}
		}
	} while(0);
	for(log_unit=4; log_unit<=log_n; log_unit++) {
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
				i_butterfly_32_x4( poly + (j+k)*unit , unit , a );
			}
		}
	}
}
