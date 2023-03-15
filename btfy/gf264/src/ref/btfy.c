

#include <stdint.h>

#include <string.h>

#include "gf264.h"

#include "cantor_to_gf264.h"

#include "btfy.h"




//////////////////////////////////


static inline
void butterfly_64( uint64_t * poly , unsigned unit ,  uint64_t a )
{
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i++) {
		poly[i] ^= gf264_mul( poly[unit_2+i] , a );
		poly[unit_2+i] ^= poly[i];
	}
}


static inline
void i_butterfly_64( uint64_t * poly , unsigned unit , uint64_t a )
{
	unsigned unit_2= unit/2;
	for(unsigned i=0;i<unit_2;i++) {
		poly[unit_2+i] ^= poly[i];
		poly[i] ^= gf264_mul( poly[unit_2+i] , a );
	}
}

/////////////////////////////////////////////////////////

static inline unsigned min(unsigned a,unsigned b) { return (a<b)?a:b; }

void btfy_64( uint64_t * poly , unsigned log_n , uint64_t scalar_a )
{
	if( 0 == log_n ) return;

	for(unsigned log_unit=log_n; log_unit>0; log_unit--) {
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));

		unsigned i = log_unit-1;
		uint64_t extra_a = cantor_to_gf264(scalar_a>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {
			unsigned diff = j^last_idx;
			last_idx = j;
			extra_a ^= cantor_to_gf264(diff<<1);
			for(unsigned k=0;k<tab_size;k++) {
				uint64_t a = cantor_to_gf264_2x[k]^extra_a;
				butterfly_64( poly + (j+k)*unit , unit , a );
			}
		}
	}
}

void ibtfy_64( uint64_t * poly , unsigned log_n , uint64_t scalar_a )
{
	if( 0 == log_n ) return;

	for(unsigned log_unit=1; log_unit<=log_n; log_unit++) {
		unsigned unit = (1<<log_unit);
		unsigned num = (1<<(log_n-log_unit));

		unsigned i = log_unit-1;
		uint64_t extra_a = cantor_to_gf264(scalar_a>>i);

		unsigned last_idx = 0;
		unsigned tab_size = min(SIZE_TBL_CANTOR2X,num);

		for(unsigned j=0;j<num;j+=tab_size) {
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
