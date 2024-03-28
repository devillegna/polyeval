
//#include "randombytes.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>



void print_u8(const unsigned char *data, unsigned len )
{
	for(unsigned i=0;i<len;i++) {
		if( 0 == (i&15) ) printf("%3d: ",i);
		printf("%02x,", data[i] );
		if( 3 == (i&3) ) printf(".");
		if( 7 == (i&7) ) printf(" ");
		if( 15 == (i&15) ) printf("\n");
	}
}


void print_u16(const uint16_t *data, unsigned len )
{
	for(unsigned i=0;i<len;i++) {
		if( 0 == (i&7) ) printf("%3d: ",i);
		printf("0x%04x,", data[i] );
		if( 3 == (i&3) ) printf(" ");
		if( 7 == (i&7) ) printf("\n");
	}
}



void print_u32(const uint32_t *d32, unsigned l32 )
{
	for(unsigned i=0;i<l32;i++) {
		if( 0 == (i&7) ) printf("%3d: ",i);
		printf("0x%08x,", d32[i] );
		if( 3 == (i&3) ) printf(" ");
		if( 7 == (i&7) ) printf("\n");
	}
}


void print_u64(const uint64_t *data, unsigned l64 )
{
	unsigned len = l64;
	for(unsigned i=0;i<len;i++) {
		if( 0 == (i&7) ) printf("%3d: ",i);
		printf("0x%016llx,", data[i] );
		if( 3 == (i&3) ) printf(" ");
		if( 7 == (i&7) ) printf("\n");
	}
}



uint8_t check_eq( const uint8_t *vec0, const uint8_t *vec1, unsigned len)
{
	uint8_t diff = 0;
	for(unsigned i=0;i<len;i++){
		diff |= vec0[i]^vec1[i];
	}
	return diff==0;
}



#include "randombytes.h"





#define LOGLEN 8
#if LOGLEN < 2
error -- length must >=  4xu64_t
#endif

#define LEN (1<<LOGLEN)

#define TEST_RUN 100


#include "gf232.h"
#include "cantor_to_gf232.h"





uint32_t polyeval_ref( const uint32_t *poly, unsigned len , uint32_t v ) {
  uint32_t r = poly[len-1];
  for(int i=(int)(len-2);i>=0;i--) {
    r = gf232_mul(r,v)^poly[i];
  }
  return r;
}


void eval_ref( uint32_t *r , const uint32_t *poly, unsigned len , uint32_t extra_a ) {
  for(unsigned i=0;i<len;i++) {
    r[i] = polyeval_ref( poly , len , cantor_to_gf232(extra_a^i) );
  }
}




#include "bc_32.h"
#include "btfy32.h"

void eval_fft( uint32_t *r , const uint32_t *poly, unsigned len , uint32_t extra_a ) {
  for(unsigned i=0;i<len;i++) { r[i] = poly[i]; }
  bc_32(r,len);
  btfy_32(r,__builtin_ctz( len ),extra_a);
}


int test_equal( uint32_t * poly0, uint32_t * eval0, uint32_t * eval1, unsigned len )
{
	printf("\npoly len: %d (uint32_t)\n", len );

	int eq = 1;
	for(int i=1;i<=TEST_RUN;i++) {
		randombytes( (uint8_t*)poly0 , len*4 );
		eval_ref( eval0 , poly0 , len , 0 );

		eval_fft( eval1 , poly0 , len , 0 );

		if( !check_eq( (uint8_t*)eval0 , (uint8_t*)eval1 , len*4 ) ) {
			printf("neq: %d.\n", i );
			print_u32( poly0 , len );
			printf("->\n");
			print_u32( eval0 , len );
			printf("->\n");
			print_u32( eval1 , len );
			eq = 0;
			break;
		}
	}

	printf("\n%s for %d tests.\n", (eq)?"OK":"ERROR" , TEST_RUN );

    return (eq)?0:-1;
}


int main(void)
{
	printf("GF(2^64) polynomial [%d] evaluation test.\n\n", LEN);

	uint32_t poly0[LEN] __attribute__ ((aligned (32)));
	uint32_t eval0[LEN] __attribute__ ((aligned (32)));
	uint32_t eval1[LEN] __attribute__ ((aligned (32)));

	for(unsigned i=0;i<=LOGLEN;i++){
		if( test_equal( poly0, eval0, eval1, 1<<i ) ) break;
	}

	return 0;
}

