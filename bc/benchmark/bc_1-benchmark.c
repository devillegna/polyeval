
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



#include "bc_1.h"
#include "randombytes.h"




#define TESTSPEED

#ifdef TESTSPEED
#include "benchmark.h"
#endif



#define TEST_RUN 100

#define LEN (16384)
//#define LEN (8192)



void test( uint8_t * inp2 , uint8_t * inp1 , uint8_t * inp0 , unsigned len_byte )
{

#ifdef TESTSPEED
struct benchmark bm0;
bm_init(&bm0);
#endif
	int eq = 1;
	for(int i=1;i<=TEST_RUN;i++) {
		randombytes( (uint8_t*)inp0 , len_byte );

		memcpy( inp1 , inp0 , len_byte );

#ifdef TESTSPEED
BENCHMARK( bm0, {
#endif
		bc_1( inp1 , len_byte );
#ifdef TESTSPEED
} );
#endif

		memcpy( inp2 , inp1 , len_byte );

		ibc_1( inp2 , len_byte );

		if( !check_eq( (uint8_t*)inp0 , (uint8_t*)inp2 , len_byte ) ) {
			printf("neq: %d.\n", i );
			print_u8( inp0 , len_byte );
			printf("->\n");
			print_u8( inp1 , len_byte );
			printf("<-\n");
			print_u8( inp2 , len_byte );
			eq = 0;
			break;
		}
#if 0
		memcpy( inp8 , inp0 , sizeof(inp0) );
		bc_1_ref( inp8 );
		if( !check_eq( (uint8_t*)inp8 , (uint8_t*)inp1 , sizeof(inp1) ) ) {
			printf("neq: %d.\n", i );
			eq = 0;
			break;
		}
#endif

	}

	printf("\n%s for (i)bc of len %d.\n", (eq)?"OK":"ERROR" , len_byte );

#ifdef TESTSPEED
char bmmsg[256];
bm_dump(bmmsg,sizeof(bmmsg),&bm0);
printf("benchmark:\n%s\n\n", bmmsg );
#endif


}




int main()
{
	printf("benchmark (i)bc() for 256bytes to %dbytes.\n\n", LEN);

	uint8_t inp0[LEN] __attribute__ ((aligned (32)));
	uint8_t inp1[LEN] __attribute__ ((aligned (32)));
	uint8_t inp2[LEN] __attribute__ ((aligned (32)));


	for(unsigned i=256; i<=LEN; i <<= 1 ) {
		test( inp2 , inp1 , inp0 , i );
	}

	return 0;
}

