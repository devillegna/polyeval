
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>


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
		printf("0x%016lx,", data[i] );
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



#include "btfy.h"
#include "randombytes.h"

#define TEST_RUN 100


#define LOG_LEN  15
#define LEN    (1<<LOG_LEN)



#define TESTSPEED

#ifdef TESTSPEED
#include "benchmark.h"
#endif




int main(void)
{

	uint64_t inp0[LEN] __attribute__ ((aligned (32)));
	uint64_t inp1[LEN] __attribute__ ((aligned (32)));
	uint64_t inp2[LEN] __attribute__ ((aligned (32)));


	uint8_t eq = 1;

#ifdef TESTSPEED
struct benchmark bm0;
bm_init(&bm0);
#endif

	printf("btfy( %d x u64 ).\n", LEN );

	for(int i=1;i<=TEST_RUN;i++) {
		randombytes( (uint8_t*)inp0 , sizeof(inp0) );
		if(0==i) {
			for(unsigned j=0;j<sizeof(inp0);j++) inp0[j]=0xffffffff;
		}
		memcpy( inp1 , inp0 , sizeof(inp0) );

#ifdef TESTSPEED
BENCHMARK( bm0, {
#endif
		btfy_64(inp1,LOG_LEN,0);
#ifdef TESTSPEED
} );
#endif

		memcpy( inp2 , inp1 , sizeof(inp1) );

		ibtfy_64(inp2,LOG_LEN,0);

		if(0==i) {
			print_u64( inp0, sizeof(inp0)/sizeof(uint64_t) );
			printf("->\n");
			print_u64( inp1 , sizeof(inp1)/sizeof(uint64_t) );
		}

		if( !check_eq( (uint8_t*)inp0 , (uint8_t*)inp2 , sizeof(inp0) ) ) {
			printf("neq: %d.\n", i );
			print_u8( (uint8_t *)inp0 , sizeof(inp0) );
			printf("->\n");
			print_u8( (uint8_t *)inp1 , sizeof(inp1) );
			printf("<-\n");
			print_u8( (uint8_t *)inp2 , sizeof(inp2) );
			eq = 0;
			break;
		}

	}

	printf("\ntest %s.\n\n", (eq)?"PASS":"FAIL");

#ifdef TESTSPEED
char bmmsg[256];
bm_dump(bmmsg,sizeof(bmmsg),&bm0);
printf("benchmark:\n%s\n", bmmsg );
#endif

	return 0;
}

