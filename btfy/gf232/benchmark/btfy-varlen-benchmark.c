
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>


static inline
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

static inline
void print_u16(const uint16_t *data, unsigned len )
{
	for(unsigned i=0;i<len;i++) {
		if( 0 == (i&7) ) printf("%3d: ",i);
		printf("0x%04x,", data[i] );
		if( 3 == (i&3) ) printf(" ");
		if( 7 == (i&7) ) printf("\n");
	}
}

static inline
void print_u32(const uint32_t *d32, unsigned l32 )
{
	for(unsigned i=0;i<l32;i++) {
		if( 0 == (i&7) ) printf("%3d: ",i);
		printf("0x%08x,", d32[i] );
		if( 3 == (i&3) ) printf(" ");
		if( 7 == (i&7) ) printf("\n");
	}
}

static inline
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


static inline
uint8_t check_eq( const uint8_t *vec0, const uint8_t *vec1, unsigned len)
{
	uint8_t diff = 0;
	for(unsigned i=0;i<len;i++){
		diff |= vec0[i]^vec1[i];
	}
	return diff==0;
}



#include "btfy32.h"
#include "randombytes.h"

#define TEST_RUN 100


#define LOG_LEN  16
#define LEN    (1<<LOG_LEN)



#define TESTSPEED

#ifdef TESTSPEED
#include "benchmark.h"
#endif



void test( uint32_t *inp0 , uint32_t *inp1 , uint32_t *inp2 , unsigned loglen )
{


	uint8_t eq = 1;

#ifdef TESTSPEED
struct benchmark bm0, bm1;
bm_init(&bm0);
bm_init(&bm1);
#endif

	unsigned len = (1<<loglen);
	printf("\nbtfy( %d x u32 = [%d]byte ).\n", len , len*4 );

	for(int i=1;i<=TEST_RUN;i++) {
		randombytes( (uint8_t*)inp0 , len*4 );
		if(0==i) {
			for(unsigned j=0;j<len;j++) inp0[j]=0xffffffff;
		}
		memcpy( inp1 , inp0 , len*4 );

#ifdef TESTSPEED
BENCHMARK( bm0, {
		btfy_32(inp1,loglen,0);
} );
#else
		btfy_32(inp1,loglen,0);
#endif

		memcpy( inp2 , inp1 , len*4 );

#ifdef TESTSPEED
BENCHMARK( bm1, {
		ibtfy_32(inp2,loglen,0);
} );
#else
		ibtfy_32(inp2,loglen,0);
#endif

		if(0==i) {
			print_u32( inp0, len );
			printf("->\n");
			print_u32( inp1 , len );
		}

		if( !check_eq( (uint8_t*)inp0 , (uint8_t*)inp2 , len*4 ) ) {
			printf("neq: %d.\n", i );
			print_u32( inp0 , len );
			printf("->\n");
			print_u32( inp1 , len );
			printf("<-\n");
			print_u32( inp2 , len );
			eq = 0;
			break;
		}

	}

	printf("\n%s for testing btfy([%d]xu64).\n", (eq)?"PASS":"FAIL" , len );

#ifdef TESTSPEED
char bmmsg[256];
bm_dump(bmmsg,sizeof(bmmsg),&bm0);
printf("bm0 :\n%s\n\n", bmmsg );
bm_dump(bmmsg,sizeof(bmmsg),&bm1);
printf("bm1 :\n%s\n\n", bmmsg );
#endif


}



int main(void)
{

	uint32_t inp0[LEN] __attribute__ ((aligned (32)));
	uint32_t inp1[LEN] __attribute__ ((aligned (32)));
	uint32_t inp2[LEN] __attribute__ ((aligned (32)));

	for(unsigned i=8;i<=LOG_LEN;i++) {
		test(inp0,inp1,inp2,i);
	}

	return 0;
}

