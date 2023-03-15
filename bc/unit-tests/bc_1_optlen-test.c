
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

#define TEST_RUN 100

//#define LEN_BIT 16384
//#define LEN_BIT 32768
#define LEN_BIT 65536
#define LEN (LEN_BIT/32)


#if 16384 == LEN_BIT
#define LOG_LEN_BIT 14
#elif 32768 == LEN_BIT
#define LOG_LEN_BIT 15
#elif 65536 == LEN_BIT
#define LOG_LEN_BIT 16
#else
error
#endif

#if 0
#include "bc_ref.h"

void bc_1_ref( uint8_t * inp )
{
  uint8_t poly[LEN_BIT];
  for(int i=0;i<LEN_BIT;i++) {
    poly[i] = (inp[i/8]>>(i&7))&1;
  }

  bc_to_fft_ref( poly , 1 , LOG_LEN_BIT );

  for(int i=0;i<(LEN_BIT/8);i++) inp[i] = 0;
  for(int i=0;i<LEN_BIT;i++) {
    inp[i/8] |= poly[i]<<(i&7);
  }
}
#endif




#define TESTSPEED

#ifdef TESTSPEED
#include "benchmark.h"
#endif


int main()
{

	uint32_t inp0[LEN] __attribute__ ((aligned (32)));
	uint32_t inp1[LEN] __attribute__ ((aligned (32)));
	uint32_t inp2[LEN] __attribute__ ((aligned (32)));

	uint32_t inp4[LEN] __attribute__ ((aligned (32)));
	uint8_t * inp8 = (uint8_t*)inp4;


	uint8_t eq = 1;

#ifdef TESTSPEED
struct benchmark bm0;
bm_init(&bm0);
#endif

	printf("bc_1_%d.\n", LEN_BIT );

	for(int i=1;i<=TEST_RUN;i++) {
		randombytes( (uint8_t*)inp0 , sizeof(inp0) );
		if(0==i) {
			for(unsigned j=0;j<sizeof(inp0);j++) inp0[j]=0xffffffff;
		}
		memcpy( inp1 , inp0 , sizeof(inp0) );

#ifdef TESTSPEED
BENCHMARK( bm0, {
#endif
#if 16384 == LEN_BIT
		bc_1_16384( inp1 );
#elif 32768 == LEN_BIT
		bc_1_32768( inp1 );
#elif 65536 == LEN_BIT
		bc_1_65536( inp1 );
#else
error
#endif
		//cvt_to_fft_256_16384( inp1 );
		//repr_s8_1_16384( inp1 );
#ifdef TESTSPEED
} );
#endif

		memcpy( inp2 , inp1 , sizeof(inp1) );

#if 16384 == LEN_BIT
		ibc_1_16384( inp2 );
#elif 32768 == LEN_BIT
		ibc_1_32768( inp2 );
#elif 65536 == LEN_BIT
		ibc_1_65536( inp2 );
#else
error.
#endif
		//icvt_to_poly_256_16384( inp2 );
		//irepr_s8_1_16384( inp2 );
		if(0==i) {
			print_u32( inp0, sizeof(inp0)/sizeof(uint32_t) );
			printf("->\n");
			print_u32( inp1 , sizeof(inp1)/sizeof(uint32_t) );
		}

		if( !check_eq( (uint8_t*)inp0 , (uint8_t*)inp2 , sizeof(inp0) ) ) {
			printf("neq: %d.\n", i );
			print_u8( inp0 , sizeof(inp0) );
			printf("->\n");
			print_u8( inp1 , sizeof(inp1) );
			printf("<-\n");
			print_u8( inp2 , sizeof(inp2) );
			eq = 0;
			break;
		}
#if 0
		memcpy( inp8 , inp0 , sizeof(inp0) );
		bc_1_ref( inp8 );
		if( !check_eq( (uint8_t*)inp8 , (uint8_t*)inp1 , sizeof(inp1) ) ) {
			printf("neq: %d.\n", i );
			printf("bc_1_ref: %d.\n", i );
			print_u8( inp8 , sizeof(inp1) );
			printf("!=\n");
			print_u8( inp1 , sizeof(inp1) );
			eq = 0;
			break;
		}
#endif

	}

	printf("\ntest %s.\n\n", (eq)?"PASS":"FAIL");

#ifdef TESTSPEED
char bmmsg[256];
bm_dump(bmmsg,sizeof(bmmsg),&bm0);
printf("benchmark:\n%s\n", bmmsg );
#endif

	return 0;
}

