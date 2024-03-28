
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
#include "utils_hash.h"


#define TEST_RUN 10


#define LOG_LEN  16
#define LEN    (1<<LOG_LEN)




int main(void)
{

	uint32_t inp0[LEN] __attribute__ ((aligned (32)));
	uint32_t inp1[LEN] __attribute__ ((aligned (32)));
	uint32_t inp2[LEN] __attribute__ ((aligned (32)));

	uint8_t eq = 1;

	printf("btfy( %d x u32 ).\n", LEN );

	for(int i=0;i<TEST_RUN;i++) {
		randombytes( (uint8_t*)inp0 , sizeof(inp0) );
		memcpy( inp1 , inp0 , sizeof(inp0) );

		btfy_32(inp1,LOG_LEN,0);

		memcpy( inp2 , inp1 , sizeof(inp1) );

		ibtfy_32(inp2,LOG_LEN,0);

		if(0==i) {
			uint32_t hh[8];
			hash_ctx h0;
			if( hash_init(&h0) ) { printf("hash_init() 1 fail.\n"); continue; }
			hash_update(&h0,(uint8_t*)inp0,sizeof(inp0));
			hash_final_digest((uint8_t*)hh,32,&h0);
			print_u32( hh , 4 );
			printf("->\n");
			if( hash_init(&h0) ) { printf("hash_init() 2 fail.\n"); continue; }
			hash_update(&h0,(uint8_t*)inp1,sizeof(inp1));
			hash_final_digest((uint8_t*)hh,32,&h0);
			print_u32( hh , 4 );
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

	printf("\n%d test %s.\n\n", TEST_RUN , (eq)?"PASS":"FAIL");

	return 0;
}

