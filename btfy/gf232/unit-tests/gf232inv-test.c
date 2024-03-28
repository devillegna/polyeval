
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



void print_u32(const uint32_t *d32, unsigned l32 )
{
	for(unsigned i=0;i<l32;i++) {
		if( 0 == (i&7) ) printf("%3d: ",i);
		printf("0x%08x,", d32[i] );
		if( 3 == (i&3) ) printf(" ");
		if( 7 == (i&7) ) printf("\n");
	}
}




#include "gf232.h"
#include "randombytes.h"


#define TEST_RUN 1000000




int main(void)
{

	uint8_t pass = 1;

	printf("r x gf232_inv( r ) == 1 ?? \n");

	for(int i=0;i<TEST_RUN;i++) {
		uint32_t r;
		randombytes( (uint8_t*)&r , sizeof(r) );

		uint32_t r_inv = gf232_inv( r );

		uint32_t c = gf232_mul( r , r_inv );

		if( 1 != c ) {
			printf("!!!(%d).\n", i );
			printf("r: %x , r_inv: %x -x-> %x\n", r , r_inv , c );
			if( 0 == c ) continue;
			pass = 0;
			break;
		}

	}

	printf("\n%d test %s.\n\n", TEST_RUN , (pass)?"PASS":"FAIL");

	return 0;
}

