
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


#include "dencoder.h"

#define LEN64   (16)

int main(void)
{
	uint64_t inp[LEN64];
	uint64_t out0[LEN64*2];
	uint64_t out1[LEN64*2];


	for(int i=0;i<5;i++) {
        for(int j=0;j<LEN64;j++) inp[j] = 794*(i+7)*(j+11);

        bc_encode_64( out0 , inp , LEN64 );

        ibc_decode_64( out1 , out0 , LEN64*2 );


        int r = check_eq( inp , out1 , LEN64*8 );
		printf("[%d] low:  %d\n", i , r );
		for(int j=0;j<LEN64;j++) {
			if (0 != out1[LEN64+j]) { r=0; break; }
		}
		printf("[%d] high: %d\n", i , r );
	}

	return 0;
}

