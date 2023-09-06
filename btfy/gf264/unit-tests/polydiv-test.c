
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>


void print_u8(const unsigned char *data, unsigned len )
{
	for(unsigned i=0;i<len;i++) {
		if( 0 == (i&15) ) printf("%3d: ",i);
		printf("%x,", data[i] );
		if( 3 == (i&3) ) printf(".");
		if( 7 == (i&7) ) printf(" ");
		if( 15 == (i&15) ) printf("\n");
	}
}

void print_u64(const uint64_t *data, unsigned len )
{
	for(unsigned i=0;i<len;i++) {
		if( 0 == (i&15) ) printf("%3d: ",i);
		printf("%x,", data[i] );
		if( 3 == (i&3) ) printf(".");
		if( 7 == (i&7) ) printf(" ");
		if( 15 == (i&15) ) printf("\n");
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



#include "polydiv.h"


#define MAX_SI 1


#define LOG_LEN  16
#define LEN    (1<<LOG_LEN)




int main(void)
{

	uint64_t inp0[LEN] __attribute__ ((aligned (32)));
	uint64_t inp1[LEN] __attribute__ ((aligned (32)));
	uint64_t inp2[LEN] __attribute__ ((aligned (32)));

	uint8_t eq = 1;

	for(int l=2;l<LOG_LEN;l++) {
		unsigned len = 1<<l;
		printf("polydiv( %d x u64 ).\n", len );

		for(int i=0;i<=MAX_SI;i++) {
			printf(" / s_%d\n", i );
			for(int k=0;k<len;k++) inp0[k] = 1;
			memcpy( inp1 , inp0 , sizeof(uint64_t)*len );

			polydiv( inp1 , len , i );

			memcpy( inp2 , inp1 , sizeof(uint64_t)*len );

			ipolydiv( inp2 , len , i );

			if( len <= 256 ) {
				print_u64( inp0 , len );
				printf("->\n");
				print_u64( inp1 , len );
				printf("<-\n");
				print_u64( inp2 , len );
				printf("\n");
			}

			if( !check_eq( (uint8_t*)inp0 , (uint8_t*)inp2 , sizeof(inp0) ) ) {
				printf("neq: [%d] / s_%d.\n", len , i );

				print_u64( inp0 , len );
				printf("->\n");
				print_u64( inp1 , len );
				printf("<-\n");
				print_u64( inp2 , len );
				printf("\n");

				eq = 0;
				break;
			}
		}
		if( 0==eq ) break;
	}

	printf("\n%s.\n\n", (eq)?"ALL PASS":"FAIL");

	return 0;
}

