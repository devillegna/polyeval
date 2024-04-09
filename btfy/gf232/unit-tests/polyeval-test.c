
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


#include "gf232.h"
#include "cantor_to_gf232.h"
#include "btfy32.h"
#include "randombytes.h"


#define TEST_RUN 10


#define LEN    (32)
#define LOG_LEN  (5)


static uint32_t eval_LEN( const uint32_t * poly , uint32_t point )
{
  uint32_t s[LOG_LEN];
  for(int i=0;i<LOG_LEN;i++) { s[i] = cantor_to_gf232( point>>i ); }

  uint32_t r = poly[0];
  for(int i=1;i<LEN;i++) {
    uint32_t t = poly[i];
    for(int j=0;j<LOG_LEN;j++) { if((1<<j)&i) t = gf232_mul( t , s[j] ); }
    r ^= t;
  }
  return r;
}




int main(void)
{

	uint32_t inp0[LEN] __attribute__ ((aligned (32)));
	uint32_t btfy[LEN] __attribute__ ((aligned (32)));
	uint32_t eval[LEN] __attribute__ ((aligned (32)));

	uint8_t eq = 1;

	printf("btfy( %d x u32 ) = eval( %d ) ??\n", LEN , LEN );

	for(int i=0;i<TEST_RUN;i++) {
		randombytes( (uint8_t*)inp0 , sizeof(inp0) );

		memcpy( btfy , inp0 , sizeof(inp0) );
		btfy_32(btfy,LOG_LEN,0);

		for(unsigned j=0;j<LEN;j++) eval[j] = eval_LEN( inp0 , j );

		if( !check_eq( (uint8_t*)btfy , (uint8_t*)eval , sizeof(inp0) ) ) {
			printf("neq: %d.\n", i );
			print_u32( inp0 , LEN );
			printf("-btfy->\n");
			print_u32( btfy , LEN );
			printf("-eval->\n");
			print_u32( eval , LEN );
			eq = 0;
			break;
		}

	}

	printf("\n%d test %s.\n\n", TEST_RUN , (eq)?"PASS":"FAIL");

	return 0;
}

