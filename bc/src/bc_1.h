#ifndef _BC_1_H_
#define _BC_1_H_

#include "stdint.h"


//
// libaray for basis conversion
// computation unit: 1 bit
//


void bc_1_256( void *poly , unsigned n_256bit );

void ibc_1_256( void *poly , unsigned n_256bit );


/////////////////////////////////////////

// n_byte >= 32
void bc_1( void * poly , unsigned n_byte );

void ibc_1( void * poly , unsigned n_byte );



/////////////////////////////////////////

void bc_1_16384( uint32_t *poly );

void ibc_1_16384( uint32_t *poly );


void bc_1_32768( uint32_t *poly );

void ibc_1_32768( uint32_t *poly );


void bc_1_65536( uint32_t *poly );

void ibc_1_65536( uint32_t *poly );






#endif
