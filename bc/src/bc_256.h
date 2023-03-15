
#ifndef _BC_256_H_
#define _BC_256_H_

#include <stdint.h>

//
// libaray for basis conversion
// computation unit: 256 bits
//



void bc_256( void *poly , unsigned n_256 );

void ibc_256( void *poly , unsigned n_256 );

//////////


void bc_256_16384( uint32_t *poly );

void ibc_256_16384( uint32_t *poly );

void bc_256_32768( uint32_t *poly );

void ibc_256_32768( uint32_t *poly );

void bc_256_65536( uint32_t *poly );

void ibc_256_65536( uint32_t *poly );



#endif
