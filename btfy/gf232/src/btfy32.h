
#ifndef _BTFY32_H_
#define _BTFY32_H_


#include <stdint.h>


#ifdef  __cplusplus
extern  "C" {
#endif



void btfy_32( uint32_t * poly , unsigned log_polylen , uint32_t idx_offset );

void ibtfy_32( uint32_t * poly , unsigned log_polylen , uint32_t idx_offset );



#ifdef  __cplusplus
}
#endif


#endif
