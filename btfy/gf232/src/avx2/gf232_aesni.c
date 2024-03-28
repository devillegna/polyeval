
#include "gf232.h"
#include "gf232_aesni.h"


uint32_t gf232_mul( uint32_t a , uint32_t b ) { return gf232_mul_u32_sse( a, b ); }

