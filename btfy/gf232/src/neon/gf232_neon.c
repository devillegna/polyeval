
#include "gf232.h"
#include "gf232_neon.h"

uint32_t gf232_mul(uint32_t a , uint32_t b ) { return _gf232_mulx1_neon( a, b ); }
