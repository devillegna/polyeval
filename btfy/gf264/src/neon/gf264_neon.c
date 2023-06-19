#include "gf264.h"
#include "gf264_neon.h"

uint64_t gf264_mul( uint64_t a , uint64_t b ) { return _gf264_mulx1_neon( a, b ); }

