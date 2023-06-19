
#include "gf264.h"
#include "gf264_aesni.h"

uint64_t gf264_mul( uint64_t a , uint64_t b ) { return gf2ext64_mul_u64( a, b ); }

