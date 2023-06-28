
#include "gf264.h"


// extended GCD
// algorithm from https://doi.org/10.13154/tches.v2019.i3.340-398
uint64_t gf264_inv(uint64_t a)
{
    uint64_t f = 0x1b;
    uint64_t g = a << 1;
    uint64_t g0 = a >> 63;
    int64_t delta = 1;

    uint64_t r = 0;
    uint64_t v = 0;

    for (int i = 0; i < 127; i++)
    {
        uint64_t minus_delta = -delta;
        uint64_t swap = (minus_delta >> 63) & g0; // get sign bit
        // f0 is always 1

        // update delta
        delta = ((-swap) & (minus_delta << 1)) + delta + 1;

        // update f, g, v, r
        uint64_t vr = v ^ r;
        r ^= v & (-g0);
        v ^= vr & (-swap);
        // v0  = swap;

        uint64_t fg = f ^ g;
        g ^= f & (-g0);
        f ^= fg & (-swap);

        v = (v >> 1) | (swap << 63);
        g0 = g >> 63;
        g <<= 1;
    }
    return v;
}


