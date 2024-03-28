
#include "gf232.h"


// extended GCD
// algorithm from https://doi.org/10.13154/tches.v2019.i3.340-398
uint32_t gf232_inv(uint32_t a)
{
    uint32_t f = 0x8d;
    uint32_t g = a << 1;
    uint32_t g0 = a >> 31;
    int delta = 1;

    uint32_t r = 0;
    uint32_t v = 0;

    for (int i = 0; i < 63; i++)
    {
        uint32_t minus_delta = -delta;
        uint32_t swap = (minus_delta >> 31) & g0; // get sign bit
        // f0 is always 1

        // update delta
        delta = ((-swap) & (minus_delta << 1)) + delta + 1;

        // update f, g, v, r
        uint32_t vr = v ^ r;
        r ^= v & (-g0);
        v ^= vr & (-swap);
        // v0  = swap;

        uint32_t fg = f ^ g;
        g ^= f & (-g0);
        f ^= fg & (-swap);

        v = (v >> 1) | (swap << 31);
        g0 = g >> 31;
        g <<= 1;
    }
    return v;
}


