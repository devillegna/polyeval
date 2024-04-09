

#include "bc_512.h"

#include <arm_neon.h>

static inline
int choose_si( int loglen )
{
  int si = 1<<0;
  for( int i=1; (1<<i) < loglen ; i++ ) {
    si = 1<<i;
  }
  return si;
}

/////////////////

static inline
void div_blk( uint32_t *poly, int si_h, int si_l, int polylen )
{
  int deg_diff = si_h-si_l;
  for(int i=polylen-1;i>=si_h;i--) {
    uint32x4_t src0 = vld1q_u32( poly+i*16 );
    uint32x4_t src1 = vld1q_u32( poly+i*16+4 );
    uint32x4_t src2 = vld1q_u32( poly+i*16+8 );
    uint32x4_t src3 = vld1q_u32( poly+i*16+12 );
    uint32x4_t dest0 = vld1q_u32( poly+(i-deg_diff)*16 );
    uint32x4_t dest1 = vld1q_u32( poly+(i-deg_diff)*16+4 );
    uint32x4_t dest2 = vld1q_u32( poly+(i-deg_diff)*16+8 );
    uint32x4_t dest3 = vld1q_u32( poly+(i-deg_diff)*16+12 );
    dest0 ^= src0;
    dest1 ^= src1;
    dest2 ^= src2;
    dest3 ^= src3;
    vst1q_u32( poly+(i-deg_diff)*16   , dest0 );
    vst1q_u32( poly+(i-deg_diff)*16+4 , dest1 );
    vst1q_u32( poly+(i-deg_diff)*16+8 , dest2 );
    vst1q_u32( poly+(i-deg_diff)*16+12 , dest3 );
  }
}

static
void rep_in_si( uint32_t *data, int datalen, int logsize_blk, int polyloglen_blk,  int si  )
{
  for(int i=polyloglen_blk-1;i>=si;i--) {
    int polylen = (1<<(i+logsize_blk+1));
    int si_h = (1<<(i+logsize_blk));
    int si_l = (1<<(i+logsize_blk-si));
    for(int j=0;j<datalen;j+=polylen) div_blk( data+j*16 , si_h, si_l, polylen );
  }
}

static
void cvt( uint32_t *data, int datalen, int logsize_blk, int polyloglen_blk )
{
  if( 1 >= polyloglen_blk ) return;
  int si = choose_si(polyloglen_blk);
  rep_in_si( data, datalen , logsize_blk , polyloglen_blk , si );

  cvt( data , datalen , logsize_blk , si );
  cvt( data , datalen , logsize_blk+si , polyloglen_blk-si );
}

void bc_512( void * poly, unsigned n_512 ) {  cvt( (uint32_t*)poly , n_512 , 0 , __builtin_ctz(n_512) ); }

////////

static inline
void idiv_blk( uint32_t *poly, int si_h, int si_l, int polylen )
{
  int deg_diff = si_h-si_l;
  for(int i=si_h;i<polylen;i++) {
    uint32x4_t src0 = vld1q_u32( poly+i*16 );
    uint32x4_t src1 = vld1q_u32( poly+i*16+4 );
    uint32x4_t src2 = vld1q_u32( poly+i*16+8 );
    uint32x4_t src3 = vld1q_u32( poly+i*16+12 );
    uint32x4_t dest0 = vld1q_u32( poly+(i-deg_diff)*16 );
    uint32x4_t dest1 = vld1q_u32( poly+(i-deg_diff)*16+4 );
    uint32x4_t dest2 = vld1q_u32( poly+(i-deg_diff)*16+8 );
    uint32x4_t dest3 = vld1q_u32( poly+(i-deg_diff)*16+12 );
    dest0 ^= src0;
    dest1 ^= src1;
    dest2 ^= src2;
    dest3 ^= src3;
    vst1q_u32( poly+(i-deg_diff)*16   , dest0 );
    vst1q_u32( poly+(i-deg_diff)*16+4 , dest1 );
    vst1q_u32( poly+(i-deg_diff)*16+8 , dest2 );
    vst1q_u32( poly+(i-deg_diff)*16+12 , dest3 );
  }
}

static
void irep_in_si( uint32_t *data, int datalen, int logsize_blk, int polyloglen_blk,  int si  )
{
  for(int i=si;i<polyloglen_blk;i++) {
    int polylen = (1<<(i+logsize_blk+1));
    int si_h = (1<<(i+logsize_blk));
    int si_l = (1<<(i+logsize_blk-si));
    for(int j=0;j<datalen;j+=polylen) idiv_blk( data+j*16 , si_h, si_l, polylen );
  }
}

static
void icvt( uint32_t *data, int datalen, int logsize_blk, int polyloglen_blk )
{
  if( 1 >= polyloglen_blk ) return;
  int si = choose_si(polyloglen_blk);

  icvt( data , datalen , logsize_blk+si , polyloglen_blk-si );
  icvt( data , datalen , logsize_blk , si );

  irep_in_si( data, datalen , logsize_blk , polyloglen_blk , si );
}

void ibc_512( void * poly, unsigned n_512 ) { icvt( (uint32_t*)poly , n_512 , 0 , __builtin_ctz(n_512) ); }

