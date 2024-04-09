
#include "bc_512.h"



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
      poly[(i-deg_diff)*16]   ^= poly[i*16];
      poly[(i-deg_diff)*16+1] ^= poly[i*16+1];
      poly[(i-deg_diff)*16+2] ^= poly[i*16+2];
      poly[(i-deg_diff)*16+3] ^= poly[i*16+3];
      poly[(i-deg_diff)*16+4] ^= poly[i*16+4];
      poly[(i-deg_diff)*16+5] ^= poly[i*16+5];
      poly[(i-deg_diff)*16+6] ^= poly[i*16+6];
      poly[(i-deg_diff)*16+7] ^= poly[i*16+7];
      poly[(i-deg_diff)*16+8] ^= poly[i*16+8];
      poly[(i-deg_diff)*16+9] ^= poly[i*16+9];
      poly[(i-deg_diff)*16+10] ^= poly[i*16+10];
      poly[(i-deg_diff)*16+11] ^= poly[i*16+11];
      poly[(i-deg_diff)*16+12] ^= poly[i*16+12];
      poly[(i-deg_diff)*16+13] ^= poly[i*16+13];
      poly[(i-deg_diff)*16+14] ^= poly[i*16+14];
      poly[(i-deg_diff)*16+15] ^= poly[i*16+15];
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

void bc_512( void *poly , unsigned n_512 ) {  cvt( (uint32_t *)poly , n_512 , 0 , __builtin_ctz(n_512)  );}

////////

static inline
void idiv_blk( uint32_t *poly, int si_h, int si_l, int polylen )
{
  int deg_diff = si_h-si_l;
  for(int i=si_h;i<polylen;i++) {
      poly[(i-deg_diff)*16]   ^= poly[i*16];
      poly[(i-deg_diff)*16+1] ^= poly[i*16+1];
      poly[(i-deg_diff)*16+2] ^= poly[i*16+2];
      poly[(i-deg_diff)*16+3] ^= poly[i*16+3];
      poly[(i-deg_diff)*16+4] ^= poly[i*16+4];
      poly[(i-deg_diff)*16+5] ^= poly[i*16+5];
      poly[(i-deg_diff)*16+6] ^= poly[i*16+6];
      poly[(i-deg_diff)*16+7] ^= poly[i*16+7];
      poly[(i-deg_diff)*16+8] ^= poly[i*16+8];
      poly[(i-deg_diff)*16+9] ^= poly[i*16+9];
      poly[(i-deg_diff)*16+10] ^= poly[i*16+10];
      poly[(i-deg_diff)*16+11] ^= poly[i*16+11];
      poly[(i-deg_diff)*16+12] ^= poly[i*16+12];
      poly[(i-deg_diff)*16+13] ^= poly[i*16+13];
      poly[(i-deg_diff)*16+14] ^= poly[i*16+14];
      poly[(i-deg_diff)*16+15] ^= poly[i*16+15];
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

void ibc_512( void *poly , unsigned n_512 ) { icvt( (uint32_t *)poly , n_512 , 0 , __builtin_ctz(n_512)  ); }

