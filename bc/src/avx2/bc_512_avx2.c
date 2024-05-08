
#include "bc_512.h"


#include <emmintrin.h>
#include <immintrin.h>




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
void div_blk( __m256i *poly, int si_h, int si_l, int polylen )
{
  int deg_diff = si_h-si_l;
  for(int i=polylen-1;i>=si_h;i--) {
#if 1
      __m256i d0 = _mm256_loadu_si256(poly+(i-deg_diff)*2);
      __m256i d1 = _mm256_loadu_si256(poly+(i-deg_diff)*2+1);
      __m256i p0 = _mm256_loadu_si256(poly+i*2);
      __m256i p1 = _mm256_loadu_si256(poly+i*2+1);
      d0 ^= p0;
      d1 ^= p1;
      _mm256_storeu_si256( poly+(i-deg_diff)*2 , d0 );
      _mm256_storeu_si256( poly+(i-deg_diff)*2+1 , d1 );
      //poly[(i-deg_diff)*2]   ^= poly[i*2];
      //poly[(i-deg_diff)*2+1] ^= poly[i*2+1];
#else
      poly[(i-deg_diff)*2]   ^= poly[i*2];
      poly[(i-deg_diff)*2+1] ^= poly[i*2+1];
#endif
  }
}

static
void rep_in_si( __m256i *data, int datalen, int logsize_blk, int polyloglen_blk,  int si  )
{
  for(int i=polyloglen_blk-1;i>=si;i--) {
    int polylen = (1<<(i+logsize_blk+1));
    int si_h = (1<<(i+logsize_blk));
    int si_l = (1<<(i+logsize_blk-si));
    for(int j=0;j<datalen;j+=polylen) div_blk( data+j*2 , si_h, si_l, polylen );
  }
}

static
void cvt( __m256i *data, int datalen, int logsize_blk, int polyloglen_blk )
{
  if( 1 >= polyloglen_blk ) return;
  int si = choose_si(polyloglen_blk);
  rep_in_si( data, datalen , logsize_blk , polyloglen_blk , si );

  cvt( data , datalen , logsize_blk , si );
  cvt( data , datalen , logsize_blk+si , polyloglen_blk-si );
}




void bc_512( void * poly, unsigned n_512 ) {  cvt( (__m256i*) poly , n_512 , 0 , __builtin_ctz(n_512)  ); }


////////

static inline
void idiv_blk( __m256i *poly, int si_h, int si_l, int polylen )
{
  int deg_diff = si_h-si_l;
  for(int i=si_h;i<polylen;i++) {
#if 1
      __m256i d0 = _mm256_loadu_si256(poly+(i-deg_diff)*2);
      __m256i d1 = _mm256_loadu_si256(poly+(i-deg_diff)*2+1);
      __m256i p0 = _mm256_loadu_si256(poly+i*2);
      __m256i p1 = _mm256_loadu_si256(poly+i*2+1);
      d0 ^= p0;
      d1 ^= p1;
      _mm256_storeu_si256( poly+(i-deg_diff)*2 , d0 );
      _mm256_storeu_si256( poly+(i-deg_diff)*2+1 , d1 );
      //poly[(i-deg_diff)*2]   ^= poly[i*2];
      //poly[(i-deg_diff)*2+1] ^= poly[i*2+1];
#else
      poly[(i-deg_diff)*2]   ^= poly[i*2];
      poly[(i-deg_diff)*2+1] ^= poly[i*2+1];
#endif
  }
}

static
void irep_in_si( __m256i *data, int datalen, int logsize_blk, int polyloglen_blk,  int si  )
{
  for(int i=si;i<polyloglen_blk;i++) {
    int polylen = (1<<(i+logsize_blk+1));
    int si_h = (1<<(i+logsize_blk));
    int si_l = (1<<(i+logsize_blk-si));
    for(int j=0;j<datalen;j+=polylen) idiv_blk( data+j*2 , si_h, si_l, polylen );
  }
}

static
void icvt( __m256i *data, int datalen, int logsize_blk, int polyloglen_blk )
{
  if( 1 >= polyloglen_blk ) return;
  int si = choose_si(polyloglen_blk);

  icvt( data , datalen , logsize_blk+si , polyloglen_blk-si );
  icvt( data , datalen , logsize_blk , si );

  irep_in_si( data, datalen , logsize_blk , polyloglen_blk , si );
}



void ibc_512( void * poly, unsigned n_512 ) { icvt( (__m256i*) poly , n_512 , 0 , __builtin_ctz(n_512) ); }



