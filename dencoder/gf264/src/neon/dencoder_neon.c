
#include "dencoder.h"

static uint64_t _enc_mat[32] = {
1,0xed57ce778f0d6244,0xb66864e6ec14b4d2,0x9cedebc39773a213,
0x54e5bf3774b3f850,0xaa41ec72aca4d601,0x4943d209fdc449b8,0x6e9a7c0e6d89dc84,
0xeb42e79f91f49f8c,0x5dc314ada2848ba2,0x503cdedab981d32c,0x21d93ab94919cfd3,
0x4553b41bda2fbec8,0xcb47a62d21df0b4d,0x4639691aa527dbba,0xe032e7a43c3d150f,
0xb0f502e4cd60039a,0xf8287767660ad2b4,0x1dbf8d7727fc654e,0x19b12a51dacba79e,
0xe691d7c0794c614d,0x28cf418faf97a29f,0xfc5c2bca45c931c,0x4cb433e14292028e,
0x189b29840c130e7f,0x315f35290c617710,0x1f4a9de7438026b,0x6aea5df783dd08ef,
0xb36ae8e530273635,0x837b54ac6006165d,0x49216693bb420e26,0xfd50e246ba0fe1a8,
};

static uint64_t _dec_mat[64] = {
1,0x8961b8a7841bfc6,0xfd56dc6306e45ec6,0x77332469d0404105,
0xd110eb4a970feaf5,0x9c94931de38226fe,0x52bf8c5ba0978896,0xbce2b823bebf0166,
0x7af3356b5d818e41,0x1e7b94c935b698d0,0xbcb73d9aabdd6ddf,0x5f75970481f366a8,
0xb33f7074d9b6e3fe,0x2adb64d29a93ef96,0x59dab92a69394e60,0x3c8004790328459c,
0xadcb0001cdaa3d05,0x9bb38fd76d5e0c0b,0xabcadd426cb54b4d,0xf042a8f93e8a48b9,
0x47b25aab146986f6,0xe750955be42931a1,0x9bf1d4a1535a0d17,0x91564f9a406d933c,
0xded3254b42d33aaa,0xbd1f32c59e84e45d,0xfa018660e3165396,0xa2f911a00f5baa7b,
0xe958695bc8926acd,0xc8c10606a8ec8d6b,0x8e6cfef54b2e9ac6,0x88007d210e507cf4,
0x17a89b9796a48397,0x1bdd636a99f218be,0x7cb821d83481a449,0xff39584957728d46,
0x9cd7c0f7b4d02de7,0x779fca946b676edb,0x5235a7728b715838,0x70ef0482c08984a5,
0xbb1c3dbfe2b5150b,0xee767356051fe0e1,0x1b5c98d75cc93fe8,0x4c25ee1fb1c8ceb6,
0xb7310079067f9e51,0x4be2a5ae359e2594,0xfe86864026631c3,0x5b4cd59958434683,
0xe550b2539a848707,0x6e50d13b6a7a5997,0x7fd8fb00a5aef99a,0xcdbd262ac16e806b,
0xc021427b4f7646e0,0xd7f262605a9bcdd0,0xb046630a7efeddc,0xa869f3dfcbe24fab,
0xa00150323be9b1de,0x1294d8724ae3d0d,0x4ab517f07f03a4b1,0xfdc5af3bb370e266,
0x2aff63efa8285238,0xfba0ace42227357a,0xc208107e1c3e5ac1,0xb2dbb5cb68d0dc74,
};



// sizeof output = 2 x sizeof input // assert(n_u64>=4)
static void encode_64_ref( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
{
    unsigned step = n_u64/4;
    for(unsigned i=0;i<step;i++){
        uint64_t rr[8] = {0};
        const uint8_t * ptr = ((const uint8_t *)fx)+i;
        for(int j=0;j<32;j++){
            uint8_t v = ptr[0];
            for(int k=0;k<8;k++) rr[k] ^= (-((v>>k)&1))&_enc_mat[j];
            ptr += step;
        }
        for(int k=0;k<8;k++){
            rfx[k] = rr[k];
        }
        rfx += 8;
    }
}


static void transpose_8x8( uint8_t * mat0 )
{
    uint8_t mat1[8];
    mat1[0] = (mat0[0]&0xf)|(mat0[4]<<4);
    mat1[1] = (mat0[1]&0xf)|(mat0[5]<<4);
    mat1[2] = (mat0[2]&0xf)|(mat0[6]<<4);
    mat1[3] = (mat0[3]&0xf)|(mat0[7]<<4);
    mat1[4] = (mat0[4]&0xf0)|(mat0[0]>>4);
    mat1[5] = (mat0[5]&0xf0)|(mat0[1]>>4);
    mat1[6] = (mat0[6]&0xf0)|(mat0[2]>>4);
    mat1[7] = (mat0[7]&0xf0)|(mat0[3]>>4);
    uint8_t mat2[8];
    mat2[0] = (mat1[0]&0x33)|((mat1[2]&0x33)<<2);
    mat2[1] = (mat1[1]&0x33)|((mat1[3]&0x33)<<2);
    mat2[2] = (mat1[2]&0xcc)|((mat1[0]&0xcc)>>2);
    mat2[3] = (mat1[3]&0xcc)|((mat1[1]&0xcc)>>2);
    mat2[4] = (mat1[4]&0x33)|((mat1[6]&0x33)<<2);
    mat2[5] = (mat1[5]&0x33)|((mat1[7]&0x33)<<2);
    mat2[6] = (mat1[6]&0xcc)|((mat1[4]&0xcc)>>2);
    mat2[7] = (mat1[7]&0xcc)|((mat1[5]&0xcc)>>2);

    mat0[0] = (mat2[0]&0x55)|((mat2[1]&0x55)<<1);
    mat0[1] = (mat2[1]&0xaa)|((mat2[0]&0xaa)>>1);
    mat0[2] = (mat2[2]&0x55)|((mat2[3]&0x55)<<1);
    mat0[3] = (mat2[3]&0xaa)|((mat2[2]&0xaa)>>1);
    mat0[4] = (mat2[4]&0x55)|((mat2[5]&0x55)<<1);
    mat0[5] = (mat2[5]&0xaa)|((mat2[4]&0xaa)>>1);
    mat0[6] = (mat2[6]&0x55)|((mat2[7]&0x55)<<1);
    mat0[7] = (mat2[7]&0xaa)|((mat2[6]&0xaa)>>1);

}


// sizeof output = sizeof input // assert(n_u64>=8)
static void decode_64_ref( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
{
    unsigned step = n_u64/8;
    for(unsigned i=0;i<step;i++){
        const uint64_t * inp = fx+i*8;
        uint64_t rr[8] = {0};
        for(int j=0;j<64;j++) {
            for(int k=0;k<8;k++) { rr[k] ^= (-((inp[k]>>j)&1))&_dec_mat[j]; }
        }

        // output
        uint8_t * ptr = ((uint8_t *)rfx)+i;
        for(int k=0;k<8;k++) {
            uint8_t mat0[8];
            uint8_t * p8 = ((uint8_t*)&rr[0])+k;
            for(int j=0;j<8;j++) { mat0[j] = p8[j*8]; }
            transpose_8x8(mat0);
            for(int j=0;j<8;j++){
                ptr[0] = mat0[j];
                ptr += step;
            }
        }
    }
}




#include "arm_neon.h"

static inline void _bit_transpose_8x8_x16( uint8x16_t * mat0 )
{
    uint8x16_t _0x0f = vdupq_n_u8(0xf);
    uint8x16_t _0xf0 = vdupq_n_u8(0xf0);
    uint8x16_t _0x33 = vdupq_n_u8(0x33);
    uint8x16_t _0xcc = vdupq_n_u8(0xcc);
    uint8x16_t _0x55 = vdupq_n_u8(0x55);
    uint8x16_t _0xaa = vdupq_n_u8(0xaa);

    uint8x16_t mat1[8];
    mat1[0] = (mat0[0]&_0x0f)|(mat0[4]<<4);
    mat1[1] = (mat0[1]&_0x0f)|(mat0[5]<<4);
    mat1[2] = (mat0[2]&_0x0f)|(mat0[6]<<4);
    mat1[3] = (mat0[3]&_0x0f)|(mat0[7]<<4);
    mat1[4] = (mat0[4]&_0xf0)|(mat0[0]>>4);
    mat1[5] = (mat0[5]&_0xf0)|(mat0[1]>>4);
    mat1[6] = (mat0[6]&_0xf0)|(mat0[2]>>4);
    mat1[7] = (mat0[7]&_0xf0)|(mat0[3]>>4);
    uint8x16_t mat2[8];
    mat2[0] = (mat1[0]&_0x33)|((mat1[2]&_0x33)<<2);
    mat2[1] = (mat1[1]&_0x33)|((mat1[3]&_0x33)<<2);
    mat2[2] = (mat1[2]&_0xcc)|((mat1[0]&_0xcc)>>2);
    mat2[3] = (mat1[3]&_0xcc)|((mat1[1]&_0xcc)>>2);
    mat2[4] = (mat1[4]&_0x33)|((mat1[6]&_0x33)<<2);
    mat2[5] = (mat1[5]&_0x33)|((mat1[7]&_0x33)<<2);
    mat2[6] = (mat1[6]&_0xcc)|((mat1[4]&_0xcc)>>2);
    mat2[7] = (mat1[7]&_0xcc)|((mat1[5]&_0xcc)>>2);

    mat0[0] = (mat2[0]&_0x55)|((mat2[1]&_0x55)<<1);
    mat0[1] = (mat2[1]&_0xaa)|((mat2[0]&_0xaa)>>1);
    mat0[2] = (mat2[2]&_0x55)|((mat2[3]&_0x55)<<1);
    mat0[3] = (mat2[3]&_0xaa)|((mat2[2]&_0xaa)>>1);
    mat0[4] = (mat2[4]&_0x55)|((mat2[5]&_0x55)<<1);
    mat0[5] = (mat2[5]&_0xaa)|((mat2[4]&_0xaa)>>1);
    mat0[6] = (mat2[6]&_0x55)|((mat2[7]&_0x55)<<1);
    mat0[7] = (mat2[7]&_0xaa)|((mat2[6]&_0xaa)>>1);

}


// lower 64-bits store even terms. higher 64s stroe odd terms.
static inline void _byte_transpose_8x8_x2( uint8x16_t * row )
{
    uint8x16x2_t t0 = vtrnq_u8(row[0], row[1]);
    uint8x16x2_t t1 = vtrnq_u8(row[2], row[3]);
    uint8x16x2_t t2 = vtrnq_u8(row[4], row[5]);
    uint8x16x2_t t3 = vtrnq_u8(row[6], row[7]);

    uint16x8x2_t t4 = vtrnq_u16(vreinterpretq_u16_u8(t0.val[0]), vreinterpretq_u16_u8(t1.val[0]));
    uint16x8x2_t t5 = vtrnq_u16(vreinterpretq_u16_u8(t0.val[1]), vreinterpretq_u16_u8(t1.val[1]));
    uint16x8x2_t t6 = vtrnq_u16(vreinterpretq_u16_u8(t2.val[0]), vreinterpretq_u16_u8(t3.val[0]));
    uint16x8x2_t t7 = vtrnq_u16(vreinterpretq_u16_u8(t2.val[1]), vreinterpretq_u16_u8(t3.val[1]));

    uint32x4x2_t t8  = vtrnq_u32(vreinterpretq_u32_u16(t4.val[0]), vreinterpretq_u32_u16(t6.val[0]));
    uint32x4x2_t t9  = vtrnq_u32(vreinterpretq_u32_u16(t4.val[1]), vreinterpretq_u32_u16(t6.val[1]));
    uint32x4x2_t t10 = vtrnq_u32(vreinterpretq_u32_u16(t5.val[0]), vreinterpretq_u32_u16(t7.val[0]));
    uint32x4x2_t t11 = vtrnq_u32(vreinterpretq_u32_u16(t5.val[1]), vreinterpretq_u32_u16(t7.val[1]));

    row[0] = vreinterpretq_u8_u32(  t8.val[0] );
    row[1] = vreinterpretq_u8_u32( t10.val[0] );
    row[2] = vreinterpretq_u8_u32(  t9.val[0] );
    row[3] = vreinterpretq_u8_u32( t11.val[0] );
    row[4] = vreinterpretq_u8_u32(  t8.val[1] );
    row[5] = vreinterpretq_u8_u32( t10.val[1] );
    row[6] = vreinterpretq_u8_u32(  t9.val[1] );
    row[7] = vreinterpretq_u8_u32( t11.val[1] );
}

static inline void _byte_transpose_4x4( uint8x16_t * row )
{
    uint8x16x2_t t0 = vtrnq_u8(row[0], row[1]);
    uint8x16x2_t t1 = vtrnq_u8(row[2], row[3]);

    uint16x8x2_t t2 = vtrnq_u16(vreinterpretq_u16_u8(t0.val[0]), vreinterpretq_u16_u8(t1.val[0]));
    uint16x8x2_t t3 = vtrnq_u16(vreinterpretq_u16_u8(t0.val[1]), vreinterpretq_u16_u8(t1.val[1]));

    row[0] = vreinterpretq_u8_u16( t2.val[0] );
    row[1] = vreinterpretq_u8_u16( t3.val[0] );
    row[2] = vreinterpretq_u8_u16( t2.val[1] );
    row[3] = vreinterpretq_u8_u16( t3.val[1] );
}





static void gather_tr_32_x32( uint32_t * dest , const uint8_t * src , unsigned step )
{
    uint32_t mat32[8*4];
    uint8_t * mat = (uint8_t*) mat32;
    uint8x16_t tmp0[8];
    for(int i=0;i<32;i++){
        uint8_t b0 = src[i*step+0];
        uint8_t b1 = src[i*step+1];
        uint8_t b2 = src[i*step+2];
        uint8_t b3 = src[i*step+3];
        mat[(i&7)*16+(i/8)+0]=b0;
        mat[(i&7)*16+(i/8)+8]=b1;
        mat[(i&7)*16+(i/8)+4]=b2;
        mat[(i&7)*16+(i/8)+12]=b3;
    }
    for(int i=0;i<8;i++) { tmp0[i] = vld1q_u8(mat+i*16); }
    _bit_transpose_8x8_x16(tmp0);
    uint32x4x4_t tmp1;
    uint32x4x4_t tmp2;
    tmp1.val[0] = vtrn1q_u32(tmp0[0],tmp0[4]);
    tmp1.val[1] = vtrn1q_u32(tmp0[1],tmp0[5]);
    tmp1.val[2] = vtrn1q_u32(tmp0[2],tmp0[6]);
    tmp1.val[3] = vtrn1q_u32(tmp0[3],tmp0[7]);
    tmp2.val[0] = vtrn2q_u32(tmp0[0],tmp0[4]);
    tmp2.val[1] = vtrn2q_u32(tmp0[1],tmp0[5]);
    tmp2.val[2] = vtrn2q_u32(tmp0[2],tmp0[6]);
    tmp2.val[3] = vtrn2q_u32(tmp0[3],tmp0[7]);
    vst4q_u32( dest , tmp1 );
    vst4q_u32( dest+16 , tmp2 );
}



extern const uint8_t enc_4bit_tbls[1024];

static void linearmap_32_x32( uint64_t * des , const uint32_t *src )
{
    uint8x16_t reorder_idx = {0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15};
    uint8x16_t m_0xf = vdupq_n_u8(0xf);
    uint8x16_t temp0[4],temp1[4];
    uint8x16_t dest0[8],dest1[8];

    for(int i=0;i<4;i++) temp0[i] = vld1q_u8( (const uint8_t*)(src + 4*i) );
    for(int i=0;i<4;i++) temp1[i] = vld1q_u8( (const uint8_t*)(src + 4*i + 16) );
    _byte_transpose_4x4(temp0);
    _byte_transpose_4x4(temp1);
    for(int i=0;i<4;i++) { // reorder to fit byte_transpose_8x8_x2()
        temp0[i] = vqtbl1q_u8(temp0[i],reorder_idx);
        temp1[i] = vqtbl1q_u8(temp1[i],reorder_idx);
    }
    uint8x16_t nib_00, nib_01, nib_10, nib_11;
    nib_00 = temp0[0]&m_0xf;
    nib_01 = temp0[0]>>4;
    nib_10 = temp1[0]&m_0xf;
    nib_11 = temp1[0]>>4;
    for(int k=0;k<8;k++) {
        uint8x16_t tab = vld1q_u8( enc_4bit_tbls+k*16 );
        dest0[k] = vqtbl1q_u8( tab , nib_00 );
        dest1[k] = vqtbl1q_u8( tab , nib_10 );
    }
    for(int k=0;k<8;k++) {
        uint8x16_t tab = vld1q_u8( enc_4bit_tbls+(8+k)*16 );
        dest0[k] ^= vqtbl1q_u8( tab , nib_01 );
        dest1[k] ^= vqtbl1q_u8( tab , nib_11 );
    }
    for(int l=1;l<4;l++) {
        nib_00 = temp0[l]&m_0xf;
        nib_01 = temp0[l]>>4;
        nib_10 = temp1[l]&m_0xf;
        nib_11 = temp1[l]>>4;
        for(int k=0;k<8;k++) {
            uint8x16_t tab = vld1q_u8( enc_4bit_tbls+(l*16+k)*16 );
            dest0[k] ^= vqtbl1q_u8( tab , nib_00 );
            dest1[k] ^= vqtbl1q_u8( tab , nib_10 );
        }
        for(int k=0;k<8;k++) {
            uint8x16_t tab = vld1q_u8( enc_4bit_tbls+(l*16+8+k)*16 );
            dest0[k] ^= vqtbl1q_u8( tab , nib_01 );
            dest1[k] ^= vqtbl1q_u8( tab , nib_11 );
        }
    }
    _byte_transpose_8x8_x2( dest0 );
    for(int i=0;i<8;i++) { vst1q_u64( des + i*2 , vreinterpretq_u64_u8(dest0[i]) ); }
    des += 16;
    _byte_transpose_8x8_x2( dest1 );
    for(int i=0;i<8;i++) { vst1q_u64( des + i*2 , vreinterpretq_u64_u8(dest1[i]) ); }
}



// sizeof output = 2 x sizeof input // assert(n_u64>=4)
void encode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
{
    if ( n_u64 < 16 ) { encode_64_ref(rfx,fx,n_u64); return; }
    unsigned step = n_u64/4;
    for(unsigned i=0;i<step;i+=4){
        uint32_t temp[32];
        gather_tr_32_x32( temp , ((const uint8_t *)fx)+i , step );
        linearmap_32_x32( rfx+i*8 , temp );
    }
}




extern const uint8_t dec_4bit_tbls[2048];

static void linearmap_64_x16( uint64_t * des , const uint64_t *src )
{
    uint8x16_t m_0xf = vdupq_n_u8(0xf);
    uint8x16_t temp[8];
    uint8x16_t dest[8];
    for(int j=0;j<8;j++) temp[j] = vld1q_u8( (const uint8_t*)(src + j*2) );
    _byte_transpose_8x8_x2( temp );
    uint8x16_t nib_0 = temp[0]&m_0xf;
    uint8x16_t nib_1 = vshrq_n_u8(temp[0],4);
    for(int k=0;k<8;k++) dest[k]  = vqtbl1q_u8( vld1q_u8( dec_4bit_tbls+k*16 ) , nib_0 );
    for(int k=0;k<8;k++) dest[k] ^= vqtbl1q_u8( vld1q_u8( dec_4bit_tbls+(8+k)*16 ) , nib_1 );
    for(int l=1;l<8;l++) {
        nib_0 = temp[l]&m_0xf;
        nib_1 = vshrq_n_u8(temp[l],4);
        for(int k=0;k<8;k++) dest[k] ^= vqtbl1q_u8( vld1q_u8( dec_4bit_tbls+(l*16+k)*16 ) , nib_0 );
        for(int k=0;k<8;k++) dest[k] ^= vqtbl1q_u8( vld1q_u8( dec_4bit_tbls+(l*16+8+k)*16 ) , nib_1 );
    }
    _byte_transpose_8x8_x2( dest );
    for(int j=0;j<8;j++) vst1q_u64( des + j*2 , vreinterpretq_u64_u8(dest[j]) );
}


static void tr_scatter_64_x16( uint8_t * rfx , const uint64_t * inp64 , unsigned step )
{
    uint8x16_t temp0[8];
    uint8x16_t temp1[8];
    for(int j=0;j<8;j++) temp0[j] = vld1q_u8( (const uint8_t*)(inp64 + j*2) );
    temp1[0] = vzip1q_u64(vreinterpretq_u64_u8(temp0[0]),vreinterpretq_u64_u8(temp0[4]));
    temp1[1] = vzip2q_u64(vreinterpretq_u64_u8(temp0[0]),vreinterpretq_u64_u8(temp0[4]));
    temp1[2] = vzip1q_u64(vreinterpretq_u64_u8(temp0[1]),vreinterpretq_u64_u8(temp0[5]));
    temp1[3] = vzip2q_u64(vreinterpretq_u64_u8(temp0[1]),vreinterpretq_u64_u8(temp0[5]));
    temp1[4] = vzip1q_u64(vreinterpretq_u64_u8(temp0[2]),vreinterpretq_u64_u8(temp0[6]));
    temp1[5] = vzip2q_u64(vreinterpretq_u64_u8(temp0[2]),vreinterpretq_u64_u8(temp0[6]));
    temp1[6] = vzip1q_u64(vreinterpretq_u64_u8(temp0[3]),vreinterpretq_u64_u8(temp0[7]));
    temp1[7] = vzip2q_u64(vreinterpretq_u64_u8(temp0[3]),vreinterpretq_u64_u8(temp0[7]));
    _bit_transpose_8x8_x16( temp1 );

    // it's a shame that there is no scatter intruction in neon.
    uint8_t m0[16];
    for(int j=0;j<8;j++){
        vst1q_u8( m0 , temp1[j] );
        for(int l=0;l<8;l++) {
            rfx[(j+l*8)*step] = m0[l];
            rfx[(j+l*8)*step+1] = m0[l+8];
        }
    }
}



void decode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
{
    if ( n_u64 < 16 ) { decode_64_ref(rfx,fx,n_u64); return; }

    uint64_t temp[16];
    for(unsigned i=0;i<n_u64;i+=16){
        const uint64_t * inp = fx+i;
        linearmap_64_x16( temp , inp );

        tr_scatter_64_x16( ((uint8_t*)rfx)+(i/8) , temp , n_u64/8 );
    }
}


static const uint8_t enc_4bit_tbls[8*8*16] = {
0x00,0x01,0x44,0x45,0xd2,0xd3,0x96,0x97,0x13,0x12,0x57,0x56,0xc1,0xc0,0x85,0x84,
0x00,0x00,0x62,0x62,0xb4,0xb4,0xd6,0xd6,0xa2,0xa2,0xc0,0xc0,0x16,0x16,0x74,0x74,
0x00,0x00,0x0d,0x0d,0x14,0x14,0x19,0x19,0x73,0x73,0x7e,0x7e,0x67,0x67,0x6a,0x6a,
0x00,0x00,0x8f,0x8f,0xec,0xec,0x63,0x63,0x97,0x97,0x18,0x18,0x7b,0x7b,0xf4,0xf4,
0x00,0x00,0x77,0x77,0xe6,0xe6,0x91,0x91,0xc3,0xc3,0xb4,0xb4,0x25,0x25,0x52,0x52,
0x00,0x00,0xce,0xce,0x64,0x64,0xaa,0xaa,0xeb,0xeb,0x25,0x25,0x8f,0x8f,0x41,0x41,
0x00,0x00,0x57,0x57,0x68,0x68,0x3f,0x3f,0xed,0xed,0xba,0xba,0x85,0x85,0xd2,0xd2,
0x00,0x00,0xed,0xed,0xb6,0xb6,0x5b,0x5b,0x9c,0x9c,0x71,0x71,0x2a,0x2a,0xc7,0xc7,
0x00,0x50,0x01,0x51,0xb8,0xe8,0xb9,0xe9,0x84,0xd4,0x85,0xd5,0x3c,0x6c,0x3d,0x6d,
0x00,0xf8,0xd6,0x2e,0x49,0xb1,0x9f,0x67,0xdc,0x24,0x0a,0xf2,0x95,0x6d,0x43,0xbb,
0x00,0xb3,0xa4,0x17,0xc4,0x77,0x60,0xd3,0x89,0x3a,0x2d,0x9e,0x4d,0xfe,0xe9,0x5a,
0x00,0x74,0xac,0xd8,0xfd,0x89,0x51,0x25,0x6d,0x19,0xc1,0xb5,0x90,0xe4,0x3c,0x48,
0x00,0x37,0x72,0x45,0x09,0x3e,0x7b,0x4c,0x0e,0x39,0x7c,0x4b,0x07,0x30,0x75,0x42,
0x00,0xbf,0xec,0x53,0xd2,0x6d,0x3e,0x81,0x7c,0xc3,0x90,0x2f,0xae,0x11,0x42,0xfd,
0x00,0xe5,0x41,0xa4,0x43,0xa6,0x02,0xe7,0x9a,0x7f,0xdb,0x3e,0xd9,0x3c,0x98,0x7d,
0x00,0x54,0xaa,0xfe,0x49,0x1d,0xe3,0xb7,0x6e,0x3a,0xc4,0x90,0x27,0x73,0x8d,0xd9,
0x00,0x8c,0xa2,0x2e,0x2c,0xa0,0x8e,0x02,0xd3,0x5f,0x71,0xfd,0xff,0x73,0x5d,0xd1,
0x00,0x9f,0x8b,0x14,0xd3,0x4c,0x58,0xc7,0xcf,0x50,0x44,0xdb,0x1c,0x83,0x97,0x08,
0x00,0xf4,0x84,0x70,0x81,0x75,0x05,0xf1,0x19,0xed,0x9d,0x69,0x98,0x6c,0x1c,0xe8,
0x00,0x91,0xa2,0x33,0xb9,0x28,0x1b,0x8a,0x49,0xd8,0xeb,0x7a,0xf0,0x61,0x52,0xc3,
0x00,0x9f,0xad,0x32,0xda,0x45,0x77,0xe8,0xb9,0x26,0x14,0x8b,0x63,0xfc,0xce,0x51,
0x00,0xe7,0x14,0xf3,0xde,0x39,0xca,0x2d,0x3a,0xdd,0x2e,0xc9,0xe4,0x03,0xf0,0x17,
0x00,0x42,0xc3,0x81,0x3c,0x7e,0xff,0xbd,0xd9,0x9b,0x1a,0x58,0xe5,0xa7,0x26,0x64,
0x00,0xeb,0x5d,0xb6,0x50,0xbb,0x0d,0xe6,0x21,0xca,0x7c,0x97,0x71,0x9a,0x2c,0xc7,
0x00,0xc8,0x4d,0x85,0xba,0x72,0xf7,0x3f,0x0f,0xc7,0x42,0x8a,0xb5,0x7d,0xf8,0x30,
0x00,0xbe,0x0b,0xb5,0xdb,0x65,0xd0,0x6e,0x15,0xab,0x1e,0xa0,0xce,0x70,0xc5,0x7b,
0x00,0x2f,0xdf,0xf0,0x27,0x08,0xf8,0xd7,0x3d,0x12,0xe2,0xcd,0x1a,0x35,0xc5,0xea,
0x00,0xda,0x21,0xfb,0xa5,0x7f,0x84,0x5e,0x3c,0xe6,0x1d,0xc7,0x99,0x43,0xb8,0x62,
0x00,0x1b,0x2d,0x36,0x1a,0x01,0x37,0x2c,0xa4,0xbf,0x89,0x92,0xbe,0xa5,0x93,0x88,
0x00,0xb4,0xa6,0x12,0x69,0xdd,0xcf,0x7b,0xe7,0x53,0x41,0xf5,0x8e,0x3a,0x28,0x9c,
0x00,0x53,0x47,0x14,0x39,0x6a,0x7e,0x2d,0x32,0x61,0x75,0x26,0x0b,0x58,0x4c,0x1f,
0x00,0x45,0xcb,0x8e,0x46,0x03,0x8d,0xc8,0xe0,0xa5,0x2b,0x6e,0xa6,0xe3,0x6d,0x28,
0x00,0x9a,0xb4,0x2e,0x4e,0xd4,0xfa,0x60,0x9e,0x04,0x2a,0xb0,0xd0,0x4a,0x64,0xfe,
0x00,0x03,0xd2,0xd1,0x65,0x66,0xb7,0xb4,0xa7,0xa4,0x75,0x76,0xc2,0xc1,0x10,0x13,
0x00,0x60,0x0a,0x6a,0xfc,0x9c,0xf6,0x96,0xcb,0xab,0xc1,0xa1,0x37,0x57,0x3d,0x5d,
0x00,0xcd,0x66,0xab,0x27,0xea,0x41,0x8c,0xda,0x17,0xbc,0x71,0xfd,0x30,0x9b,0x56,
0x00,0xe4,0x67,0x83,0x77,0x93,0x10,0xf4,0x51,0xb5,0x36,0xd2,0x26,0xc2,0x41,0xa5,
0x00,0x02,0x77,0x75,0x8d,0x8f,0xfa,0xf8,0x2a,0x28,0x5d,0x5f,0xa7,0xa5,0xd0,0xd2,
0x00,0xf5,0x28,0xdd,0xbf,0x4a,0x97,0x62,0xb1,0x44,0x99,0x6c,0x0e,0xfb,0x26,0xd3,
0x00,0xb0,0xf8,0x48,0x1d,0xad,0xe5,0x55,0x19,0xa9,0xe1,0x51,0x04,0xb4,0xfc,0x4c,
0x00,0x4d,0x9f,0xd2,0x1c,0x51,0x83,0xce,0x8e,0xc3,0x11,0x5c,0x92,0xdf,0x0d,0x40,
0x00,0x61,0xa2,0xc3,0x93,0xf2,0x31,0x50,0x02,0x63,0xa0,0xc1,0x91,0xf0,0x33,0x52,
0x00,0x4c,0x97,0xdb,0x5c,0x10,0xcb,0x87,0x92,0xde,0x05,0x49,0xce,0x82,0x59,0x15,
0x00,0x79,0xaf,0xd6,0xa4,0xdd,0x0b,0x72,0x42,0x3b,0xed,0x94,0xe6,0x9f,0x49,0x30,
0x00,0xc0,0x8f,0x4f,0xbc,0x7c,0x33,0xf3,0xe1,0x21,0x6e,0xae,0x5d,0x9d,0xd2,0x12,
0x00,0xd7,0x41,0x96,0xc2,0x15,0x83,0x54,0x33,0xe4,0x72,0xa5,0xf1,0x26,0xb0,0x67,
0x00,0x91,0xcf,0x5e,0xc5,0x54,0x0a,0x9b,0xb4,0x25,0x7b,0xea,0x71,0xe0,0xbe,0x2f,
0x00,0xe6,0x28,0xce,0x0f,0xe9,0x27,0xc1,0x4c,0xaa,0x64,0x82,0x43,0xa5,0x6b,0x8d,
0x00,0x7f,0x10,0x6f,0x6b,0x14,0x7b,0x04,0xef,0x90,0xff,0x80,0x84,0xfb,0x94,0xeb,
0x00,0x0e,0x77,0x79,0x02,0x0c,0x75,0x7b,0x08,0x06,0x7f,0x71,0x0a,0x04,0x7d,0x73,
0x00,0x13,0x61,0x72,0x38,0x2b,0x59,0x4a,0xdd,0xce,0xbc,0xaf,0xe5,0xf6,0x84,0x97,
0x00,0x0c,0x0c,0x00,0x74,0x78,0x78,0x74,0x83,0x8f,0x8f,0x83,0xf7,0xfb,0xfb,0xf7,
0x00,0x84,0x29,0xad,0xde,0x5a,0xf7,0x73,0xf7,0x73,0xde,0x5a,0x29,0xad,0x00,0x84,
0x00,0x29,0x35,0x1c,0xa9,0x80,0x9c,0xb5,0x5d,0x74,0x68,0x41,0xf4,0xdd,0xc1,0xe8,
0x00,0x9b,0x5f,0xc4,0xf4,0x6f,0xab,0x30,0xea,0x71,0xb5,0x2e,0x1e,0x85,0x41,0xda,
0x00,0x18,0x31,0x29,0x01,0x19,0x30,0x28,0x6a,0x72,0x5b,0x43,0x6b,0x73,0x5a,0x42,
0x00,0x35,0x5d,0x68,0x26,0x13,0x7b,0x4e,0xa8,0x9d,0xf5,0xc0,0x8e,0xbb,0xd3,0xe6,
0x00,0x36,0x16,0x20,0x0e,0x38,0x18,0x2e,0xe1,0xd7,0xf7,0xc1,0xef,0xd9,0xf9,0xcf,
0x00,0x27,0x06,0x21,0x42,0x65,0x44,0x63,0x0f,0x28,0x09,0x2e,0x4d,0x6a,0x4b,0x6c,
0x00,0x30,0x60,0x50,0xbb,0x8b,0xdb,0xeb,0xba,0x8a,0xda,0xea,0x01,0x31,0x61,0x51,
0x00,0xe5,0xac,0x49,0x93,0x76,0x3f,0xda,0x46,0xa3,0xea,0x0f,0xd5,0x30,0x79,0x9c,
0x00,0xe8,0x54,0xbc,0x66,0x8e,0x32,0xda,0xe2,0x0a,0xb6,0x5e,0x84,0x6c,0xd0,0x38,
0x00,0x6a,0x7b,0x11,0x21,0x4b,0x5a,0x30,0x50,0x3a,0x2b,0x41,0x71,0x1b,0x0a,0x60,
0x00,0xb3,0x83,0x30,0x49,0xfa,0xca,0x79,0xfd,0x4e,0x7e,0xcd,0xb4,0x07,0x37,0x84,
};


static const uint8_t dec_4bit_tbls[16*8*16] = {
0x00,0x01,0xc6,0xc7,0xc6,0xc7,0x00,0x01,0x05,0x04,0xc3,0xc2,0xc3,0xc2,0x05,0x04,
0x00,0x00,0xbf,0xbf,0x5e,0x5e,0xe1,0xe1,0x41,0x41,0xfe,0xfe,0x1f,0x1f,0xa0,0xa0,
0x00,0x00,0x41,0x41,0xe4,0xe4,0xa5,0xa5,0x40,0x40,0x01,0x01,0xa4,0xa4,0xe5,0xe5,
0x00,0x00,0x78,0x78,0x06,0x06,0x7e,0x7e,0xd0,0xd0,0xa8,0xa8,0xd6,0xd6,0xae,0xae,
0x00,0x00,0x8a,0x8a,0x63,0x63,0xe9,0xe9,0x69,0x69,0xe3,0xe3,0x0a,0x0a,0x80,0x80,
0x00,0x00,0x1b,0x1b,0xdc,0xdc,0xc7,0xc7,0x24,0x24,0x3f,0x3f,0xf8,0xf8,0xe3,0xe3,
0x00,0x00,0x96,0x96,0x56,0x56,0xc0,0xc0,0x33,0x33,0xa5,0xa5,0x65,0x65,0xf3,0xf3,
0x00,0x00,0x08,0x08,0xfd,0xfd,0xf5,0xf5,0x77,0x77,0x7f,0x7f,0x8a,0x8a,0x82,0x82,
0x00,0xf5,0xfe,0x0b,0x96,0x63,0x68,0x9d,0x66,0x93,0x98,0x6d,0xf0,0x05,0x0e,0xfb,
0x00,0xea,0x26,0xcc,0x88,0x62,0xae,0x44,0x01,0xeb,0x27,0xcd,0x89,0x63,0xaf,0x45,
0x00,0x0f,0x82,0x8d,0x97,0x98,0x15,0x1a,0xbf,0xb0,0x3d,0x32,0x28,0x27,0xaa,0xa5,
0x00,0x97,0xe3,0x74,0xa0,0x37,0x43,0xd4,0xbe,0x29,0x5d,0xca,0x1e,0x89,0xfd,0x6a,
0x00,0x4a,0x1d,0x57,0x5b,0x11,0x46,0x0c,0x23,0x69,0x3e,0x74,0x78,0x32,0x65,0x2f,
0x00,0xeb,0x93,0x78,0x8c,0x67,0x1f,0xf4,0xb8,0x53,0x2b,0xc0,0x34,0xdf,0xa7,0x4c,
0x00,0x10,0x94,0x84,0xbf,0xaf,0x2b,0x3b,0xe2,0xf2,0x76,0x66,0x5d,0x4d,0xc9,0xd9,
0x00,0xd1,0x9c,0x4d,0x52,0x83,0xce,0x1f,0xbc,0x6d,0x20,0xf1,0xee,0x3f,0x72,0xa3,
0x00,0x41,0xd0,0x91,0xdf,0x9e,0x0f,0x4e,0xa8,0xe9,0x78,0x39,0x77,0x36,0xa7,0xe6,
0x00,0x8e,0x98,0x16,0x6d,0xe3,0xf5,0x7b,0x66,0xe8,0xfe,0x70,0x0b,0x85,0x93,0x1d,
0x00,0x81,0xb6,0x37,0xdd,0x5c,0x6b,0xea,0xf3,0x72,0x45,0xc4,0x2e,0xaf,0x98,0x19,
0x00,0x5d,0x35,0x68,0xab,0xf6,0x9e,0xc3,0x81,0xdc,0xb4,0xe9,0x2a,0x77,0x1f,0x42,
0x00,0x6b,0xc9,0xa2,0x9a,0xf1,0x53,0x38,0x04,0x6f,0xcd,0xa6,0x9e,0xf5,0x57,0x3c,
0x00,0x35,0x94,0xa1,0x3d,0x08,0xa9,0x9c,0x97,0xa2,0x03,0x36,0xaa,0x9f,0x3e,0x0b,
0x00,0xf3,0x7b,0x88,0xb7,0x44,0xcc,0x3f,0x75,0x86,0x0e,0xfd,0xc2,0x31,0xb9,0x4a,
0x00,0x7a,0x1e,0x64,0xbc,0xc6,0xa2,0xd8,0x5f,0x25,0x41,0x3b,0xe3,0x99,0xfd,0x87,
0x00,0xfe,0x96,0x68,0x60,0x9e,0xf6,0x08,0x9c,0x62,0x0a,0xf4,0xfc,0x02,0x6a,0x94,
0x00,0xe3,0xef,0x0c,0x4e,0xad,0xa1,0x42,0x45,0xa6,0xaa,0x49,0x0b,0xe8,0xe4,0x07,
0x00,0xb6,0x93,0x25,0x39,0x8f,0xaa,0x1c,0x28,0x9e,0xbb,0x0d,0x11,0xa7,0x82,0x34,
0x00,0xd9,0x9a,0x43,0x69,0xb0,0xf3,0x2a,0x03,0xda,0x99,0x40,0x6a,0xb3,0xf0,0x29,
0x00,0x74,0xd2,0xa6,0x2a,0x5e,0xf8,0x8c,0x79,0x0d,0xab,0xdf,0x53,0x27,0x81,0xf5,
0x00,0x70,0x64,0x14,0xb9,0xc9,0xdd,0xad,0x04,0x74,0x60,0x10,0xbd,0xcd,0xd9,0xa9,
0x00,0x3f,0xdb,0xe4,0xda,0xe5,0x01,0x3e,0x80,0xbf,0x5b,0x64,0x5a,0x65,0x81,0xbe,
0x00,0xb3,0x2a,0x99,0x59,0xea,0x73,0xc0,0x3c,0x8f,0x16,0xa5,0x65,0xd6,0x4f,0xfc,
0x00,0x05,0x0b,0x0e,0x4d,0x48,0x46,0x43,0xb9,0xbc,0xb2,0xb7,0xf4,0xf1,0xff,0xfa,
0x00,0x3d,0x0c,0x31,0x4b,0x76,0x47,0x7a,0x48,0x75,0x44,0x79,0x03,0x3e,0x0f,0x32,
0x00,0xaa,0x5e,0xf4,0xb5,0x1f,0xeb,0x41,0x8a,0x20,0xd4,0x7e,0x3f,0x95,0x61,0xcb,
0x00,0xcd,0x6d,0xa0,0x6c,0xa1,0x01,0xcc,0x3e,0xf3,0x53,0x9e,0x52,0x9f,0x3f,0xf2,
0x00,0x01,0xd7,0xd6,0x42,0x43,0x95,0x94,0xf9,0xf8,0x2e,0x2f,0xbb,0xba,0x6c,0x6d,
0x00,0x00,0x8f,0x8f,0xdd,0xdd,0x52,0x52,0xa8,0xa8,0x27,0x27,0x75,0x75,0xfa,0xfa,
0x00,0xcb,0xb3,0x78,0xca,0x01,0x79,0xb2,0x42,0x89,0xf1,0x3a,0x88,0x43,0x3b,0xf0,
0x00,0xad,0x9b,0x36,0xab,0x06,0x30,0x9d,0xf0,0x5d,0x6b,0xc6,0x5b,0xf6,0xc0,0x6d,
0x00,0xf6,0xa1,0x57,0x17,0xe1,0xb6,0x40,0x3c,0xca,0x9d,0x6b,0x2b,0xdd,0x8a,0x7c,
0x00,0x86,0x31,0xb7,0x0d,0x8b,0x3c,0xba,0x93,0x15,0xa2,0x24,0x9e,0x18,0xaf,0x29,
0x00,0x69,0x29,0x40,0x5a,0x33,0x73,0x1a,0x6d,0x04,0x44,0x2d,0x37,0x5e,0x1e,0x77,
0x00,0x14,0xe4,0xf0,0x53,0x47,0xb7,0xa3,0x40,0x54,0xa4,0xb0,0x13,0x07,0xf7,0xe3,
0x00,0xab,0x5b,0xf0,0xa1,0x0a,0xfa,0x51,0x9a,0x31,0xc1,0x6a,0x3b,0x90,0x60,0xcb,
0x00,0x5a,0x95,0xcf,0xd4,0x8e,0x41,0x1b,0x4f,0x15,0xda,0x80,0x9b,0xc1,0x0e,0x54,
0x00,0xb2,0x50,0xe2,0xf1,0x43,0xa1,0x13,0x56,0xe4,0x06,0xb4,0xa7,0x15,0xf7,0x45,
0x00,0x47,0xe7,0xa0,0x9b,0xdc,0x7c,0x3b,0x91,0xd6,0x76,0x31,0x0a,0x4d,0xed,0xaa,
0x00,0xaa,0x5d,0xf7,0x96,0x3c,0xcb,0x61,0x7b,0xd1,0x26,0x8c,0xed,0x47,0xb0,0x1a,
0x00,0x3a,0xe4,0xde,0x53,0x69,0xb7,0x8d,0xaa,0x90,0x4e,0x74,0xf9,0xc3,0x1d,0x27,
0x00,0xd3,0x84,0x57,0x16,0xc5,0x92,0x41,0x5b,0x88,0xdf,0x0c,0x4d,0x9e,0xc9,0x1a,
0x00,0x42,0x9e,0xdc,0xe3,0xa1,0x7d,0x3f,0x0f,0x4d,0x91,0xd3,0xec,0xae,0x72,0x30,
0x00,0x4b,0xc5,0x8e,0x60,0x2b,0xa5,0xee,0xa0,0xeb,0x65,0x2e,0xc0,0x8b,0x05,0x4e,
0x00,0x25,0x32,0x17,0x86,0xa3,0xb4,0x91,0x11,0x34,0x23,0x06,0x97,0xb2,0xa5,0x80,
0x00,0xd3,0x1f,0xcc,0x01,0xd2,0x1e,0xcd,0xf9,0x2a,0xe6,0x35,0xf8,0x2b,0xe7,0x34,
0x00,0xde,0xbd,0x63,0xfa,0x24,0x47,0x99,0xa2,0x7c,0x1f,0xc1,0x58,0x86,0xe5,0x3b,
0x00,0xcd,0x6b,0xa6,0xc6,0x0b,0xad,0x60,0xf4,0x39,0x9f,0x52,0x32,0xff,0x59,0x94,
0x00,0x6a,0x8d,0xe7,0x9a,0xf0,0x17,0x7d,0x7c,0x16,0xf1,0x9b,0xe6,0x8c,0x6b,0x01,
0x00,0x92,0xec,0x7e,0x2e,0xbc,0xc2,0x50,0x50,0xc2,0xbc,0x2e,0x7e,0xec,0x92,0x00,
0x00,0xc8,0xa8,0x60,0x4b,0x83,0xe3,0x2b,0x0e,0xc6,0xa6,0x6e,0x45,0x8d,0xed,0x25,
0x00,0x5b,0x06,0x5d,0xf5,0xae,0xf3,0xa8,0x21,0x7a,0x27,0x7c,0xd4,0x8f,0xd2,0x89,
0x00,0x69,0x06,0x6f,0xfe,0x97,0xf8,0x91,0x7d,0x14,0x7b,0x12,0x83,0xea,0x85,0xec,
0x00,0x58,0xc1,0x99,0x6c,0x34,0xad,0xf5,0x00,0x58,0xc1,0x99,0x6c,0x34,0xad,0xf5,
0x00,0xe9,0xc8,0x21,0x8e,0x67,0x46,0xaf,0x88,0x61,0x40,0xa9,0x06,0xef,0xce,0x27,
0x00,0x97,0xbe,0x29,0x49,0xde,0xf7,0x60,0x46,0xd1,0xf8,0x6f,0x0f,0x98,0xb1,0x26,
0x00,0x83,0x18,0x9b,0xa4,0x27,0xbc,0x3f,0x8d,0x0e,0x95,0x16,0x29,0xaa,0x31,0xb2,
0x00,0xa4,0xf2,0x56,0x81,0x25,0x73,0xd7,0x72,0xd6,0x80,0x24,0xf3,0x57,0x01,0xa5,
0x00,0x96,0x99,0x0f,0x34,0xa2,0xad,0x3b,0x57,0xc1,0xce,0x58,0x63,0xf5,0xfa,0x6c,
0x00,0x97,0x6a,0xfd,0xd8,0x4f,0xb2,0x25,0x49,0xde,0x23,0xb4,0x91,0x06,0xfb,0x6c,
0x00,0x9b,0x63,0xf8,0x21,0xba,0x42,0xd9,0x58,0xc3,0x3b,0xa0,0x79,0xe2,0x1a,0x81,
0x00,0xa8,0xdd,0x75,0xb8,0x10,0x65,0xcd,0x39,0x91,0xe4,0x4c,0x81,0x29,0x5c,0xf4,
0x00,0x17,0x1b,0x0c,0x7c,0x6b,0x67,0x70,0xff,0xe8,0xe4,0xf3,0x83,0x94,0x98,0x8f,
0x00,0xe7,0xdb,0x3c,0x38,0xdf,0xe3,0x04,0xa5,0x42,0x7e,0x99,0x9d,0x7a,0x46,0xa1,
0x00,0x2d,0x6e,0x43,0x58,0x75,0x36,0x1b,0x84,0xa9,0xea,0xc7,0xdc,0xf1,0xb2,0x9f,
0x00,0xd0,0x67,0xb7,0x71,0xa1,0x16,0xc6,0x89,0x59,0xee,0x3e,0xf8,0x28,0x9f,0x4f,
0x00,0xb4,0x6b,0xdf,0x8b,0x3f,0xe0,0x54,0xc0,0x74,0xab,0x1f,0x4b,0xff,0x20,0x94,
0x00,0xf7,0x94,0x63,0x72,0x85,0xe6,0x11,0x82,0x75,0x16,0xe1,0xf0,0x07,0x64,0x93,
0x00,0xc0,0xca,0x0a,0xa7,0x67,0x6d,0xad,0x04,0xc4,0xce,0x0e,0xa3,0x63,0x69,0xa9,
0x00,0xd7,0x9f,0x48,0x35,0xe2,0xaa,0x7d,0xef,0x38,0x70,0xa7,0xda,0x0d,0x45,0x92,
0x00,0x9c,0x77,0xeb,0x52,0xce,0x25,0xb9,0x70,0xec,0x07,0x9b,0x22,0xbe,0x55,0xc9,
0x00,0x0b,0xe1,0xea,0xe8,0xe3,0x09,0x02,0xb6,0xbd,0x57,0x5c,0x5e,0x55,0xbf,0xb4,
0x00,0x15,0xe0,0xf5,0x3f,0x2a,0xdf,0xca,0xce,0xdb,0x2e,0x3b,0xf1,0xe4,0x11,0x04,
0x00,0xb5,0x1f,0xaa,0xc9,0x7c,0xd6,0x63,0xc8,0x7d,0xd7,0x62,0x01,0xb4,0x1e,0xab,
0x00,0xe2,0x05,0xe7,0x5c,0xbe,0x59,0xbb,0xb1,0x53,0xb4,0x56,0xed,0x0f,0xe8,0x0a,
0x00,0xbf,0x56,0xe9,0xd7,0x68,0x81,0x3e,0x1f,0xa0,0x49,0xf6,0xc8,0x77,0x9e,0x21,
0x00,0x3d,0x73,0x4e,0x98,0xa5,0xeb,0xd6,0xee,0xd3,0x9d,0xa0,0x76,0x4b,0x05,0x38,
0x00,0x1c,0x76,0x6a,0x5c,0x40,0x2a,0x36,0x25,0x39,0x53,0x4f,0x79,0x65,0x0f,0x13,
0x00,0xbb,0xee,0x55,0x1b,0xa0,0xf5,0x4e,0x4c,0xf7,0xa2,0x19,0x57,0xec,0xb9,0x02,
0x00,0x51,0x94,0xc5,0xc3,0x92,0x57,0x06,0x83,0xd2,0x17,0x46,0x40,0x11,0xd4,0x85,
0x00,0x9e,0x25,0xbb,0x31,0xaf,0x14,0x8a,0x46,0xd8,0x63,0xfd,0x77,0xe9,0x52,0xcc,
0x00,0x7f,0x9e,0xe1,0x66,0x19,0xf8,0x87,0x43,0x3c,0xdd,0xa2,0x25,0x5a,0xbb,0xc4,
0x00,0x06,0x35,0x33,0x02,0x04,0x37,0x31,0x58,0x5e,0x6d,0x6b,0x5a,0x5c,0x6f,0x69,
0x00,0x79,0xae,0xd7,0x64,0x1d,0xca,0xb3,0x99,0xe0,0x37,0x4e,0xfd,0x84,0x53,0x2a,
0x00,0x00,0xa5,0xa5,0x68,0x68,0xcd,0xcd,0xd5,0xd5,0x70,0x70,0xbd,0xbd,0x18,0x18,
0x00,0x31,0xe2,0xd3,0xe8,0xd9,0x0a,0x3b,0x4c,0x7d,0xae,0x9f,0xa4,0x95,0x46,0x77,
0x00,0xb7,0x4b,0xfc,0x0f,0xb8,0x44,0xf3,0x5b,0xec,0x10,0xa7,0x54,0xe3,0x1f,0xa8,
0x00,0x07,0x97,0x90,0x9a,0x9d,0x0d,0x0a,0x6b,0x6c,0xfc,0xfb,0xf1,0xf6,0x66,0x61,
0x00,0x87,0x59,0xde,0xf9,0x7e,0xa0,0x27,0x80,0x07,0xd9,0x5e,0x79,0xfe,0x20,0xa7,
0x00,0x84,0x7a,0xfe,0xae,0x2a,0xd4,0x50,0x6e,0xea,0x14,0x90,0xc0,0x44,0xba,0x3e,
0x00,0x9a,0x6a,0xf0,0xa5,0x3f,0xcf,0x55,0xc1,0x5b,0xab,0x31,0x64,0xfe,0x0e,0x94,
0x00,0x53,0x3b,0x68,0x00,0x53,0x3b,0x68,0x2a,0x79,0x11,0x42,0x2a,0x79,0x11,0x42,
0x00,0xb2,0xd1,0x63,0xfb,0x49,0x2a,0x98,0x26,0x94,0xf7,0x45,0xdd,0x6f,0x0c,0xbe,
0x00,0x50,0x50,0x00,0xd8,0x88,0x88,0xd8,0xbd,0xed,0xed,0xbd,0x65,0x35,0x35,0x65,
0x00,0xe5,0x6e,0x8b,0x7f,0x9a,0x11,0xf4,0xcd,0x28,0xa3,0x46,0xb2,0x57,0xdc,0x39,
0x00,0xe0,0xd0,0x30,0xdc,0x3c,0x0c,0xec,0xab,0x4b,0x7b,0x9b,0x77,0x97,0xa7,0x47,
0x00,0x46,0xcd,0x8b,0xed,0xab,0x20,0x66,0x4f,0x09,0x82,0xc4,0xa2,0xe4,0x6f,0x29,
0x00,0x76,0x9b,0xed,0xef,0x99,0x74,0x02,0xe2,0x94,0x79,0x0f,0x0d,0x7b,0x96,0xe0,
0x00,0x4f,0x5a,0x15,0xa7,0xe8,0xfd,0xb2,0xcb,0x84,0x91,0xde,0x6c,0x23,0x36,0x79,
0x00,0x7b,0x60,0x1b,0x30,0x4b,0x50,0x2b,0xdf,0xa4,0xbf,0xc4,0xef,0x94,0x8f,0xf4,
0x00,0x42,0x62,0x20,0x66,0x24,0x04,0x46,0xf3,0xb1,0x91,0xd3,0x95,0xd7,0xf7,0xb5,
0x00,0x21,0xf2,0xd3,0x04,0x25,0xf6,0xd7,0x69,0x48,0x9b,0xba,0x6d,0x4c,0x9f,0xbe,
0x00,0xc0,0xd7,0x17,0x0b,0xcb,0xdc,0x1c,0xa8,0x68,0x7f,0xbf,0xa3,0x63,0x74,0xb4,
0x00,0xde,0x0d,0xd3,0xb1,0x6f,0xbc,0x62,0x66,0xb8,0x6b,0xb5,0xd7,0x09,0xda,0x04,
0x00,0xb1,0x3d,0x8c,0xa4,0x15,0x99,0x28,0xe2,0x53,0xdf,0x6e,0x46,0xf7,0x7b,0xca,
0x00,0xe9,0xae,0x47,0x03,0xea,0xad,0x44,0x70,0x99,0xde,0x37,0x73,0x9a,0xdd,0x34,
0x00,0x3b,0x24,0x1f,0x7f,0x44,0x5b,0x60,0xb3,0x88,0x97,0xac,0xcc,0xf7,0xe8,0xd3,
0x00,0x32,0x87,0xb5,0xf0,0xc2,0x77,0x45,0x3b,0x09,0xbc,0x8e,0xcb,0xf9,0x4c,0x7e,
0x00,0x50,0x4d,0x1d,0x17,0x47,0x5a,0x0a,0xaf,0xff,0xe2,0xb2,0xb8,0xe8,0xf5,0xa5,
0x00,0x01,0x29,0x28,0xb5,0xb4,0x9c,0x9d,0xc5,0xc4,0xec,0xed,0x70,0x71,0x59,0x58,
0x00,0xa0,0x01,0xa1,0x4a,0xea,0x4b,0xeb,0xfd,0x5d,0xfc,0x5c,0xb7,0x17,0xb6,0x16,
0x00,0x38,0x7a,0x42,0xc1,0xf9,0xbb,0x83,0x74,0x4c,0x0e,0x36,0xb5,0x8d,0xcf,0xf7,
0x00,0x52,0x35,0x67,0x5a,0x08,0x6f,0x3d,0xdc,0x8e,0xe9,0xbb,0x86,0xd4,0xb3,0xe1,
0x00,0x28,0x27,0x0f,0x3e,0x16,0x19,0x31,0xd0,0xf8,0xf7,0xdf,0xee,0xc6,0xc9,0xe1,
0x00,0xa8,0x22,0x8a,0x1c,0xb4,0x3e,0x96,0x68,0xc0,0x4a,0xe2,0x74,0xdc,0x56,0xfe,
0x00,0xef,0xe4,0x0b,0x7e,0x91,0x9a,0x75,0xcb,0x24,0x2f,0xc0,0xb5,0x5a,0x51,0xbe,
0x00,0x63,0xac,0xcf,0x10,0x73,0xbc,0xdf,0xb5,0xd6,0x19,0x7a,0xa5,0xc6,0x09,0x6a,
0x00,0xff,0xa0,0x5f,0x08,0xf7,0xa8,0x57,0xdb,0x24,0x7b,0x84,0xd3,0x2c,0x73,0x8c,
0x00,0x2a,0xfb,0xd1,0xc2,0xe8,0x39,0x13,0xb2,0x98,0x49,0x63,0x70,0x5a,0x8b,0xa1,
};