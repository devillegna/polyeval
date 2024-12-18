
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
void encode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
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
void decode_64( uint64_t * rfx , const uint64_t * fx , unsigned n_u64 )
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


