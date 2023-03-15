/// @file utils_hash.h
/// @brief the interface for adapting hash functions.
///
///
#ifndef _UTILS_HASH2_H_
#define _UTILS_HASH2_H_


#include <openssl/evp.h>


#ifdef  __cplusplus
extern  "C" {
#endif


typedef struct hash_ctx {
    EVP_MD_CTX *x;
} hash_ctx;


#ifdef  __cplusplus
}
#endif


#define _HASH_SHAKE128_

static inline
int hash_init( hash_ctx *ctx )
{
    ctx->x = EVP_MD_CTX_create();
    if (!ctx->x) {
        return -1;
    }

    #if defined(_HASH_SHAKE128_)
    int ok = EVP_DigestInit_ex(ctx->x, EVP_shake128(), NULL);
    #else
    int ok = EVP_DigestInit_ex(ctx->x, EVP_shake256(), NULL);
    #endif
    return (ok) ? 0 : -1;
}

static inline
int hash_update( hash_ctx *ctx, const unsigned char *mesg, size_t mlen )
{
    int ok = EVP_DigestUpdate(ctx->x, mesg, mlen);
    return (ok) ? 0 : -1;
}


static inline
int hash_final_digest( unsigned char *digest, size_t dlen, hash_ctx *ctx )     // free ctx
{
    int ok = EVP_DigestFinalXOF(ctx->x, digest, dlen);
    EVP_MD_CTX_destroy(ctx->x);
    return (ok) ? 0 : -1;
}




#endif // _UTILS_HASH_H_

