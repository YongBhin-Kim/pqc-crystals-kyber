#include <stdint.h>
#include "params.h"
#include "poly.h"
// #include "ntt.h"     // mul
#include "reduce.h"     // reduce
#include "cbd.h"        // cbd (used in sample vector)
#include "symmetric.h"  // shake128_xof and shake256_prf
#include "stdint.h"

#include "print.h"

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials
*
* Arguments: - poly *res:      pointer to output polynomial
*            - const poly *p1: pointer to first input polynomial
*            - const poly *p2: pointer to second input polynomial
**************************************************/
void poly_add(poly *res, poly *p1, poly *p2) {

    for (unsigned int i=0; i<N; i++) {
        res->coeffs[i] = p1->coeffs[i] + p2->coeffs[i];
    }
    poly_reduce(res); 
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials
*
* Arguments: - poly *res:      pointer to output polynomial
*            - const poly *p1: pointer to first input polynomial
*            - const poly *p2: pointer to second input polynomial
**************************************************/
void poly_sub(poly *res, poly *p1, poly *p2) {

    for (unsigned int i=0; i<N; i++)
        res->coeffs[i] = p1->coeffs[i] - p2->coeffs[i];
    poly_reduce(res); // 빼도 될지
}

/*************************************************
* Name:        poly_csubq
*
* Description: Applies conditional subtraction of q to each coefficient
*              of a polynomial. For details of conditional subtraction
*              of q see comments in reduce.c
*
* Arguments:   - poly *res: pointer to input/output polynomial
    
    return value is 'coeffs - Q' if coeffs >= Q, 
                        else 'coeffs'
**************************************************/
void poly_csubq(poly *res) {

    for (unsigned int i=0; i<N; i++)
        res->coeffs[i] = csubq(res->coeffs[i]);
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies reduction to all coefficients of a polynomial
*              for details of the reduction see comments in reduce.c
*
* Arguments:   - poly *res: pointer to input/output polynomial

Kyber uses barrett reduce
**************************************************/
void poly_reduce(poly *res) {

    for (unsigned int i=0; i<N; i++)
        res->coeffs[i] = simple_reduce(res->coeffs[i]);
}

/*************************************************
* Name:        poly_mul
*
* Description: multiplation two polynomials
*
* Arguments:   - poly *res: pointer to output polynomial
               - poly *p1:  pointer to first  input polynomial
               - poly *p1:  pointer to second input polynomial
Kyber uses ntt
**************************************************/
void poly_mul(poly *res, poly *p1, poly *p2) {

    int32_t tempcoef = 0;
    int16_t temppoly[N * 2] = { 0, };

    poly_reduce(p1);
    poly_reduce(p2);

    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {

            tempcoef = (int32_t) p1->coeffs[i] * (int32_t) p2->coeffs[j];
            tempcoef = tempcoef + (int32_t) temppoly[i+j]; // < 32-bits
            temppoly[i+j] =  (int16_t) (simple_reduce32(tempcoef));
        }
    }

    for (int i=N; i<N*2; i++)
        res->coeffs[i-N] = simple_reduce(temppoly[i-N] - temppoly[i]);
}

/*************************************************
* Name:        poly_tobytes (Encode_12)
*
* Description: Serialization of a polynomial
*
* Arguments:   - uint8_t *res: pointer to output byte array
*                            (needs space for POLYBYTES bytes)
*              - poly *p:    pointer to input polynomial

 pk <- Encode_12 ( b[0] )
 pk is 12*K*N/8-bytes and concate publicseed

 sk <- Encode_12 ( s[0] )
 sk is 12*K*N/8-bytes
**************************************************/
void poly_tobytes(uint8_t res[POLYBYTES], poly *p) {

    uint16_t t0, t1;

    poly_csubq(p);
    for (uint8_t i=0; i<N/2; i++) {
        t0 = p->coeffs[2*i+0];
        t1 = p->coeffs[2*i+1];

        res[3*i + 0] = (t0 >> 0);
        res[3*i + 1] = (t0 >> 8) | (t1 << 4);
        res[3*i + 2] = (t1 >> 4);
    }
}

/*************************************************
* Name:        poly_frombytes (Decode_12)
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *ers:          pointer to output polynomial
*              - const uint8_t *a:   pointer to input byte array
*                                   (of POLYBYTES bytes)
    b <- Decode_12 ( pk )
    s <- Decode_12 ( sk )
**************************************************/
void poly_frombytes(poly *res, const uint8_t a[POLYBYTES]) {
    
    for (uint8_t i=0; i<N/2; i++) {
        res->coeffs[2*i+0] = ((a[3*i+0]>>0) | ((uint16_t)a[3*i+1]<<8)) & 0b111111111111;
        res->coeffs[2*i+1] = ((a[3*i+1]>>4) | ((uint16_t)a[3*i+2]<<4)) & 0b111111111111;
    }
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *res:            pointer to output polynomial
*              - const uint8_t *m: pointer to input message

Decompress_Q( Decode_1(m, 1) )

    256-element coeffs <- 32-byte message 
    if m_indexbit is 0
        then coeffs_index is 0
    else if m_indexbit is 1
        then coeffs_index is Q/2 == 1665
**************************************************/
void poly_frommsg(poly *res, const uint8_t m[SYMBYTES]) {

    for (unsigned int i=0; i<N/8; i++) {
        for (unsigned int j=0; j<8; j++) { 
            res->coeffs[8*i + j] = 0;
            if ( (m[i] >> j) & 0x01 )
                res->coeffs[8*i+j] = (Q+1)/2;
        }
    }
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *res: pointer to output message
*              - poly *p:      pointer to input polynomial

    32-byte message <- 256-element coeffs
    If coeffs_index is near by 0
        then m_indexbit is 0
    else if coeffs_index is near by Q/2
        then m_indexbit is 1
**************************************************/
void poly_tomsg(uint8_t res[SYMBYTES], poly *p) {
    
    poly_csubq(p);
    uint16_t t;

    for (unsigned int i=0; i<N/8; i++) {
        res[i] = 0;
        for (unsigned int j=0; j<8; j++) {
            t = ((((uint16_t)p->coeffs[8*i+j] << 1) + Q/2)/Q) & 1;
            res[i] |= t << j;
        }
    }
}

/*************************************************
* Name:        poly_compress == Encode_4(Compress_Q(ct2, 4))
*
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - uint8_t *res: pointer to output byte array
*                            (of length POLYCOMPRESSEDBYTES)
*              - poly *p:    pointer to input polynomial

    
    Used in compress 'ct_v'
    byte array 'res' (128-bytes) <- poly 'ct_v' (512-bytes) (in Kyber512)
    
    c2 <- Encode_4( Compress_Q(ct_v, 4) )                   (for Kyber512, 'd_v' is 4)
        (4-bits entropy <- 12-bits entropy)
**************************************************/
void poly_compress(uint8_t res[POLYCOMPRESSEDBYTES], poly *p) { 

    uint8_t t[8]; // depended on (LCM(d_v, 8) / d_v) (parallel 8-times)

    poly_csubq(p);

#if (POLYCOMPRESSEDBYTES == 128)    // Kyber512/768
    for ( unsigned int i=0; i<N/8; i++ ) {
        for ( unsigned int j=0; j<8; j++ ) {
            /* t <- Compress_Q( ct_v, d_v ) */
            t[j] = ((( (uint16_t)p->coeffs[8*i+j] << 4 ) + Q/2)/Q) & 15;
        }
        
        /* ct <- Encode_d_v( t ) */
        res[0] = t[0] | (t[1] << 4);
        res[1] = t[2] | (t[3] << 4);
        res[2] = t[4] | (t[5] << 4);
        res[3] = t[6] | (t[7] << 4); 
        res += 4;
    }
#elif (POLYCOMPRESSEDBYTES == 160)  // Kyber1024
    for (uint8_t i=0; i<N/8; i++) {
        for (uint8_t j=0; j<8; j++)
            t[j] = ((((uint32_t)p->coeffs[8*i+j] << 5) + Q/2)/Q) & 0b11111;
        res[0] = (t[0] >> 0) | (t[1] << 5);
        res[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
        res[2] = (t[3] >> 1) | (t[4] << 4);
        res[3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
        res[4] = (t[6] >> 2) | (t[7] << 3);

        rest += 5;
    }
#else
#error "POLYCOMPRESSEDBYTES needs to be in Kyber512/768 = 128, Kyber1024 = 160 becase dv = 4 or 5"
#endif
}

/*************************************************
* Name:        poly_decompress == Encode_10(Compress_Q(ct2, 10))
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
*
* Arguments:   - poly *res:          pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of length POLYCOMPRESSEDBYTES bytes)

    Function poly.decompress
    
    Used in decompress ct_v
    poly ct_v (512-bytes) <- byte array (128-bytes) (in Kyber512)
    
    ct_v <- Decompress_Q( Decode_4(c2, 4) )         (for Kyber512, 'd_v' is 4)
**************************************************/
void poly_decompress(poly *res, const uint8_t a[POLYCOMPRESSEDBYTES]) {

    // depended on (LCM(d_u, 8) / d_u)

#if (POLYCOMPRESSEDBYTES == 128)
    for ( unsigned int i=0; i<N/2; i++ ) {
        /*  ct_v <-  Decompress_Q( Decode_d_u(ct), d_v ) */
        res->coeffs[2*i+0] = (((uint16_t)(a[0] & 15)*Q) + 8) >> 4; //    [ a7 a6 a5 a4 || a3 a2 a1 a0  ]
        res->coeffs[2*i+1] = (((uint16_t)(a[0] >> 4)*Q) + 8) >> 4; // -> [   coeffs1   ,    coeffs0    ]
        a += 1;
    }
#else
#error "POLYCOMPRESSEDBYTES needs to be in Kyber512/768 = 128, Kyber1024 = 160 becase dv = 4 or 5"
#endif

}

/*************************************************
* Name:        poly_getnoise_eta1 (PRF, CBD)
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA1
*
* Arguments:   - poly *res:            pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce

Case : ETA is 3
**************************************************/
void poly_getnoise_eta1(poly *res, const uint8_t seed[SYMBYTES], uint8_t nonce) {
    
    uint8_t buf[ETA1 * N/4];
    kyber_shake256_prf(buf, sizeof(buf), seed, nonce);
    cbd_eta1(res, buf);
}

/*************************************************
* Name:        poly_getnoise_eta2
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
*
* Arguments:   - poly *res:           pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce

Case : ETA is 2
**************************************************/
void poly_getnoise_eta2(poly *res, const uint8_t seed[SYMBYTES], uint8_t nonce) {
    
    uint8_t buf[ETA2 * N/4];
    kyber_shake256_prf(buf, sizeof(buf), seed, nonce);
    cbd_eta2(res, buf);
}
