#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

#include "print.h"
#include <string.h>

/*************************************************
* Name:        polyvec_add
*
* Description: Add vectors of polynomials
*
* Arguments: - polyvec *res:       pointer to output vector of polynomials
*            - const polyvec *pv1: pointer to first input vector of polynomials
*            - const polyvec *pv2: pointer to second input vector of polynomials
**************************************************/
void polyvec_add(polyvec *res, polyvec *pv1, polyvec *pv2) {

    for ( unsigned int i=0; i<K; i++ )
        poly_add(&res->vec[i], &pv1->vec[i], &pv2->vec[i]);
}

/*************************************************
* Name:        polyvec_csubq
*
* Description: Applies conditional subtraction of q to each coefficient
*              of each element of a vector of polynomials
*              for details of conditional subtraction of q see comments in
*              reduce.c
*
* Arguments:   - poly *res: pointer to input/output polynomial
**************************************************/
void polyvec_csubq(polyvec *res) {

    for ( unsigned int i=0; i<K; i++ )
        poly_csubq(&res->vec[i]);
}

/*************************************************
* Name:        polyvec_reduce
*
* Description: Applies modulo reduction of Q to each coefficient
*              of each element of a vector of polynomials
*              for details of reduce of Q see comments in
*              reduce.c
*
* Arguments:   - poly *res: pointer to input/output polynomial

Kyber uses barrett reduce
**************************************************/
void polyvec_reduce(polyvec *res) {

    for ( unsigned int i=0; i<K; i++ )
        poly_reduce(&res->vec[i]);
}

/*************************************************
* Name:        polyvec_mul
*
* Description: multiplation two polynomial vectors
*
* Arguments:   - poly *res:     pointer to output polynomial
               - polyvec *pv1:  pointer to first  input polynomial vector
               - polyvec *pv2:  pointer to second input polynomial vector
Kyber uses ntt
**************************************************/
void polyvec_mul(poly *res, polyvec *pv1, polyvec *pv2) {

    poly temppoly;

    for ( unsigned int i=0; i<K; i++ ) {
        poly_mul(&temppoly, &pv1->vec[i], &pv2->vec[i]);
        poly_add(res, res, &temppoly);
    }
}

/*************************************************
* Name:        polyvec_tobytes (Encode_12)
*
* Description: Serialize vector of polynomials
*
* Arguments:   - uint8_t *res: pointer to output byte array
*                            (needs space for POLYVECBYTES)
*              - polyvec *pv: pointer to input vector of polynomials

used to pack public key 'pk' and secret key 'sk' forms of entry of polynomial vector 'b', 's'

pk <- Encode_12 ( b[0] )
pk is 12*K*N/8-bytes and concate publicseed

sk <- Encode_12 ( s[0] )
sk is 12*K*N/8-bytes
**************************************************/
void polyvec_tobytes(uint8_t res[POLYVECBYTES], polyvec *pv) {

    for ( unsigned int i=0; i<K; i++ )
        poly_tobytes(res + (POLYBYTES*i), &pv->vec[i]);

}

/*************************************************
* Name:        polyvec_frombytes (Decode_12)
*
* Description: De-serialize vector of polynomials;
*              inverse of polyvec_tobytes
*
* Arguments:   - uint8_t *res:      pointer to output byte array
*              - const polyvec *a:  pointer to input vector of polynomials
*                                   (of length POLYVECBYTES)

used to unpack public key 'b' and secret key 's' forms of entry of byte array 'pk', 'sk'

b <- Decode_12 ( pk )
s <- Decode_12 ( sk )
**************************************************/
void polyvec_frombytes(polyvec *res, const uint8_t a[POLYVECBYTES]) {

    for ( unsigned int i=0; i<K; i++ )
        poly_frombytes(&res->vec[i], a + (POLYBYTES*i));

}

/*************************************************
* Name:        polyvec_compress (Encode and Compress)
*
* Description: Compress and serialize vector of polynomials
*
* Arguments:   - uint8_t *res: pointer to output byte array
*                            (needs space for POLYVECCOMPRESSEDBYTES)
*              - polyvec *pv: pointer to input vector of polynomials

    Used in compress ct_u
    byte array C1 (640-bytes) <- polyvec ct_u (1024-bytes)
    
    C1 <- Encode_d_u( Compress_Q( ct_u, d_u ) )
        (d_u-bits entropy <- 12-bits entropy)
**************************************************/
void polyvec_compress(uint8_t res[POLYVECCOMPRESSEDBYTES], polyvec *pv) {

    uint16_t t[4]; // depended on (LCM(d_u, 8) / d_u)

    polyvec_csubq(pv);

    for ( unsigned int i=0; i<K; i++ ) {
        for ( unsigned int j=0; j<N/4; j++ ) {
            for ( unsigned int k=0; k<4; k++ ) {
                /* t <- Compress_Q( ct_u, d_u ) */
                t[k] = ((( (uint32_t)(pv->vec[i].coeffs[4*j+k]) << 10 ) + Q/2)/Q) & 0b1111111111;
            }
            
            /* ct <- Encode_d_u( t ) */
            res[0] = (t[0] >> 0);
            res[1] = (t[0] >> 8) | (t[1] << 2);
            res[2] = (t[1] >> 6) | (t[2] << 4);
            res[3] = (t[2] >> 4) | (t[3] << 6);
            res[4] = (t[3] >> 2);
            res += 5;
        }
    }

}

/*************************************************
* Name:        polyvec_decompress (Decode and Decompress)
*
* Description: De-serialize and decompress vector of polynomials;
*              approximate inverse of polyvec_compress
*
* Arguments:   - polyvec *res:       pointer to output vector of polynomials
*              - const uint8_t *a:   pointer to input byte array
*                                    (of length POLYVECCOMPRESSEDBYTES)

    Used in decompress C1
    polyvec ct_u (1024-bytes) <- byte array C1 (640-bytes in kyber512) 
    
    ct_u <- Decompress_Q( Decode_10(C1), 10 ) (kyber512 : d_u = 10)
**************************************************/
void polyvec_decompress(polyvec *res, const uint8_t a[POLYVECCOMPRESSEDBYTES]) { 
        
    uint16_t t[4]; // depended on (LCM(d_u, 8) / d_u)

    for ( unsigned int i=0; i<K; i++ ) {
        for ( unsigned int j=0; j<N/4; j++ ) {
            /* t <- Decode_d_u( ct ) */
            t[0] = (a[0] >> 0) | ((uint16_t)a[1] << 8);
            t[1] = (a[1] >> 2) | ((uint16_t)a[2] << 6);
            t[2] = (a[2] >> 4) | ((uint16_t)a[3] << 4);
            t[3] = (a[3] >> 6) | ((uint16_t)a[4] << 2);
            a += 5;

            for ( unsigned int k=0; k<4; k++ )
                /* ct_u <- Decompress_Q( t, d_u ) */
                res->vec[i].coeffs[4*j + k] = ((uint32_t)(t[k] & 0x3ff)*Q + 512) >> 10;
        }
    }
}

/*************************************************
* Name:        samplevector_eta1
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter eta=3
*
* Arguments:   - polyvec *res:            pointer to output polynomial vector
*              - const uint8_t *seed:     pointer to input byte array 
*              - uint8_t *nonce:          pointer to input unsigned 8-bits integer array

Sample vector 's, e, r'
**************************************************/
void samplevector_eta1(polyvec *res, const uint8_t seed[SYMBYTES], uint8_t *nonce) {

    for ( unsigned int i=0; i<K; i++ )
        poly_getnoise_eta1(&res->vec[i], seed, (*nonce)++);
}

/*************************************************
* Name:        samplevector_eta2
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter eta=2
*
* Arguments:   - polyvec *res:            pointer to output polynomial vector
*              - const uint8_t *seed:     pointer to input byte array 
*              - uint8_t *nonce:          pointer to input unsigned 8-bits integer array

Sample vector 'e1, e2'
**************************************************/
void samplevector_eta2(polyvec *res, const uint8_t seed[SYMBYTES], uint8_t *nonce) {

    for ( unsigned int i=0; i<K; i++ )
        poly_getnoise_eta2(&res->vec[i], seed, (*nonce)++);
}
