#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "params.h"     // Define parameters
#include "pke.h"
#include "poly.h"
#include "polyvec.h"
#include "rng.h"        // AES-256 CTR
#include "fips202.h"    // SHA3-512
#include "symmetric.h"  // XOF : SHAKE-128,  PRF : SHAKE-256  (Used in sample matrix and sample vector)
#include "print.h"      // Print

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *res:          pointer to the output serialized public key
*              polyvec *pv:         pointer to the input public-key polyvec ('b')
*              const uint8_t *seed: pointer to the input public seed ('publicseed')

b : polynomial (1024-bytes)

pk <- ( Encode_12(b) || publicseed )
      ( 768-bytes    || 32-bytes )

**************************************************/
static void pack_pk(uint8_t res[PUBLICKEYBYTES], polyvec *pv, const uint8_t seed[SYMBYTES]) {
    polyvec_tobytes(res, pv);
    memcpy(res+POLYVECBYTES, seed, SYMBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *res:  pointer to output serialized secret key
*              - polyvec *pv: pointer to input vector of polynomials (secret key) ('s')

s : polynomial (1024-bytes)

sk <- ( Encode_12(s) )
      ( 768-bytes    )
**************************************************/
static void pack_sk(uint8_t res[SECRETKEYBYTES], polyvec *pv) {
    polyvec_tobytes(res, pv);
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk:             pointer to output public-key
*                                         polynomial vector
*              - uint8_t *seed:           pointer to output seed to generate
*                                         matrix A ('publicseed')
*              - const uint8_t *a: pointer to input serialized public key

pk <- ( Encode_12(b) || publicseed )
      ( 768-bytes    || 32-bytes )

publicseed <- pk + 768 (call by address)

b <- Decode_12( pk )
b : polynomial (1024-bytes)
**************************************************/
static void unpack_pk(polyvec *res1, uint8_t res2[SYMBYTES], const uint8_t a[PUBLICKEYBYTES]) {
    polyvec_frombytes(res1, a);
    memcpy(res2, a+POLYVECBYTES, SYMBYTES);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key;
*              inverse of pack_sk
*
* Arguments:   - polyvec *res:             pointer to output vector of
*                                         polynomials (secret key) ('s')
*              - const uint8_t *a: pointer to input serialized secret key ('sk')

sk <- ( Encode_12(s) )
      ( 768-bytes    )

s <- Decode_12( sk )
s : polynomial (1024-bytes)
**************************************************/
static void unpack_sk(polyvec *res, const uint8_t a[SECRETKEYBYTES]) {
    polyvec_frombytes(res, a);
}

/*************************************************
* Name:        pack_ct
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *res: pointer to the output serialized ciphertext
*              polyvec *pv:     pointer to the input vector of polynomials ('ct_u')
*              poly *p:      pointer to the input polynomial ('ct_v')

C1 <- Encode_d_u( Compress_Q( ct_u, d_u ) ) 
C2 <- Encode_d_v( Compress_Q( ct_v, d_v ) )

ct <- ( C1 || C2 ) // parameter 'r'

C1 : (K * (512 * 10/16) -bytes)   (Depended on d_u) (d_u == 10)
C2 : (512 * 4/16 -bytes)          (Depended on d_v) (d_v == 4)
ct : 512 * 4/16 + K * (512 * 10/16) - bytes    // param 'CIPHERTEXTBYTES'
 r <- ( C1 || C2 )
**************************************************/
static void pack_ct(uint8_t res[CIPHERTEXTBYTES], polyvec *pv, poly *p) {
    polyvec_compress(res, pv);
    poly_compress(res + POLYVECCOMPRESSEDBYTES, p);
}

/*************************************************
* Name:        unpack_ct
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *res1:       pointer to the output vector of polynomials ('ct_u')
*              - poly *res2:          pointer to the output polynomial  ('ct_v')
*              - const uint8_t *a: pointer to the input serialized ciphertext ('ct')

C1, C2 <- ct 

ct_u <- Decompress_Q( Decode_d_u( C1 ), d_u ) // polyvec_decompress
ct_v <- Decompress_Q( Decode_d_v( C2 ), d_v ) // poly_decompress

ct_u is param 'b'
ct_v is param 'v'

**************************************************/
static void unpack_ct(polyvec *res1, poly *res2, const uint8_t a[CIPHERTEXTBYTES]) {
    polyvec_decompress(res1, a);
    poly_decompress(res2, a + POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform 
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod Q
*
* Arguments:   - int16_t *res:          pointer to output buffer
*                                     (uniform mod q)
*              - const uint8_t *buf:  pointer to input buffer
*                                     (assumed to be uniform random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)

Function : Parse 
See more detail in Document 'KYBer'

**************************************************/
static unsigned int rej_uniform(int16_t *res, uint8_t *buf, unsigned int buflen) {

    unsigned int cnt, pos;
    uint16_t val0, val1;

    cnt = pos = 0;

    while ( cnt < N && pos + 3 <= buflen ) {
        val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0x11f;
        val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xfff;
        pos  += 3;

        if ( val0 < Q )
            res[cnt++] = val0;
        if ( cnt < N  && val1 < Q )
            res[cnt++] = val1;
    }

    return cnt;
}

/*************************************************
* Name:        samplematrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *res:        pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed ('publicseed')
*              - int transposed:      boolean deciding whether A or A^T
*                                     is generated

See more detail in Document 'KYBer'

**************************************************/
#define GEN_MATRIX_NBLOCKS ((12*N/8*(1 << 12)/Q + SHAKE128_RATE)/SHAKE128_RATE)
void samplematrix(polyvec *res, const uint8_t seed[SYMBYTES], int transposed) {
    
    unsigned int cnt;
    unsigned int buflen;
    uint8_t buf[GEN_MATRIX_NBLOCKS * SHAKE128_RATE + 2]; // Fixed : 506 (Depended on N, Q)
    keccak_state state;

    for (int i=0; i<K; i++) {
        for (int j=0; j<K; j++) {
            if (transposed)
                kyber_shake128_absorb(&state, seed, i, j);
            else
                kyber_shake128_absorb(&state, seed, j, i);
            
            shake128_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
            buflen = GEN_MATRIX_NBLOCKS*SHAKE128_RATE; // 504
            cnt = rej_uniform(res[i].vec[j].coeffs, buf, buflen);
            if ( cnt < N )
                printf(" Fatal error : error to generate public matrix 'A' ! \n ");
            // else
                // printf(" Success to create 'poly' which is an element of public matrix 'A' \n");
        }
    }
}

/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length SECRETKEYBYTES bytes)

---
Function : Kyber.KeyGen

Input  : -
Output : Key pair (pk, sk)

Generate random 32-bytes array 'publicseed' and 'noiseseed'

Sample matrix 'A' which element is sampled in uniform distribution
Sample vector 's' and 'e' use gaussian distribution

b <- As + e

pk <- ( Encode_12( b ) || publicseed )
sk <- Encode_12( s )

return pk, sk
---
**************************************************/
void keygen(uint8_t pk[PUBLICKEYBYTES], uint8_t sk[SECRETKEYBYTES]) {
    uint8_t nonce = 0;

    uint8_t buf[2*SYMBYTES];
    const uint8_t *publicseed = buf;
    const uint8_t *noiseseed  = buf + SYMBYTES;
    polyvec A[K], b, e, s;
    memset(&b, 0, sizeof(polyvec));

    randombytes(buf, SYMBYTES);                 // AES-256 CTR

    sha3_512(buf, buf, SYMBYTES);               // SHA3-512 (hash_g)

    samplematrix(A, publicseed, 0);             // SHAKE-128 and rejection sampling (Parse)

    samplevector_eta1(&s, noiseseed, &nonce);   // cbd, SHAKE-256 (prf), nonce = 0 1

    samplevector_eta1(&e, noiseseed, &nonce);   // cbd, SHAKE-256 (prf), nonce = 2 3

    for ( unsigned int i=0; i<K; i++ )          // b  <- A * s
        polyvec_mul(&b.vec[i], &A[i], &s);

    polyvec_add(&b, &b, &e);                    // b  <- A * s + e

    pack_pk(pk, &b, publicseed);                // pk <- ( Encode_12(b) || publicseed )
    pack_sk(sk, &s);                            // sk <- ( Encode_12(s) )

}

/*************************************************
* Name:        encrypt
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *ct:           pointer to output ciphertext
*                                      (of length CIPHERTEXTBYTES bytes)
*              - const uint8_t *m:     pointer to input message
*                                      (of length MSGBYTES bytes)
*              - const uint8_t *pk:    pointer to input public key
*                                      (of length PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins
*                                      used as seed (of length KYBER_SYMBYTES)
*                                      to deterministically generate all
*                                      randomness

---
Input  : 32-bytes message 'm'
Input  : 32-bytes random coin 'coins'
Input  : Public key 'pk' (Kyber512 : 768-bytes)
Output : Ciphertext 'ct' (Kyber512 : 768-bytes)


Generate random 32-bytes array 'publicseed' and 'noiseseed'

Sample matrix 'A' which element is sampled in uniform distribution
Sample vector 's' and 'e' use gaussian distribution

b <- As + e

pk <- ( Encode_12( b ) || publicseed )
sk <- Encode_12( s )

return pk, sk
---
**************************************************/
void encrypt(uint8_t ct[CIPHERTEXTBYTES], const uint8_t m[MSGBYTES], const uint8_t pk[PUBLICKEYBYTES], const uint8_t coins[SYMBYTES]) {
    uint8_t nonce = 0;
    polyvec AT[K], b, r, e1, ct_u; 
    poly e2, ct_v, mp;
    uint8_t publicseed[SYMBYTES];

    memset(&ct_u, 0, sizeof(ct_u)); memset(&ct_v, 0, sizeof(ct_v));

    unpack_pk(&b, publicseed, pk);              // b, publicseed <- ( Decode_12(pk) || publicseed )
    
    samplematrix(AT, publicseed, 1);            // A^T  <- Parse( XOF(publicseed) )

    samplevector_eta1(&r, publicseed, &nonce);  // r    <- cbd( prf(publicseed, nonce) )
    
    samplevector_eta2(&e1, coins, &nonce);      // e1   <- cbd( prf(publicseed, nonce) )
    
    poly_getnoise_eta1(&e2, coins, nonce);      // e    <- cbd( prf(publicseed, nonce) )
    
    for ( unsigned int i=0; i<K; i++ )          // ct_u <- A^T * r
        polyvec_mul(&ct_u.vec[i], &AT[i], &r);
    
    polyvec_add(&ct_u, &ct_u, &e1);             // ct_u <- A^T * r + e1

    polyvec_mul(&ct_v, &b, &r);                 // ct_v <- b^T * r
    poly_add(&ct_v, &ct_v, &e2);                // ct_v <- b^T * r + e2
    poly_frommsg(&mp, m);                       // mp   <- Decompress(m)
    poly_add(&ct_v, &ct_v, &mp);                // ct_v <- b^T * r + e2 + mp
    pack_ct(ct, &ct_u, &ct_v);                  // ct   <- ( Encode_d_u(Compress(ct_u, d_u)) || Encode_d_v(Compress(ct_v, d_v)) )
                                                // ct = ( C1 || C2 )
    
}

/*************************************************
* Name:        decrypt
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m:        pointer to output decrypted message
*                                   (of length MSGBYTES)
*              - const uint8_t *ct: pointer to input ciphertext
*                                   (of length CIPHERTEXTBYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length SECRETKEYBYTES)

---
Input  : Ciphertext 'ct' (Kyber512 : 768-bytes)
Input  : Secret key 'sk' (Kyber512 : 768-bytes)
Output : 32-bytes recovered message 'recovered'


Get first ciphertext ct_u and second ciphertext ct_v using function Decode_d and Decompress 
    ct_u <- Decompress_Q( Decode_d_u(ct) )
    ct_v <- Decompress_Q( Decode_d_v(ct + POLYVECBYTES) )

Get secret key polynomial vector 's' using function Decode_12
    s <- Decode_12( sk )

Compute decrypt process with ciphertext and secret key
    and compress_Q( mp_with_error, 1 )
        recovered <- poly_tomsg( ct_v - s^T*ct_u )

return recovered
---
**************************************************/
void decrypt(uint8_t recovered[MSGBYTES], const uint8_t ct[CIPHERTEXTBYTES], const uint8_t sk[SECRETKEYBYTES]) {

    polyvec ct_u, s;
    poly ct_v, mp, sTu;
    memset(&ct_u, 0, sizeof(ct_u)); memset(&s, 0, sizeof(s)); memset(&ct_v, 0, sizeof(ct_v)); memset(&mp, 0, sizeof(mp)); memset(&sTu, 0, sizeof(sTu));

    unpack_ct(&ct_u, &ct_v, ct);                // ct_u, ct_v <- Decompress(Decode_d_u(C1), d_u), Decompress(Decode_d_v(C2), d_v)
    unpack_sk(&s, sk);                          // s <- Decode_12(sk)

    polyvec_mul(&sTu, &s, &ct_u);
    poly_sub(&mp, &ct_v, &sTu);

    poly_tomsg(recovered, &mp);
}


int main() {

    /* Kyber512 keygen part */
    uint8_t sk[SECRETKEYBYTES], pk[PUBLICKEYBYTES];
    keygen(pk, sk);

    /* Treating message */
    uint8_t m[MSGBYTES];
    gettestmessage(m);

    /* Kyber512 encryption part */
    uint8_t ct[CIPHERTEXTBYTES], coins[SYMBYTES];
    randombytes(coins, SYMBYTES);
    encrypt(ct, m, pk, coins);

    /* Kyber512 decryption part */
    uint8_t recovered[MSGBYTES];
    decrypt(recovered, ct, sk);

    /* Print */
#if MESSAGETYPE == TESTVECTOR
    printarray(" message   = ", m, MSGBYTES);
    printstate(" recovered = ", recovered);
#elif MESSAGETYPE == STDIN
    printf(" * recovered = %s\n", recovered);
#endif

    return 0;
}

/*************************************************
* Name:        gettestmessage
*
* Description: Get const and standard input message
*
* Arguments:   - uint8_t m[MSGBYTES]:  pointer to output message
*                                   (of length MSGBYTES)

**************************************************/
void gettestmessage(uint8_t m[MSGBYTES]) {
#if MESSAGETYPE == TESTVECTOR
    // const uint8_t m[MSGBYTES] = { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};
    for ( unsigned int i=0; i<MSGBYTES; i++ )
        m[i] = i;
#else
    printf(" * Input your message(less than 32-bytes) > ");
    gets((char*) m);
    unsigned int messagelen = strlen((char*) m);
    if ( messagelen>32 ) {
        printf(" * Message space must be in 32-bytes! ");
        exit(0);
    }
    else if ( messagelen < 32 )
        m[messagelen] = 0;
#endif 
}

// EOF