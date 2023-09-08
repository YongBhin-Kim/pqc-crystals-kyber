/*
    example for 
        poly, polyvec  compress/decompress
        poly, polyvec  tobytes, frombytes
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "poly_mul.h"

#define KYBER_K 2
// #define POLY_PRESS 1
// #define POLY_BYTES 1
// #define POLYVEC_PRESS  1
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320)

typedef struct {
    poly vec[3];     // 한 행의 원소의 개수는 k 일 것임 (2/3/4 at Kyber512/768/1024) // poly들을 모은 것
} polyvec;           // 원소가 R_3329로 이루어진 행렬의 한 행 (entry of A or b which is used in Crystals-Kyber as pk), primitive polynomial is <x^256+1>


#define KYBER_POLYCOMPRESSEDBYTES 128

#define KYBER_POLYBYTES 384
// #define KYBER_NAMESPACE(s) pqcrystals_kyber512_ref##s
// #define poly_csubq KYBER_NAMESPACE(_poly_csubq)
// void poly_csubq(poly *r);

void printbits(char *msg, uint16_t t);


/*************************************************
* Name:        poly_compress
*
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (of length KYBER_POLYCOMPRESSEDBYTES)
*              - poly *a:    pointer to input polynomial
**************************************************/
void poly_compress(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a)
{
    unsigned int i,j;
    uint8_t t[8];

    // poly_csubq(a);

    for(i=0;i<KYBER_N/8;i++) {
        for(j=0;j<8;j++)
            t[j] = ((((uint16_t)a->coef[8*i+j] << 4) + KYBER_Q/2)/KYBER_Q) & 15;

        r[0] = t[0] | (t[1] << 4);
        r[1] = t[2] | (t[3] << 4);
        r[2] = t[4] | (t[5] << 4);
        r[3] = t[6] | (t[7] << 4);
        r += 4;
        
    }
}

/*************************************************
* Name:        poly_decompress
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r:          pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of length KYBER_POLYCOMPRESSEDBYTES bytes)
**************************************************/
void poly_decompress(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES])
{
    unsigned int i;

    for(i=0;i<KYBER_N/2;i++) {
        r->coef[2*i+0] = (((uint16_t)(a[0] & 15)*KYBER_Q) + 8) >> 4;
        r->coef[2*i+1] = (((uint16_t)(a[0] >> 4)*KYBER_Q) + 8) >> 4;
        a += 1;
    }
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYBYTES bytes)
*              - poly *a:    pointer to input polynomial
**************************************************/
void poly_tobytes(uint8_t r[KYBER_POLYBYTES], poly *a)
{
  unsigned int i;
  uint16_t t0, t1;


  for(i=0;i<KYBER_N/2;i++) {
    t0 = a->coef[2*i];
    t1 = a->coef[2*i+1];
    r[3*i+0] = (t0 >> 0);
    r[3*i+1] = (t0 >> 8) | (t1 << 4);
    r[3*i+2] = (t1 >> 4);
  }
}






/*************************************************
* Name:        polyvec_compress
*
* Description: Compress and serialize vector of polynomials
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYVECCOMPRESSEDBYTES)
*              - polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_compress(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES], polyvec *a)
{
  unsigned int i,j,k;

#if (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 320))
  uint16_t t[4];
  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_N/4;j++) {
      for(k=0;k<4;k++)
        t[k] = ((((uint32_t)a->vec[i].coef[4*j+k] << 10) + KYBER_Q/2) / KYBER_Q) & 0x3ff;

      r[0] = (t[0] >> 0);
      r[1] = (t[0] >> 8) | (t[1] << 2);
      r[2] = (t[1] >> 6) | (t[2] << 4);
      r[3] = (t[2] >> 4) | (t[3] << 6);
      r[4] = (t[3] >> 2);
      r += 5;
    }
  }
#endif
}


/*************************************************
* Name:        polyvec_decompress
*
* Description: De-serialize and decompress vector of polynomials;
*              approximate inverse of polyvec_compress
*
* Arguments:   - polyvec *r:       pointer to output vector of polynomials
*              - const uint8_t *a: pointer to input byte array
*                                  (of length KYBER_POLYVECCOMPRESSEDBYTES)
**************************************************/
void polyvec_decompress(polyvec *r,
                        const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES])
{
  unsigned int i,j,k;

#if (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 320))
  uint16_t t[4];
  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_N/4;j++) {
      t[0] = (a[0] >> 0) | ((uint16_t)a[1] << 8);
      t[1] = (a[1] >> 2) | ((uint16_t)a[2] << 6);
      t[2] = (a[2] >> 4) | ((uint16_t)a[3] << 4);
      t[3] = (a[3] >> 6) | ((uint16_t)a[4] << 2);
      a += 5;

      for(k=0;k<4;k++)
        r->vec[i].coef[4*j+k] = ((uint32_t)(t[k] & 0x3FF)*KYBER_Q + 512) >> 10;
    }
  }
#else
#error "KYBER_POLYVECCOMPRESSEDBYTES needs to be in {320*KYBER_K, 352*KYBER_K}"
#endif
}

// main function
int main() {

#ifdef POLY_PRESS
    /*
        Poly compress & decompress example
    */
    uint8_t r[KYBER_POLYCOMPRESSEDBYTES];
    poly a; 
    memset(&r, 0, sizeof(r));
    memset(&a, 0, sizeof(a));

    for (int i=0; i<KYBER_N; i++) { 
        if ( i%2 == 0 )
            a.coef[i] = i; 
        else 
            a.coef[i] = i*2; // a <- 0101... 01 (256 bits)
    }

    // pack ciphertext // r <- compress(a)
    poly_compress(r, &a);

    // unpack ciphertext // r <- decompress(a)
    poly a2; 
    memset(&a2, 0, sizeof(a2));
    poly_decompress(&a2, r);

    printf(" [ Compress polynomial ] r <- poly_compress(a) \n\n");
    printf("a = \n");
    for (int i=0; i<KYBER_N; i++) {
        printf(" %3d", a.coef[i]);
        if (i % 16 == 15)
            printf("\n");
    }
    printf("\ncompress(a) = r = \n");
    for (int i=0; i<KYBER_POLYCOMPRESSEDBYTES; i++) {
        printf(" %3d", r[i]);
        if (i % 16 == 15)
            printf("\n");
    }
    printf("\ndecompress(r) = a2 = \n");
    for (int i=0; i<KYBER_N; i++) {
        printf(" %3d", a2.coef[i]);
        if (i % 16 == 15)
            printf("\n");
    }
    printf("\n");
#endif
#ifdef POLY_BYTES
    uint8_t r2[KYBER_POLYBYTES];
    memset(&r2, 0, sizeof(r2));
    poly_tobytes(r2, &a2); // uint8 array <- poly
    printf("\npoly_tobytes(a2) = r2 = \n");
    for (int i=0; i<KYBER_POLYBYTES; i++) {
        printf(" %3d", r2[i]);
        if (i % 16 == 15)
            printf("\n");
    }
#endif
#ifdef POLYVEC_PRESS
    polyvec a3;
    uint8_t r3[KYBER_POLYVECCOMPRESSEDBYTES];
    memset(&a3, 0, sizeof(a3));
    memset(&r3, 0, sizeof(r3));

    polyvec_compress(r3, &a3);

    printf("%d\n", KYBER_Q/2);
#endif

    // uint16_t t = 0b0000110100000001; // Q   -> 0
    // uint16_t t = 0b0000110011111111; // Q-1 -> 11..1 (10 bits)
    // uint16_t t = 0b0000110011111110; // Q-2 -> 11..1 (10 bits)
    // uint16_t t = 0b0000110011111101; // Q-3 -> 11..1 (10 bits)

    // uint16_t t = 0b0000110011111100; // Q-4 -> 11...0 (10 bits)
    // uint16_t t = 0b0000110011111011; // Q-5 -> 11...0 (10 bits)
    // uint16_t t = 0b0000110011111010; // Q-6 -> 11...0 (10 bits)
    // uint16_t t = 0b0000110011111001; // Q-7 -> 11...0 (10 bits)
    // uint16_t t = 0b0000110011111000; // Q-8 -> 11..01 (10 bits)

    // 0~1 : 0
    // 2~4 : 1
    // 3~8 : 2
    // 9~ : 3

    uint16_t t = 0b0000000000000000; // 0 -> 
    uint16_t compress_t, decompress_t;
    int prob_success = 0; //
    for (int i=0; i<KYBER_Q; i++) {
        t = i;
        compress_t = ((((uint32_t)t << 10) + KYBER_Q/2)/ KYBER_Q)&0x3ff;
        decompress_t = ((uint32_t)(compress_t & 0x3FF)*KYBER_Q + 512) >> 10;
        printf(" (%4d)", i);
        printbits("t      ", t);
        printbits("compress(t)   ", compress_t);
        printbits("decompress(t) ", decompress_t);
        printf("\n");
        if ( t == decompress_t ) {
            prob_success++;
        }
    }
    printf("prob of ( t == decompress(compress(t)) ) : %d percent\n", prob_success*100/KYBER_Q);


    // uint16_t compress_t = ((((uint32_t)t << 10) + KYBER_Q/2)/ KYBER_Q)&0x3ff;
    // printbits("t             ", t);
    // printbits("compress(t)   ", compress_t);


    // uint16_t decompress_t = ((uint32_t)(t & 0x3FF)*KYBER_Q + 512) >> 10;
    // printbits("decompress(t) ", decompress_t);

    return 0;
}

void printbits(char *msg, uint16_t t) {
    printf(" %s = 0b ", msg);
    for (int i=0; i<16; i++) {
        printf("%d", (t>>(15-i)) & 0x0001);
    }
    printf("\n");
}

// EOF