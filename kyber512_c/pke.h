#ifndef PKE_H
#define PKE_H

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

static void pack_pk(uint8_t res[PUBLICKEYBYTES], polyvec *p, const uint8_t seed[SYMBYTES]);
static void pack_sk(uint8_t res[SECRETKEYBYTES], polyvec *p);
static void unpack_pk(polyvec *res1, uint8_t res2[SYMBYTES], const uint8_t a[PUBLICKEYBYTES]);
static void unpack_sk(polyvec *res, const uint8_t a[SECRETKEYBYTES]);
static void pack_ct(uint8_t res[CIPHERTEXTBYTES], polyvec *pv, poly *p);
static void unpack_ct(polyvec *res1, poly *res2, const uint8_t a[CIPHERTEXTBYTES]);
static unsigned int rej_uniform(int16_t *res, uint8_t *buf, unsigned int buflen);
void samplematrix(polyvec *res, const uint8_t seed[SYMBYTES], int transposed);

void keygen(uint8_t pk[PUBLICKEYBYTES], uint8_t sk[SECRETKEYBYTES]);
void enc(uint8_t ct[CIPHERTEXTBYTES], const uint8_t m[MSGBYTES], const uint8_t pk[PUBLICKEYBYTES], const uint8_t coins[SYMBYTES]);
void dec(uint8_t recovered[MSGBYTES], const uint8_t ct[CIPHERTEXTBYTES], const uint8_t sk[SECRETKEYBYTES]);

void gettestmessage(uint8_t m[MSGBYTES]);
#endif