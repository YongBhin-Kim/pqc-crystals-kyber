#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

typedef struct {
    poly vec[K];
} polyvec;

void polyvec_add(polyvec *res, polyvec *pv1, polyvec *pv2);
void polyvec_csubq(polyvec *res);
void polyvec_reduce(polyvec *res);
void polyvec_mul(poly *res, polyvec *pv1, polyvec *pv2);
void polyvec_tobytes(uint8_t res[POLYVECBYTES], polyvec *pv);
void polyvec_frombytes(polyvec *res, const uint8_t a[POLYVECBYTES]);
void polyvec_compress(uint8_t res[POLYVECCOMPRESSEDBYTES], polyvec *pv);
void polyvec_decompress(polyvec *res, const uint8_t a[POLYVECCOMPRESSEDBYTES]);
void samplevector_eta1(polyvec *res, const uint8_t seed[SYMBYTES], uint8_t *nonce);
void samplevector_eta2(polyvec *res, const uint8_t seed[SYMBYTES], uint8_t *nonce);
#endif