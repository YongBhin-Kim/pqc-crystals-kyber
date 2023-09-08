#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

typedef struct {
    int16_t coeffs[N];
} poly;

void poly_add(poly *res, poly *p1, poly *p2);
void poly_sub(poly *res, poly *p1, poly *p2);
void poly_csubq(poly *res);
void poly_reduce(poly *res);
void poly_mul(poly *res, poly *p1, poly *p2);
void poly_tobytes(uint8_t res[POLYBYTES], poly *p);
void poly_frombytes(poly *res, const uint8_t a[POLYBYTES]);
void poly_frommsg(poly *res, const uint8_t m[SYMBYTES]);
void poly_tomsg(uint8_t res[SYMBYTES], poly *p);
void poly_compress(uint8_t res[POLYCOMPRESSEDBYTES], poly *p);
void poly_decompress(poly *res, const uint8_t a[POLYCOMPRESSEDBYTES]);
void poly_getnoise_eta1(poly *res, const uint8_t seed[SYMBYTES], uint8_t nonce);
void poly_getnoise_eta2(poly *res, const uint8_t seed[SYMBYTES], uint8_t nonce);
#endif