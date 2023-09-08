/* Polynomial multiplications over Ring (Yongbhin Kim) */


#include "poly_mul.h"

/*

    poly c = (poly a) * (poly b) mod Q
    schoolbook multiplication : O(N^2)
*/
void schoolbook_mul(poly* res, const poly* a, const poly* b) {
    uint16_t c = KYBER_N;
    uint16_t q = KYBER_Q;
    uint32_t tempcoef = 0;
    uint16_t temppoly[KYBER_N * 2] = { 0, };

    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++) {
            tempcoef = a->coef[i] * b->coef[j];
            tempcoef = tempcoef + temppoly[i + j];
            temppoly[i + j] = tempcoef % KYBER_Q;
        }
    }

    for (int i=KYBER_N; i<KYBER_N*2; i++) { // 0 ~ 4
        res->coef[i-KYBER_N] = (KYBER_Q + temppoly[i-KYBER_N] - temppoly[i]) % KYBER_Q;
    }

}

// EOF