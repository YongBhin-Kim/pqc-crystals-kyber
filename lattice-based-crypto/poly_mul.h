/* Polynomial multiplications over Ring (Yongbhin Kim) */

#ifndef _POLY_MUL_
#define _POLY_MUL_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>


#define KYBER_N 256
#define KYBER_Q 3329


typedef struct {
    uint16_t coef[KYBER_N]; // R_Q의 원소 (coefficient of poly)
} poly;


// Declare Function
void schoolbook_mul(poly* res, const poly* a, const poly* b);

#endif

// EOF