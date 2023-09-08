/* 2. Polynomial multiplications over Ring - example (Yongbhin Kim) */

/* 
    Polynomial  Quotient ring := Rq
*/
#include "poly_mul.h"


// main function
int main(void) {
    poly a, b, c;
    memset(&a, 0, sizeof(poly));
    memset(&b, 0, sizeof(poly));
    memset(&c, 0, sizeof(poly));

    for (int i = 0; i < KYBER_N; i++) {
        // a.coef[i] = i + 1;
        // b.coef[i] = i * 2;
        a.coef[i] = 1;       // 예시용  a := 1 + x + x^2 + x^3 + ... + x^255
        b.coef[i] = 1;       // 예시용  b := 1 + x + x^2 + x^3 + ... + x^255
    }


    schoolbook_mul(&c, &a, &b);

    printf(" [ Do polymul ] : c = a * b \n\n");
    for (int i = 0; i < KYBER_N; i++)
        printf("   c[%d] = %d\n", i, c.coef[i]);
    return 0;
}

// EOF