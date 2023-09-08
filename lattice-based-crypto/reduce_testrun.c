
#include <stdio.h>
#include <stdlib.h>

#define KYBER_Q 3329
#define KYBER_N 256
#define KYBER_K 2

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q in {0,...,q}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {0,...,q} congruent to a modulo q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
  int16_t t;
  const int16_t v = ((1U << 26) + KYBER_Q/2)/KYBER_Q;

  t  = (int32_t)v*a >> 26;
  t *= KYBER_Q;
  return a - t;
} 

int16_t simple_reduce(int16_t a) {

  if ( a>>15 ) // If a < 0
    return (a + KYBER_Q); // a <- a % Q
  else
    return a; // a <- a % Q
    
}

int main() {
    int16_t a = -1000; // 1000 1111 1111
    int16_t res1, res2;
    printf("%d mod Q = %d\n", a, a % 3329);

    res1 = simple_reduce(a);
    res2 = barrett_reduce(a);
    printf("%d mod Q = %d\n", a, (a % KYBER_Q));
    printf("%d mod Q = %d\n", a, res1);
    printf("%d mod Q = %d\n", a, res2);

}