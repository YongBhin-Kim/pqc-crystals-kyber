#include <stdint.h>
#include "params.h"
#include "reduce.h"


/*************************************************
* Name:        csubq
*
* Description: Modulo reduction ( Case : x > Q )
**************************************************/
int16_t csubq(int16_t x) {

    int16_t res;
    res = x - Q; 
    res += (res >> 15) & Q;
    return res;
}

/*************************************************
* Name:        simple_reduce
*
* Description: Modulo reduction ( Case : x >= Q or x < 0 )
**************************************************/
int16_t simple_reduce(int16_t x) {

    while ( x>>15 )
        x += Q;  // x <- x + Q
    while ( x >= Q ) 
        x -= Q;

    return x;
}


#include <stdio.h> 
/*************************************************
* Name:        simple_reduce32
*
* Description: Modulo reduction 
                input value is bigger than 16-bits (used in multiplication over ring)
**************************************************/
int16_t simple_reduce32(int32_t x) {
    while ( x>>31 )
        x += Q;  // x <- x + Q
    while ( x >= Q ) 
        x -= Q;

    return x;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod Q in {0,...,Q-1}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {0,...,Q-1} congruent to a modulo Q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
    int16_t t;
    const int16_t v = ((1U << 26) + Q/2)/Q; // 20159

    t  = (int32_t)v*a >> 26;
    
    printf("%d\n", t);
    t *= Q;
    return a - t;
}
