//
// Function for print byte array, 'poly', 'polyvec' and 'polymat'
//

#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

/* If you want to use an already defined test vector as a message, define this value(MESSAGETYPE) as TESTVECTOR, 
    and if you want a message as standard input, define this value(MESSAGETYPE) as STDIN*/
#define TESTVECTOR 0
#define STDIN 1
#define MESSAGETYPE STDIN 

/* Debug options : If 1, debug mode operated (print the values ​​used and discarded inside the function) */
#define _DEBUG 1

/* If you want to check error tolerance gaps then modify below : */
#define CHECK_ERROR_RANGE 0

/* Determine print type : HEX or DECIMAL */
#define DECIMAL 10
#define HEX 16
#define PRINTTYPE DECIMAL // check if you want

void printarray(char *msg, const uint8_t *state, int bytelen);
void printpoly(char *msg, poly *p);
void printvector(char *msg, polyvec *pv);
void printmatrix(char *msg, polyvec pv[K]);
#endif