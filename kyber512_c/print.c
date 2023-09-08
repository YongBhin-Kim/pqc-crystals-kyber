//
// Function for print byte array, 'poly', 'polyvec' and 'polymat'
//


#include "print.h"


/* Print byte array */
void printarray(char *msg, const uint8_t *state, int bytelen) {
    printf("    * %s \n", msg);
    printf("     [ ");

    if ( PRINTTYPE == DECIMAL ) {
        for (int i=0; i<bytelen; i++) {

            printf("%3d ", state[i]);
            if ( (bytelen > 32) && ((i%32) == 31) ) 
                printf("\n       ");
        }
    }
    else if ( PRINTTYPE == HEX ) {
        for (int i=0; i<bytelen; i++) {

            printf("%2x ", state[i]);
            if ( (bytelen > 32) && ((i%32) == 31) ) 
                printf("\n       ");
        }
    }
    printf(" ] \n");
    printf("\n");
}

/* Print polynomial */
void printpoly(char *msg, poly *p) {

    printf("    * %s \n", msg);
    printf("     [ ");
    if ( PRINTTYPE == DECIMAL ) {
        for (int j=0; j<N; j++) {
            printf("%4d ", p->coeffs[j]);
            if ( j%16==15 )
                printf("\n        ");
        }
    }
    else if ( PRINTTYPE == HEX ) {
        for (int j=0; j<N; j++) {
            printf("%3x ", p->coeffs[j]);
            if ( j%16==15 )
                printf("\n        ");
        }
    }

    printf(" ] \n");
    printf("\n");
}

/* Print polynomial vector */
void printvector(char *msg, polyvec *pv) {

    printf("    * %s \n", msg);
    if ( PRINTTYPE == DECIMAL ) {
        for (int i=0; i<K; i++) {
            printf("     [ ");
            for (int j=0; j<N; j++) {
                printf("%4d ", pv->vec[i].coeffs[j]);
                if ( j%16==15 )
                    printf("\n        ");
            }
            printf(" ] \n");
        }
    }
    else if ( PRINTTYPE == HEX ) {
        for (int i=0; i<K; i++) {
            printf("     [ ");
            for (int j=0; j<N; j++) {
                printf("%3x ", pv->vec[i].coeffs[j]);
                if ( j%16==15 )
                    printf("\n        ");
            }
            printf(" ] \n");
        }
    }
    printf("\n");
}

/* Print polynomial vector array(K x K matrix) */
void printmatrix(char *msg, polyvec pv[K]) {

    printf("    * %s \n", msg);



    printf("     [ ");
    if ( PRINTTYPE == DECIMAL ) {
        for (int k=0; k<K; k++ ) {
            for (int i=0; i<K; i++) {
                printf("  [ ");
                for (int j=0; j<N; j++) { 
                    printf("%4d ", pv[k].vec[i].coeffs[j]);
                    if ( j%16==15 )
                        printf("\n        ");
                }
                printf(" ] \n      ");
            }
        }
    }
    else if ( PRINTTYPE == HEX ) {
        for (int k=0; k<K; k++ ) {
            for (int i=0; i<K; i++) {
                printf(" [ ");
                for (int j=0; j<N; j++) { 
                    printf("%3x ", pv[k].vec[i].coeffs[j]);
                    if ( j%16==15 )
                        printf("\n        ");
                }
                printf(" ] \n");
            }
        }
    }

    printf("  ] \n");
    printf("\n");
}