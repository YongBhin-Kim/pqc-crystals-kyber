/* 1. LWE example (Yongbhin Kim) */


// order : lwe_example -> polymul -> polyvec_mul_example -> module_lwe_example


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define LWE

/* Define LWE-example cryptosystem parameters */
#define N 20        // columnsize of matrix
#define Q 79        // q which is used in Z_q
#define E 2         // Determine small error vector (range: 1~2)

// #define N 20     // 예시용
// #define Q 119    // 예시용
// #define E 3      // 예시용 small error vector 확장 --> Q의 범위 재적용

#define SUCCESS  1  // recover success
#define FAIL    -1  // recover failed


int lwe() {
    /*************************/
    /* Keypair generate part */
    /*************************/

    srand(time(NULL));

    uint8_t public_A[N]   = { 0x00, };      // nx1 maxrix
    uint8_t public_b[N]   = { 0x00, };      // nx1 vector
    uint8_t secret_err[N] = { 0x00, };      // nx1 vector
    uint8_t secret_key;                     // 1x1 vector

    // Generate Param : [public_A] and [secret_err] and [secret_key]
    for (int i=0; i<N; i++) {
        public_A[i] = rand() % Q;           // Range of public_A is [0, 78]
        secret_err[i] = rand() % E +1;      // Range of secret_err is [1, 2]
    }

    secret_key = rand() % E;                // Range of secret_key is [0, 78]
    uint8_t msg = rand() & 0x01;            // Let message in {0, 1}

    // Compute public b := As + e (public_A * secret_key + secret_err)
    for (int i=0; i<N; i++)
        public_b[i] = (public_A[i] * secret_key + secret_err[i]) % Q; 
    


    /*************************/
    /* Encryption part       */
    /*************************/

    // Compute ct_u := (sigma_public_A) mod Q 
    //                  r = ( 1 1 1 .. 1 ), e1 = 0
    //             ct_v := ((sigma_public_b) + Q/2 * msg) mod Q
    //                  r = ( 1 1 1 .. 1 ), e2 = 0
    uint8_t ct_u = 0x00;
    uint8_t ct_v = 0x00;
    for (int i=0; i<N; i++) {
        ct_u += public_A[i];    // rA
        ct_v += public_b[i];    // rb
        ct_u %= Q;
        ct_v %= Q;

    }
    ct_v += floor(Q / 2) * msg;
    ct_v %= Q;



    /*************************/
    /* Decryption part       */
    /*************************/
    uint8_t temp = 0x00;
    if ( ct_v > secret_key*ct_u ) {
        temp = ct_v - secret_key * ct_u;
    }
    else {
        temp = ct_v - secret_key * ct_u + Q;
        temp %= Q;
    }

    uint8_t recovered;
    if ( temp <= Q/2 )
        recovered = 0;
    else
        recovered = 1;

    /* print */
#if 1
    printf(" * [ public A (rand)             ] : ");
    for (int i=0; i<N; i++)
        printf("%3d ", public_A[i]);
    
    printf("\n * [ secret e (rand)           ] : ");
    for (int i=0; i<N; i++)
        printf("%3d ", secret_err[i]);
    
    printf("\n * [ secret k (rand)           ] : %3d", secret_key);
    
    printf("\n * [ public b (deterministic ) ] : ");
    for (int i=0; i<N; i++)
        printf("%3d ", public_b[i]);
    
    printf("\n----------------------------------------------------------------------------------------------------------");
    printf("\n * [ msg                       ] : %3d", msg);
    printf("\n * [ ct_u     (deterministic ) ] : %3d", ct_u);
    printf("\n * [ ct_v     (deterministic ) ] : %3d", ct_v);
    printf("\n * [ v-su     (deterministic ) ] : %3d", temp);
    printf("\n * [ recover  (deterministic ) ] : %3d\n", recovered);
    printf("----------------------------------------------------------------------------------------------------------\n");
#endif
    int ret = FAIL;
    if ( msg == recovered )
        ret = SUCCESS;
    
    return ret;
        
    
}

// main function
int main() {

    int ret;
    printf("==========================================================================================================\n");
    ret = lwe();
    if ( ret == SUCCESS )
        printf("* Success to recover!\n");
    else if ( ret == FAIL )
        printf("* Failed to recover..\n");
    else {
        printf("* Error\n");
        exit(1);
    }
    printf("==========================================================================================================\n");

    return 0;
}

// EOF