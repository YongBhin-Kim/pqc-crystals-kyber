#ifndef MODULE_RLWE_EXAMPLE
#define MODULE_RLWE_EXAMPLE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
// #include "poly_mul.h"

#define K 2
#define Q 3329
#define N 256
#define ETA1 3

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

typedef struct {
    int16_t coef[N]; // R_Q의 원소 (coefficient of poly)
} poly;

typedef struct {
    poly vec[K];     // 한 행의 원소의 개수는 k 일 것임 (2/3/4 at Kyber512/768/1024) // poly들을 모은 것
} polyvec;           // 원소가 R_3329로 이루어진 행렬의 한 행 (element of A or b which is used in Crystals-Kyber as pk), primitive polynomial is <x^256+1>


/* 다항식 몫환 위에서의 연산 정의 */
int16_t barrett_reduce(int16_t a); // (23.2.1 추가)
int16_t csubq(int16_t a); // Conditionally subtract Q (23.2.1 추가)
void poly_csubq(poly *res, poly *a); // Conditionally subtract Q over R_Q(23.2.1 추가)
void poly_add(poly *res, poly *a, const poly *b); // 덧셈
void poly_sub(poly *res, poly *a, const poly *b); // 뺄셈1
void poly_sub2(poly *res, poly *a, const poly *b); // 뺄셈2
void polyvec_add(polyvec *res, polyvec *a, const polyvec *b); // 벡터 덧셈
void poly_mul(poly* res, const poly* a, const poly* b); // 곱셈(O(N^2))
void polyvec_mul(poly *res, const polyvec *a, const polyvec *b); // 벡터 곱셈
void mat_polyvec_mul(polyvec *res, polyvec a[K], polyvec *b); // 행렬 * 벡터 연산
void trans_matrix(polyvec *res, polyvec *a); // 행렬 전치

/* 샘플링과 관련된 함수 */
static uint32_t load32_littleendian(const uint8_t x[4]); // 리틀엔디안 변환
static void cbd(poly *res, const uint8_t buf[ETA1*N/4]); // 중심이항분포
void prf(uint8_t res[ETA1*N/4], uint8_t seed[32], uint8_t nonce); // 가짜난수생성
void samplevector(polyvec *res, uint8_t seed[32], uint8_t *nonce); // 벡터 샘플링
void toy_xof(uint8_t res[512], uint8_t seed[32]); // 확장함수 1
void xof(uint8_t res[512], uint8_t seed[32]); // 확장함수 2
// void parse()
void samplematrix(polyvec *res, uint8_t seed[32]); // 행렬 샘플링 1
void updateseed(uint8_t res[32], uint8_t seed[32]); // samplematrix2를 위해 seed를 갱신하는 함수
void samplematrix2(polyvec *res, uint8_t seed[32]); // 행렬 샘플링 2
void get_randomarray(uint8_t res[32]); // 랜덤 배열 생성

/* 메시지 다루는 함수 */
void poly_frommsg(poly *res, const uint8_t m[32]); // 메시지를 다항식 몫환의 원소로 변환
void poly_tomsg(uint8_t res[32], poly *mp); // 다항식 몫환의 원소를 메시지로 변환

/* Module-RLWE Scheme 함수 */
// void keygen(polyvec *b, uint8_t publicseed[32], polyvec *s); // 키 생성 함수
// void enc(polyvec *ct_u, poly* ct_v, uint8_t *m, polyvec *b, uint8_t publicseed[32], uint8_t coin[32]); // 암호화 함수
// void dec(polyvec *ct_u, poly *ct_v, polyvec *s, uint8_t recovered[32]); // 복호화 함수
polyvec keygen(polyvec *b, uint8_t publicseed[32], polyvec *s); // + return e
polyvec enc(polyvec *ct_u, poly* ct_v, uint8_t *m, polyvec *b, uint8_t publicseed[32], uint8_t coin[32]); // + return r
void dec(polyvec *ct_u, poly *ct_v, polyvec *s, uint8_t recovered[32], polyvec *e, polyvec *r);

void module_rlwe_test(); // Module-RLWE 테스트 함수

/* 출력 함수 */
void printstate(char *msg, uint8_t state[32]); // 배열 출력 함수
void printpoly(char *msg, poly *polynomial); // 다항식 몫환의 원소 출력 함수
void printvector(char *msg, polyvec *vector); // 벡터 출력 함수
void printmatrix(char *msg, polyvec vector[K]); // 행렬 출력 함수
#endif