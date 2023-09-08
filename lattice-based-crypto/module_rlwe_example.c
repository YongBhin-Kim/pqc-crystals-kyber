/* Module-RLWE(Simple Kyber512) example (Yongbhin Kim) 

    1. All functions output the argument assigned to the res parameter as an address. In other words, res is an output that does not affect input

    2. 
        Let A in R_Q^{K*K}
        Let b in R_Q^K
        Tr(A) := Transpose matrix of A
        Tr(b) := Transpose vector of b
*/

/*
[23.1.27 세미나 이후 수정사항 + 23.1.29 세미나]
    매크로 가시성 수정 (modify visibility of macro)

    samplematrix -> samplematrix2 (XOF function Implementing...)
    xof함수 구현 중.. shake256 끌어다 쓸지?
    주석 보완 (전부 영문)
    message m, poly mp 변환 과정과 그 의미 설명 (Q/2 +- Q/4) 그리고 중심이항분포(cbd)를 쓰는 이유
*/

/*
[23.1.29 세미나 이후 수정사항]
    문서의 벡터 표기를 강조로 변환
*/

#include "module_rlwe_example.h"

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q in {0,...,q-1}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {0,...,q-1} congruent to a modulo q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
    int16_t t;
    const int16_t v = ((1U << 26) + Q/2)/Q;

    t  = (int32_t)v*a >> 26;
    t *= Q;
    return a - t;
}

/*************************************************
* Name:        simple_reduce
*
* Description: simple reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q in {0,...,q-1}
*               (It means a <- a % Q)
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {0,...,q-1} congruent to a modulo q.
**************************************************/
int16_t simple_reduce(int16_t a) {

    if ( a>>15 ) // if a < 0
        a += Q;  // a <- a + Q 
    return a;
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies Simple or Barrett reduction to all coefficients of a polynomial
*
* Arguments:   - poly *a: pointer to input polynomial
*              - poly *res: pointer to output polynomial
**************************************************/
void poly_reduce(poly *res, poly *a) {
    unsigned int i;
    for(i=0;i<N;i++)
        res->coef[i] = simple_reduce(a->coef[i]);
}

/*************************************************
* Name:        polyvec_reduce
*
* Description: Applies Simple or Barrett reduction to all polynomial of a polyvec
*
* Arguments:   - poly *a: pointer to input polynomial vector
*              - poly *res: pointer to output polynomial vector
**************************************************/
void polyvec_reduce(polyvec *res, polyvec *a) {
    for (int i=0; i<K; i++)
        poly_reduce(&(res->vec[i]), &(a->vec[i]));
}

/* Result(polynomial) = polynomial + polynomial (polynomial in R_Q)*/
void poly_add(poly *res, poly *a, const poly *b) {
    poly ap, bp;
    poly_reduce(&ap, a); // (23.2.1. 감산 추가)
    poly_reduce(&bp, b);

    for (int i = 0; i < N; i++)
        res->coef[i] = (ap.coef[i] + bp.coef[i]) % Q;
}

/* Polynomial subtraction method 1
     Result(polynomial) = polynomial - polynomial (polynomial in R_Q) */
void poly_sub(poly *res, poly *a, const poly *b) {

    
    /* Method 1 */
    for (int i=0; i<N; i++)
        res->coef[i] = (a->coef[i] + Q - b->coef[i]) % Q;
    
    /* Method 2 */
    // poly temp;
    // for (int i=0; i<N; i++) {
        // temp.coef[i] = a->coef[i] - b->coef[i];
    // }
    // poly_csubq(res, &temp);
}

/* Polynomial subtraction method 2 
     Result(polynomial) = polynomial - polynomial (polynomial in R_Q) */
void poly_sub2(poly *res, poly *a, const poly *b) {
    poly ap, bp;
    poly_csubq(&ap, a);
    poly_csubq(&bp, b);

    for (int i=0; i<N; i++) {
        if (ap.coef[i] >= bp.coef[i])
            res->coef[i] = (ap.coef[i] - bp.coef[i]);     // Q isn't need becapuse coef of a is less than Q and coef of b is less than coef of a
        else
            res->coef[i] = (ap.coef[i] + Q - bp.coef[i]); // a->coef < b->coef < Q
    }

}

/* Result(Polynomial vector) = vector + vector (vector in R_Q^K)*/
void polyvec_add(polyvec *res, polyvec *a, const polyvec *b) {
    
    for (int i=0; i<K; i++)
        poly_add(&res->vec[i], &a->vec[i], &b->vec[i]);
}

/*
    Result (polynomial) = polynomial * polynomial (polynomial in R_Q)
        poly_mul uses schoolbook multiplication algorithm: O(N^2)
*/
void poly_mul(poly* res, const poly* a, const poly* b) {
    int32_t tempcoef = 0;
    int16_t temppoly[N * 2] = { 0, };
    poly ap, bp;

    poly_csubq(&ap, a);
    poly_csubq(&bp, b);

    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            tempcoef = ap.coef[i] * bp.coef[j];
            tempcoef = tempcoef + temppoly[i + j]; // Required because i + j appears multiple times
            temppoly[i + j] = tempcoef % Q;
        }
    }


    /* Reduction method 1 */
    // for (int i=N; i<N*2; i++) { 
    //     res->coef[i-N] = (Q + temppoly[i-N] - temppoly[i]) % Q; // 수정완료 -> temppoly가 구조체로 선언하지 않아서 sub함수 이용 못 함, mod Q는 반드시 필요함
    // }

    /* Reduction method 2 */
    for (int i=N; i<N*2; i++) { 
        if ( temppoly[i-N] >= temppoly[i] )
            res->coef[i-N] = temppoly[i-N] - temppoly[i];
        else
            res->coef[i-N] = temppoly[i-N] + Q - temppoly[i]; // temppoly[i-N] + Q < 2Q => 13비트로 표현 가능, 즉 16비트 자료형 안에 들어오므로 문제 x
    }



}

/* Result(polynomial) = Tr(vector) x vector (vector in R_Q^K)
    (It means rowvector x columnvector) 
*/
void polyvec_mul(poly *res, const polyvec *a, const polyvec *b) { // It is a column, but it is used as a row, i.e., a transposed vector is used as input

    poly temppoly;

    for (int i=0; i<K; i++) {
        // memset(&temppoly, 0, sizeof(poly)); // => Not needed. because in poly_mul, res doesn't affect the input
        poly_mul(&temppoly, &(a->vec[i]), &(b->vec[i]));
        poly_add(res, res, &temppoly); // res <- res + temppoly
    }
}

/* Result(Vector) = Matrix x ColumnVector */
void mat_polyvec_mul(polyvec *res, polyvec a[K], polyvec *b) { // It is a column but contains row values, that is, a transposed matrix is ​​used as input
    for (int i=0; i<K; i++)
        polyvec_mul(&(res->vec[i]), &a[i], b); // First poly(element) of res <- multiplication of first row of A and s

}

/* Tr(A) <- Transpose matrix of matrix A */
void trans_matrix(polyvec *res, polyvec *a) { 
    for (int i=0; i<K; i++ )
        for (int j=0; j<K; j++)
            for (int k=0; k<N; k++)
                res[i].vec[j].coef[k] = a[j].vec[i].coef[k];
} 


/*
    32-bit unsigned integer <- array in B^4
    (Little endian : x3x2x1x0 <- [x0][x1][x2][x3])
*/
static uint32_t load32_littleendian(const uint8_t x[4]) { 
    uint32_t res;
    res  = (uint32_t)x[0];
    res |= (uint32_t)x[1] << 8;
    res |= (uint32_t)x[2] << 16;
    res |= (uint32_t)x[3] << 24;
    return res;
}

/* Sample secret vector s, e pseudorandomly 
    Crystals-Kyber uses CBD_ETA(PRF(noiseseed, nonce))
     -> PRF : SHAKE256 (중심이항분포에서 비밀키 s와 에러벡터 e를 추출)
*/
static void cbd(poly *res, const uint8_t buf[ETA1*N/4]) { // (buf is random array) <- buf0123 b14b12b10... + b15
    unsigned int i,j;
    uint32_t t,d;
    int16_t a,b;


    for(i=0;i<N/8;i++) {
        t  = load32_littleendian(buf+4*i);
        d  = t & 0x55555555;      // t    & 0101 0101 0101 0101 .... // ..   t14 t12 t10 t8 t6 t4 t2 t0
        d += (t>>1) & 0x55555555; // t>>1 & 0101 0101 0101 0101 .... // .. + t15 t13 t11 t9 t7 t5 t3 t1

        for(j=0;j<8;j++) {
            a = (d >> (4*j+0)) & 0x3;
            b = (d >> (4*j+2)) & 0x3;
            res->coef[8*i+j] = a - b; // Kyber is drawn from the range of -ETA to ETA
            // res->coef[8*i+j] = (a+b)%(ETA1); // In this example, it is drawn from the range of 0 to ETA1-1.
        }
    }
}
static uint32_t load24_littleendian(const uint8_t x[3]) {
    uint32_t r;
    r  = (uint32_t)x[0];
    r |= (uint32_t)x[1] << 8;
    r |= (uint32_t)x[2] << 16;
    return r;
}
static void cbd3(poly *r, const uint8_t buf[3*N/4])
{
    unsigned int i,j;
    uint32_t t,d;
    int16_t a,b;

    for(i=0;i<N/4;i++) {
        t  = load24_littleendian(buf+3*i);
        d  = t & 0x00249249;
        d += (t>>1) & 0x00249249;
        d += (t>>2) & 0x00249249;

        for(j=0;j<4;j++) {
            a = (d >> (6*j+0)) & 0x7;
            b = (d >> (6*j+3)) & 0x7;
            r->coef[4*i+j] = a - b;
        }
    }
}


/*
Pseudo Random Function (뻥튀기)
Compute res[192] using seed[32] and nonce
*/
void prf(uint8_t res[ETA1*N/4], uint8_t seed[32], uint8_t nonce) {

    int count = 0;
    for (int i=0; i<ETA1*N/4; i++) { // 
        count++;
        if ( (count % ((ETA1*N/4)/6)) == 0 ) // 6번마다 // 0번 -> 0번, ..., 192번 -> 32번
            srand((seed[(uint8_t)(count/6)] + nonce)); // rand(seed[0~31] and nonce) -> 1번째에 seed[0]+nonce, 7번째에 seed[1]+nonce, .. 187번째에 seed[31]+nonce 을 seed로 사용
        res[i] = (uint8_t)(rand()%256); // 8비트 자료형 안에 들어가야함
    }
}
void samplevector(polyvec *res, uint8_t seed[32], uint8_t *nonce) {
    // 중심이항분포에서 s와 e의 원소를 샘플링
    uint8_t buf[ETA1*N/4];
    for (int i=0; i<K; i++) {
        prf(buf, seed, (*nonce));
        cbd3(res->vec+i, buf);
        (*nonce)++; // nonce의 주소를 받아서 이 함수에서 증가하는 것으로 수정 (23.1.26) s : 0, 1 e : 2, 3 r : 0, 1
    }

#if CHECK_ERROR_RANGE == 1
    for (int i=0; i<K; i++) {
        for (int j=0; j<N; j++)
            res->vec[i].coef[j] = 2;
    }
#endif
}


/* XOF function 
    23.1.30 (Comment)
*/
void toy_xof(uint8_t res[512], uint8_t seed[32]) { 

    uint16_t i = 0;
    uint8_t k = 0;
    int j = 0;
    uint8_t temp;

    while (i < 512) { // 0 ~ 512
        temp = 0x00;
        for (j = 0; j < 8; j++) {
            temp ^= seed[(k+j)%32];
        }
        res[i] %= 128;
        for (j = 0; j < 8; j++) {
            res[i] ^= seed[(k+j)%32];
        }
        res[i] <<= 3;
        res[i] ^= (temp>>3);

        k+=15;
        k%=32;
        i++;
    }
    // printstate("res = ", res+256);
}

/* Extended function
    Generate determenistic 512-byte array using random 32-byte seed array
    (Used in generating matrix)
*/
void xof(uint8_t res[512], uint8_t seed[32]) {

    for (int i=0; i<512; i++) { // srand(seed[0]), srand(seed[1]), srand(seed[2]), ...,  srand(seed[31])
        if ((i % 16) == 0)  
            srand(seed[(uint8_t)(i/32)]);
        res[i] = (uint8_t)(rand()%256);
    }
}


/* Sample public key matrix A pseudorandomly
    Crystals-Kyber uses Parse(XOF(publicseed))
     -> parse : rej_uniform (균등 분포에서 무작위로 행렬 A의 원소를 추출)

*/
void samplematrix(polyvec *res, uint8_t seed[32]) { 
    int count = 0;
    for (int i=0; i<K; i++) { // 2
        for (int j=0; j<K; j++) {  // 2
            for (int k=0; k<N; k++) { // 256 --> 1024번
                count++;
                if ( (count % 32) == 1 )
                    srand(seed[(uint8_t)(count/32)]); // rand(seed[0~31]) -> 1번째에 seed[0], 32번째에 seed[1], .. 992번째에 seed[31]을 seed로 사용
                res[i].vec[j].coef[k] = (int16_t)(rand()%Q); // rand에 seed를 이용하여 균등 분포로 샘플링
            }
        }
    }

}

/* Update seed used in 'samplematrix2'

*/
void updateseed(uint8_t res[32], uint8_t seed[32]) {
    // for (int i=0; i<32; i++) {
    //     res[i] = seed[i] ^ seed[(i+1)%32];
    // }
    int i = 0;
    int k = 0;
    int j = 0;
    while (i < 32) { // 0 ~ 31
        for (j = 0; j < 11; j++) {
            res[i] ^= seed[(k+j)%32];
        }
        k+=7;
        k%=32;
        i++;
    }
}

/* Sample public key matrix A pseudorandomly
    Crystals-Kyber uses Parse(XOF(publicseed))
     -> parse : rej_uniform (균등 분포에서 무작위로 행렬 A의 원소를 추출)
*/
void samplematrix2(polyvec *res, uint8_t seed[32]) {
    uint8_t buf[512] = { 0x00, };
    uint16_t l = 0;
    uint8_t upseed[32] = { 0x00, };
    memcpy(upseed, seed, 32);

    for (int i=0; i<K; i++) {
        for (int j=0; j<K; j++) {
            memset(&buf, 0, sizeof(buf));
            xof(buf, upseed);
            l = 0;
            for (int k=0; k<N; k++) {
                res[i].vec[j].coef[k] = ((((uint16_t)buf[l])<<7) ^ ((uint16_t)(buf[l+1])&0x00ff)) % Q;
                res[i].vec[j].coef[k] ^= k;
                l+=2;
            }
            updateseed(upseed, upseed);
        }
    }
}


/* Get random 32-byte array (B^32) */
void get_randomarray(uint8_t res[32]) {
    for (int i=0; i<32; i++) {
        res[i] = (uint8_t) (rand() % 256);  // entry in 8-bit integer (max : 256)
    }
}

/* polynomial <- message(32-byte array) 
    if message bit is 0   
        then coef value is 0
    elif message bit is 1
        then coef value is 1665 = ceil(Q/2)
*/
void poly_frommsg(poly *res, const uint8_t m[32]) {
    unsigned int i,j;
    int16_t mask;

    for(i=0; i<N/8; i++) {
        for(j=0; j<8; j++) {
            mask = -(int16_t)((m[i] >> j)&1); // - is bit-reversing
            res->coef[8*i+j] = mask & ((Q+1)/2); // 0 or 1665 == ceil((messagebit * Q) / 2)
        }
    }
}

/* message(32-byte array) <- polynomial(2*N-byte structor) // 256비트 <- 256개의 계수
*/
void poly_tomsg(uint8_t res[32], poly *mp) {
    unsigned int i,j;
    int16_t t;

    for(i=0;i<N/8;i++) {
        res[i] = 0;
        for(j=0;j<8;j++) {
            t = ((((int16_t)mp->coef[8*i+j] << 1) + Q/2)/Q) & 1; 
            res[i] |= t << j;
        }
    }
}


/* Module-RLWE.KeyGen
 * Input : -
 * Output : Public key b in R_Q^K
 * Output : Public seed publicseed in B^32 for 'A'
 * Output : Secret key s in R_Q^K

 * In Kyber :
        pk : b || publicseed (768 + 32 bytes == 800 bytes)
        sk : s (768 bytes)
 */
polyvec keygen(polyvec *b, uint8_t publicseed[32], polyvec *s) {
    polyvec A[K], e;
    uint8_t noiseseed[32];

    /* Sample publicseed, noiseseed with ETA */
    get_randomarray(publicseed); // pubseed in B^32
    get_randomarray(noiseseed);  // pubseed in B^32

#if _DEBUG == 1

    printf(" \n============================================================================================================\n");
    printf("                                                 [Debug mode]\n");
    printf("                         Print the values used and discarded inside the function (with print the message m)\n");
    printf("                                    ( noiseseed, r, mp(for enc), mp(for dec) )\n");
    printf(" \n------------------------------------------------------------------------------------------------------------\n");
    printstate("noiseseed (random 32-byte array) = ", noiseseed);
#endif

    /* Sample public key matrix A with publicseed */
    samplematrix2(A, publicseed);

    /* Sample secret key vector s and secret noise vector e with noiseseed(eta) and nonce */
    uint8_t nonce = 0;
    samplevector(s, noiseseed, &nonce); // nonce 0, 1, .. K-1
    samplevector(&e, noiseseed, &nonce); // nonce K, K+1, .. 2K-1

#if _DEBUG == 1 
    printvector("e (in {-ETA1 ~ ETA1}) = ", s);
#endif

    /* Compute public key b := A^T * s + e */
    mat_polyvec_mul(b, A, s);
    polyvec_add(b, b, &e);
    return e;

}


/* Module-RLWE.Enc
Input : Public key b in R_Q^K
Input : Public seed publicseed in B^32
Input : Message m in B^32
Input : Random coin in B^32
Output : Ciphertext c_u in R_Q^K
Output : Ciphertext c_v in R_Q
 */
polyvec enc(polyvec *ct_u, poly* ct_v, uint8_t *m, polyvec *b, uint8_t publicseed[32], uint8_t coin[32]) {
    polyvec A[K], r, AT[K], bT;
    poly mp;

    /* Sample r in R_ETA1^K (Kyber : Sample r in {-ETA1, ..., 0, 1, ..., ETA1}) 
        중심이항분포에서 추출
    */
    uint8_t nonce = 0;
    samplevector(&r, coin, &nonce); // 0, 1
#if _DEBUG == 1
    printvector("r (in R_ETA1^K) = ", &r);
#endif
    /* Compute sampled public key matrix 'A' with publicseed 
        (find matrix 'a' using publicseed which included in public key(균등분포에서 원소를 무작위로 추출한 행렬 A in R_Q^{K*K}))
    */
    samplematrix2(A, publicseed);

#if _DEBUG == 1 
    printmatrix("A (entry in R_Q^{K*K} which is uniformly random) = ", A);
#endif

    /* Compute ct_u = A^T * r (Kyber에서는 e1_ETA2 (in R_Q^K) 가 더해집니다)
        (ct_u in R_Q^K)
    */
    trans_matrix(AT, A);
    mat_polyvec_mul(ct_u, AT, &r); 


    /* Compute ct_v = b^T * r + mp (Kyber에서는 e2_ETA2 (in R_Q^K) 가 더해집니다)
        (ct_v in R_Q) 

        (polynomial mp in R_Q) <- (32 byte array m) */
    poly_frommsg(&mp, m);
    polyvec_mul(ct_v, b, &r);
    poly_add(ct_v, ct_v, &mp);
    
#if _DEBUG == 1
    printpoly("mp (Polynomial forms of message 'm' in R_Q) = ", &mp);
#endif
    return r;
}


/* Module-RLWE.Dec
Input : Ciphertext c_u in R_Q^K
Input : Ciphertext c_v in R_Q
Input : Secret key s in R_Q^K
Output : Recovered message recovered in B^32
 */
void dec(polyvec *ct_u, poly *ct_v, polyvec *s, uint8_t recovered[32], polyvec *e, polyvec *r) {
    /* ct_v - s^T * ct_u */
    poly mp, stu;
    memset(&mp, 0, sizeof(mp));
    memset(&stu, 0, sizeof(stu));

    polyvec_mul(&stu, s, ct_u); // stu <- s^T * u
    poly_sub(&mp, ct_v, &stu); // mp <- ct_v - s^T * ct_u

    poly decryptionerror;
    memset(&decryptionerror, 0, sizeof(poly)); polyvec_mul(&decryptionerror, e,r);
#if _DEBUG2 == 1
    printpoly("etr (Decryption error e^T * r in RQ) = ", &decryptionerror);
#elif _DEBUG == 1
    printpoly("mp (Polynomial forms of recovered message 'recovered' in R_Q) = ", &mp);
#endif

    /* Rounding */
    poly_tomsg(recovered, &mp); // recovered(32-byte array) <- mp(poly)

}


/* Module-RLWE test function */
void module_rlwe_test() {

    srand(time(NULL));



    /* Genereate key pair
        * public key : A(메모리 절약을 위해 A를 만들 수 있는 publicseed로 대체), b
        * secret key : s
        * */
    polyvec b, s;
    uint8_t publicseed[32];
    memset(&b, 0, sizeof(b));
    memset(&s, 0, sizeof(s));
    memset(&publicseed, 0, sizeof(publicseed));

    polyvec e; // 23.1.27 임시추가
    e = keygen(&b, publicseed, &s);

    /* Setting message */
#if MESSAGETYPE == TESTVECTOR
    // uint8_t m[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09,
    //                 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 
    //                 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d,
    //                 0x1e, 0x1f };
    uint8_t m[32] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00 };
#else
    uint8_t m[32];
    printf(" * Input your message(less than 32-byte) > ");
    gets(m);
    m[strlen(m)] = 0;
#endif 

#if CHECK_ERROR_RANGE == 1
    // uint8_t m[32] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
                    // 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
                    // 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                    // 0xff, 0xff };
    uint8_t m[32] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00 };
#endif
    printstate(" Message m (32-byte array) = ", m);



    /* Encryption
        * Message : m (in B^32)
        * Public seed : publicseed (in B^32)
        * Public key : b (in R_Q^K)
        * Ciphertext : ct_u (in R_Q^K), ct_v (in R_Q)*/
    polyvec ct_u;
    poly ct_v;
    uint8_t coin[32];
    memset(&ct_u, 0, sizeof(ct_u));
    memset(&ct_v, 0, sizeof(ct_v));

    get_randomarray(coin);

    polyvec r; // 23.1.27 임시추가
    r = enc(&ct_u, &ct_v, m, &b, publicseed, coin);



    
    /* Decryption
     * Ciphertext : ct_u (in R_Q^K), ct_v (in R_Q)
     * Secret key : s (in R_Q^K)
     * Decrypted message : recovered (in B^32)
    */
    uint8_t recovered[32];
    memset(&recovered, 0, sizeof(recovered));
    dec(&ct_u, &ct_v, &s, recovered, &e, &r);




 // message 헥사값으로 출력하느 ㄴ옵션 추가


    /* Print our state into stdout 
        (여기에는 함수 내부에서 쓰고 버리는 값을 제외하고 출력합니다.)
        (UI 매크로의 값을 1로 설정하면 함수 내부에서 쓰고 버리는 값 또한 출력합니다.)
    */
    printf(" \n\n\n============================================================================================================\n");
    printf("                          Print the input and output values of a function 'keygen', 'enc', 'dec' \n");
    
    printf(" \n------------------------------------------------ [ KeyGen ] ------------------------------------------------\n");
    printf("                                            ( b, s, publicseed ) \n");
    printf(" \n------------------------------------------------------------------------------------------------------------\n");
    printvector("Public key b (in R_Q^K) = ", &b);
    printvector("Secret key s (in R_Q^K) = ", &s);
    printstate("Public seed publicseed (in B^32) = ", publicseed);

    printf(" \n------------------------------------------------- [ Enc ] --------------------------------------------------\n");
    printf("                                            ( coin, ct_u, ct_v ) \n");
    printf(" \n------------------------------------------------------------------------------------------------------------\n");
    printstate("Random coin (in B^32) = ", coin);
    printvector("Ciphertext ct_u (in R_Q^K) = ", &ct_u);
    printpoly("Ciphertext ct_v (in R_Q) = ", &ct_v);

    printf(" \n------------------------------------------------- [ Dec ] --------------------------------------------------\n");
    printf("                                                 ( recovered ) \n");
    printf(" \n------------------------------------------------------------------------------------------------------------\n");
#if MESSAGETYPE == TESTVECTOR
    printstate(" Recovered message (in B^32) = ", recovered);
#elif MESSAGETYPE == STDIN
    printf(" Message recovered (in B^32) = %s\n", recovered);
#endif
    printf(" ============================================================================================================\n\n");


}

 
 


/* main function */
int main(void) {


    /* Module-RLWE test */
    module_rlwe_test();


    return 0;
}



/* Print 32-byte array */
void printstate(char *msg, uint8_t state[32]) {
    printf("    * %s \n", msg);
    printf("     [ ");

    if ( PRINTTYPE == DECIMAL ) {
        for (int i=0; i<32; i++)
            printf("%3d ", state[i]);
    }
    else if ( PRINTTYPE == HEX ) {
        for (int i=0; i<32; i++)
            printf("%2x ", state[i]);
    }
    printf(" ] \n");
    printf("\n");
}

/* Print polynomial */
void printpoly(char *msg, poly *polynomial) {

    printf("    * %s \n", msg);
    printf("     [ ");
    if ( PRINTTYPE == DECIMAL ) {
        for (int j=0; j<N; j++) {
            printf("%4d ", polynomial->coef[j]);
            if ( j%16==15 )
                printf("\n        ");
        }
    }
    else if ( PRINTTYPE == HEX ) {
        for (int j=0; j<N; j++) {
            printf("%3x ", polynomial->coef[j]);
            if ( j%16==15 )
                printf("\n        ");
        }
    }

    printf(" ] \n");
    printf("\n");
}

/* Print polynomial vector */
void printvector(char *msg, polyvec *vector) {

    printf("    * %s \n", msg);
    if ( PRINTTYPE == DECIMAL ) {
        for (int i=0; i<K; i++) {
            printf("     [ ");
            for (int j=0; j<N; j++) {
                printf("%4d ", vector->vec[i].coef[j]);
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
                printf("%3x ", vector->vec[i].coef[j]);
                if ( j%16==15 )
                    printf("\n        ");
            }
            printf(" ] \n");
        }
    }
    printf("\n");
}

/* Print polynomial vector array(k x k matrix) */
void printmatrix(char *msg, polyvec vector[K]) {

    printf("    * %s \n", msg);



    printf("     [ ");
    if ( PRINTTYPE == DECIMAL ) {
        for (int k=0; k<K; k++ ) {
            for (int i=0; i<K; i++) {
                printf("  [ ");
                for (int j=0; j<N; j++) { 
                    printf("%4d ", vector[k].vec[i].coef[j]);
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
                    printf("%3x ", vector[k].vec[i].coef[j]);
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


// EOF



/*************************************************
* Name:        csubq
*
* Description: Conditionallly subtract q
*
* Arguments:   - int16_t x: input integer
*
* Returns:     a - q if a >= q, else a
**************************************************/
int16_t csubq(int16_t a) {
    int16_t res;
    res = a - Q;
    res += (res >> 15) & Q;
    return res;
}
/*************************************************
* Name:        poly_csubq
*
* Description: Applies conditional subtraction of q to each coefficient
*              of a polynomial. For details of conditional subtraction
*              of q see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
// * Returns:     a - q if a >= q, else a // yb
void poly_csubq(poly *res, poly *a) {
    unsigned int i;
    for(i=0;i<N;i++)
        res->coef[i] = csubq(a->coef[i]);
}


// /*************************************************
// * Name:        poly_tomsg
// *
// * Description: Convert polynomial to 32-byte message
// *
// * Arguments:   - uint8_t *msg: pointer to output message
// *              - poly *a:      pointer to input polynomial
// **************************************************/
// void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *a)
// {
//   unsigned int i,j;
//   uint16_t t;

//   poly_csubq(a); // == mod Q

//   for(i=0;i<N/8;i++) {
//     msg[i] = 0;
//     for(j=0;j<8;j++) {
//       t = ((((uint16_t)a->coef[8*i+j] << 1) + Q/2)/Q) & 1;
//       msg[i] |= t << j;
//     }
//   }
// }