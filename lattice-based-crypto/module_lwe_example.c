
/* 4. Module-RLWE example (Yongbhin Kim) */

// [Key setup]
// pk : (A, b) := As + e
//      A is some element of R_3329
//      b is element of R_3329 which is computated by A*s + e
// sk : s 
//      s is some element of R_3
// secret value : s, e 
//      e is some element of R_3 and it means error vector

// [Encryotion]
// ct_u <- ∑A
// ct_v <- ∑b + ceiling(q/2)*m

// [Decryption]
// result <- ct_v - s*ct_u

// [Rounding]
// recovered <- if result > q/2 then 1 else then 0

#include "module_lwe_example.h"


typedef struct {
    poly vec[3];     // 한 행의 원소의 개수는 k 일 것임 (2/3/4 at Kyber512/768/1024) // poly들을 모은 것
} polyvec;           // 원소가 R_3329로 이루어진 행렬의 한 행 (element of A or b which is used in Crystals-Kyber as pk), primitive polynomial is <x^256+1>


void polyadd(poly *res, poly *a, const poly *b);
void polysub(poly *res, poly *a, const poly *b);
void polyvecadd(polyvec *res, polyvec *a, const polyvec *b);
void mul_row_with_vec(poly *res, const polyvec *row, const polyvec *col);
void mul_row_with_vec_test();
void mul_mat_with_vec(polyvec *res, polyvec A_row[3], polyvec *s);
void mul_mat_with_vec_test();
void mul_vec_with_mat(polyvec *res, polyvec *s, polyvec A_row[3]);
void module_lwe_test();

void printpoly(char *msg, poly *polynomial);
void printmatrix(char *msg, polyvec vector[3]);
void printvector(char *msg, polyvec *vector);


/* Result(Polynomial) = Polynomial + Polynomial */
void polyadd(poly *res, poly *a, const poly *b) {

    uint16_t c = KYBER_N;
    uint16_t q = KYBER_Q;
    uint32_t result = 0;

    for (int i = 0; i < c; i++) {
        result = (a->coef[i] + b->coef[i]);
        res->coef[i] = result % KYBER_Q;
    }
}

/* Result(Polynomial) = Polynomial - Polynomial */
void polysub(poly *res, poly *a, const poly *b) {

    uint16_t c = KYBER_N;
    uint16_t q = KYBER_Q;
    uint32_t result = 0;

    for (int i = 0; i < c; i++) {
        result = (a->coef[i] + KYBER_Q - b->coef[i]);
        res->coef[i] = result % KYBER_Q;
    }
}

/* Result(Vector) = Vector + Vector */
void polyvecadd(polyvec *res, polyvec *a, const polyvec *b) {
    polyadd(&res->vec[0], &a->vec[0], &b->vec[0]);
    polyadd(&res->vec[1], &a->vec[1], &b->vec[1]);
    polyadd(&res->vec[2], &a->vec[2], &b->vec[2]);
}

/* Result(Polynomial) = RowVector x ColumnVector */
void mul_row_with_vec(poly *res, const polyvec *row, const polyvec *col) {
    poly temppoly;

    // 이 작업을 세번 반복하면 행렬곱 완성
    for (int i=0; i<3; i++) {
        schoolbook_mul(&temppoly, &row->vec[i], &col->vec[i]);
        polyadd(res, res, &temppoly);
    }
}

/* RowVector x ColumnVector test function */
void mul_row_with_vec_test() {
    polyvec A, s, res;
    poly res_element;
    memset(&A, 0, sizeof(polyvec));
    memset(&s, 0, sizeof(polyvec));
    memset(&res, 0, sizeof(polyvec));
    memset(&res_element, 0, sizeof(poly));

    for (int i = 1; i < KYBER_N; i++) {
        A.vec[0].coef[i] = 0;       // test : [ 1, 1, 1 ] x [ 0 ]  == [ x + 1 ]
        A.vec[1].coef[i] = 0;       //                      [ 1 ]
        A.vec[2].coef[i] = 0;       //                      [ x ]
        s.vec[0].coef[i] = 0;
        s.vec[1].coef[i] = 0;
        s.vec[2].coef[i] = 0;
    }
    A.vec[0].coef[0] = 1;           // A = [ 1 0 0 ] 
    A.vec[1].coef[0] = 1;           // A = [ 1 1 0 ]
    A.vec[2].coef[0] = 1;           // A = [ 1 1 1 ]
    s.vec[0].coef[0] = 0;           // s = [ 0 ]
    s.vec[1].coef[0] = 1;           //     [ 1 ]
    s.vec[2].coef[1] = 1;           //     [ 1 ] == col(0, 1, x)

    mul_row_with_vec(&res_element, &A, &s);

    printf(" [ Do polymul ] : c = a * b \n\n");
    for (int i = 0; i < KYBER_N; i++)
        printf("   As[%d] = %d\n", i, res_element.coef[i]);
}

/* Result(Vector) = Matrix x ColumnVector */
void mul_mat_with_vec(polyvec *res, polyvec mat_row[3], polyvec *cvec) { 
    // mul one row of A 3 times
    mul_row_with_vec(&res->vec[0], &mat_row[0], cvec); // res의 첫번째 poly(element) <- A의 첫번째 row와 s의 곱
    mul_row_with_vec(&res->vec[1], &mat_row[1], cvec); //
    mul_row_with_vec(&res->vec[2], &mat_row[2], cvec); //

}

/* Result(Vector) = RowVector x Matrix */
void mul_vec_with_mat(polyvec *res, polyvec *rvec, polyvec mat_row[3]) {
    // mul one col of A 3 times
    polyvec mat_col[3];
    memset(&mat_col, 0, sizeof(mat_col));

    // Realloc mat_col <- mat_row->vec[0].coef[0] to mat_row->vec[3].coef[255]
    for (int i=0; i<3; i++ ) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<KYBER_N; k++) {
                mat_col[i].vec[j].coef[k] = mat_row[j].vec[i].coef[k];
            }
        }
    }

    mul_row_with_vec(&res->vec[0], &mat_col[0], rvec); // 
    mul_row_with_vec(&res->vec[1], &mat_col[1], rvec); //
    mul_row_with_vec(&res->vec[2], &mat_col[2], rvec); //

}

/* Result(Vector) = Matrix x ColumnVector */
void mul_mat_with_vec_test() {
    // polyvec A_row[3] == matrix A;

    polyvec A_row[3], res, s; // Let A := A_row[0], A_row[1], A_row[3], result is 1x3 vector
    memset(&A_row, 0, sizeof(A_row));
    memset(&res, 0, sizeof(polyvec));
    memset(&s, 0, sizeof(polyvec));

    // [test]    
    //           res             A          s 
    //        [ x + 1 ]     [ 1, 1, 1 ]   [ 0 ]
    //        [ x + 1 ]  <- [ 1, 1, 1 ] x [ 1 ]
    //        [ x + 1 ]     [ 1, 1, 1 ]   [ x ]

    // A = [ 1 1 1 ]
    //     [ 0 0 0 ]
    //     [ 0 0 0 ]
    A_row[0].vec[0].coef[0] = 1; // 
    A_row[0].vec[1].coef[0] = 1; // 
    A_row[0].vec[2].coef[0] = 1; // 

    // A = [ 1 1 1 ]
    //     [ 1 1 1 ]
    //     [ 0 0 0 ]
    A_row[1].vec[0].coef[0] = 1; //     
    A_row[1].vec[1].coef[0] = 1; //     
    A_row[1].vec[2].coef[0] = 1; // 

    // A = [ 1 1 1 ]
    //     [ 1 1 1 ]
    //     [ 1 1 1 ]
    A_row[2].vec[0].coef[0] = 1; //     
    A_row[2].vec[1].coef[0] = 1; //     
    A_row[2].vec[2].coef[0] = 1; // 2행2열의 값 // yb

    s.vec[0].coef[0] = 0; // s = [ 0 ]
    s.vec[1].coef[0] = 1; //     [ 1 ]
    s.vec[2].coef[1] = 1; //     [ 1 ] == col(0, 1, x)
#if 0
    printf("Matrix A\n");
    printf(" ====================================\n");
    uint32_t decimalvalue;
    for (int k=0; k<3; k++ ) {
        printf("  [ ");
        for (int i=0; i<3; i++) {
            decimalvalue = 0;
            for (int j=0; j<4; j++) { // coef index
                printf("%3d*x^%d ", A_row[k].vec[i].coef[j], j);
                // decimalvalue += pk_b.vec[i].coef[j]<<j; 
            }
            printf(", ");
        }
        printf(" ] ");
        printf("over R_q \n");
    }
#endif
    // mul one row of A 3 times
    mul_mat_with_vec(&res, A_row, &s);

    // print result
    printf(" [ Test ] : mat x vec (res = A * s) over R_3329, primitive polynimial is x^256 + 1 \n\n");
    for (int i=0; i<3; i++) {
        printf("[ ");
        for (int j=0; j<4; j++) { // coef index
            printf(" %3d*x^%d ", res.vec[i].coef[j], j);
        }
        printf(" ] over R_q\n");
    }
}

/* Module-LWE test function */
void module_lwe_test() {

    /* [ Generate key pair ] *
     *  (pk_A, pk_b), sk_s   */
    polyvec pk_A[3], pk_b, sk_s, secret_e, As;
    memset(&pk_A, 0, sizeof(pk_A));
    memset(&pk_b, 0, sizeof(polyvec));
    memset(&sk_s, 0, sizeof(polyvec));
    memset(&secret_e, 0, sizeof(polyvec));
    memset(&As, 0, sizeof(polyvec));

#if 1
    //          pk_b           pk_A       sk_s    secret_e
    //        [ x + 2 ]     [ 1, 1, 1 ]   [ 0 ]    [ 1 ]
    //        [ x + 2 ]  <- [ 1, 1, 1 ] x [ 1 ] +  [ 1 ]
    //        [ x + 2 ]     [ 1, 1, 1 ]   [ x ]    [ 1 ]
    pk_A[0].vec[0].coef[0] = 1;  // 
    pk_A[0].vec[1].coef[0] = 1;  // 
    pk_A[0].vec[2].coef[0] = 1;  // 

    pk_A[1].vec[0].coef[0] = 1;  //     
    pk_A[1].vec[1].coef[0] = 1;  //     
    pk_A[1].vec[2].coef[0] = 1;  // 

    pk_A[2].vec[0].coef[0] = 1;  //     
    pk_A[2].vec[1].coef[0] = 1;  //     
    pk_A[2].vec[2].coef[0] = 1;  // 

    sk_s.vec[0].coef[0] = 0;     // s = [ 0 ]
    sk_s.vec[1].coef[0] = 1;     //     [ 1 ]
    sk_s.vec[2].coef[1] = 1;     //     [ x ] == col(0, 1, x)

    secret_e.vec[0].coef[0] = 1; // e = [ 1 ] 
    secret_e.vec[1].coef[0] = 1; //     [ 1 ] 
    secret_e.vec[2].coef[0] = 1; //     [ 1 ] == col(1, 1, 1)

    mul_mat_with_vec(&As, pk_A, &sk_s); // As <- A x s

    polyvecadd(&pk_b, &As, &secret_e);  // b  <- As + e
    
    /* [Setting message]  */
    srand(time(NULL));
    uint16_t message = rand()&0x01;
    uint16_t temp = ((KYBER_Q+1)/2) * message;
    poly q2m;
    memset(&q2m, 0, sizeof(poly));
    q2m.coef[0]=temp; // ceiling(q/2)*m

    /* [ Encryotion ] *
     *  (ct_u, ct_v) <- ( ∑A, ∑b + ceiling(q/2)*m )  */
    polyvec ct_u, r;
    poly ct_v;
    memset(&ct_u, 0, sizeof(polyvec));
    memset(&r, 0, sizeof(polyvec));
    memset(&ct_v, 0, sizeof(ct_v));

    // ct_u := ∑A := rA, setting vector 'r' belows
    r.vec[0].coef[0] = 1; //     [ 1 ]
    r.vec[1].coef[0] = 1; // r = [ 1 ]
    r.vec[2].coef[0] = 1; //     [ 1 ]

    // Compute ct_u
    mul_vec_with_mat(&ct_u, &r, pk_A);

    // Compute ct_v := ∑b + ceiling(q/2)*m
    mul_row_with_vec(&ct_v, &r, &pk_b); // ∑b
    polyadd(&ct_v, &ct_v, &q2m);

    /* [ Decryption & Rounding ] *
     * res <- ct_v - ct_u x s 
     * if res > ceiling(q/2)*m 
     *  then recovered <- 1
     * else 
     *  then recovered <- 0 */
    
    poly result;
    uint16_t recovered;
    memset(&result, 0, sizeof(poly));
    mul_row_with_vec(&result, &ct_u, &sk_s); 
    polysub(&result, &ct_v, &result);
    if ( result.coef[0] > KYBER_Q/2 )
        recovered = 1;
    else 
        recovered = 0;





    /* Print message */
    printf(" =======================================================================================================================\n");
    printf("   * [ Module-LWE Test]\n");
    printf("     ( Operate Matrix, Vector, Polynomial over R_3329, primitive polynimial is x^256 + 1 ) \n\n");
    printf(" -----------------------------------------------------------------------------------------------------------------------\n");
    printf("   * [ message = %d ]\n", message);

    /* Print keygen */
    printf(" -----------------------------------------------------------------------------------------------------------------------\n");
    printf("\n  [ Keygen ] \n");
    printf("   [public key] \n");
    printmatrix("A = ", pk_A);
    printvector("b = ", &pk_b);

    printf("\n   [private key] \n");
    printvector("s = ", &sk_s);

    printf("\n   [secret value] \n");
    printvector("e = ", &secret_e);

    /* Print ciphertext */
    printf(" -----------------------------------------------------------------------------------------------------------------------\n");
    printf("\n  [ Encrypt ] ");
    printf("\n   [Ciphertext] \n");
    printvector("ct_u = ", &ct_u);
    printpoly("ct_v = ", &ct_v);

    /* Print decrypted & recovered */
    printf(" -----------------------------------------------------------------------------------------------------------------------\n");
    printf("\n  [ Decrypt & Rounding ] ");
    printf("\n   [ Decrypted result ] \n");
    printpoly("decrypted = ", &result);
    printf("\n   [ Recovered ] \n");
    printf("   * recovered = %d\n", recovered);
    printf(" -----------------------------------------------------------------------------------------------------------------------\n");
    if ( message == recovered )
        printf("   * Success to recover!\n");
    printf(" =======================================================================================================================\n");


#endif



}

/* main function */
int main(void) {

    /* Poly row * col test */
    // mul_row_with_vec_test(); // polyvec('row of A') * polyvec(col; vector 's')

    /* Matrix 'A' * polyvector 's' */
    // mul_mat_with_vec_test();

    /* Module-LWE test */
    module_lwe_test();


    return 0;
}

void printvector(char *msg, polyvec *vector) {

    printf("    * %s \n", msg);

    for (int i=0; i<3; i++) {
        printf("     [ ");
        for (int j=0; j<4; j++) {
            printf("%3d*x^%d ", vector->vec[i].coef[j], j);
        }
        printf(" ] over R_q \n");
    }
}

void printmatrix(char *msg, polyvec vector[3]) {

    printf("    * %s \n", msg);

    for (int k=0; k<3; k++ ) {
        printf("     [ ");
        for (int i=0; i<3; i++) {
            for (int j=0; j<4; j++) { 
                printf("%3d*x^%d ", vector[k].vec[i].coef[j], j);
            }
            printf(", ");
        }
        printf(" ] ");
        printf("over R_q \n");
    }
}

void printpoly(char *msg, poly *polynomial) {

    printf("    * %s \n", msg);

    printf("     [ ");
    for (int j=0; j<4; j++)
        printf("%3d*x^%d ", polynomial->coef[j], j);
    printf(" ] over R_q \n");
}

// EOF