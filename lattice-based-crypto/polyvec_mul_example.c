
/* 3. polyvec_mul example (Yongbhin Kim) */



// 3-1. multiplication of polyvec(1 row) and polyvec(1 col)
// 3-2. multiplication of matrix(3 row) and polyvec(1 row)
// row, col consists 3 vector (i.e. 3 entry of ring R_3329)

#include "module_lwe_example.h"

typedef struct {
    poly vec[3];     // 한 행의 원소의 개수는 k 일 것임 (2/3/4 at Kyber512/768/1024) // poly들을 모은 것
} polyvec;           // 원소가 R_3329로 이루어진 행렬의 한 행 (entry of A or b which is used in Crystals-Kyber as pk), primitive polynomial is <x^256+1>

typedef struct {
    polyvec pvec[3]; // polyvec 3개를 모아놓음. 즉 3 행을 3개 모아놓은 행렬
} matrix;

void polyadd(poly *R, poly *a, const poly *b) {

    uint16_t c = KYBER_N;
    uint16_t q = KYBER_Q;
    uint32_t result = 0;

    for (int i = 0; i < c; i++) {
        result = (a->coef[i] + b->coef[i]);
        R->coef[i] = result % KYBER_Q;
    }
}

void mul_row_with_vec(poly *res, const polyvec *row, const polyvec *col) {
    poly temppoly;

    // 이 작업을 세번 반복하면 행렬곱 완성
    for (int i=0; i<3; i++) {
        schoolbook_mul(&temppoly, &row->vec[i], &col->vec[i]);
        polyadd(res, res, &temppoly);
    }
}

void mul_row_with_vec_test() {
    polyvec A, s, res;
    poly res_entry;
    memset(&A, 0, sizeof(polyvec));
    memset(&s, 0, sizeof(polyvec));
    memset(&res, 0, sizeof(polyvec));
    memset(&res_entry, 0, sizeof(poly));

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

    mul_row_with_vec(&res_entry, &A, &s);

    printf(" [ Do polymul ] : c = a * b \n\n");
    for (int i = 0; i < KYBER_N; i++)
        printf("   As[%d] = %d\n", i, res_entry.coef[i]);
}

void mul_mat_with_vec(polyvec *res, polyvec (*A_row)[3], polyvec *s) { //poly *res, const polyvec *row, const polyvec *col) {
    // mul one row of A 3 times
    mul_row_with_vec(&res->vec[0], A_row[0], s); // res의 첫번째 poly(entry) <- A의 첫번째 row와 s의 곱
    mul_row_with_vec(&res->vec[1], A_row[0], s); //
    mul_row_with_vec(&res->vec[2], A_row[0], s); //

}

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
    A_row[2].vec[2].coef[0] = 1; // 

    s.vec[0].coef[0] = 0; // s = [ 0 ]
    s.vec[1].coef[0] = 1; //     [ 1 ]
    s.vec[2].coef[1] = 1; //     [ 1 ] == col(0, 1, x)


    // mul one row of A 3 times
    mul_mat_with_vec(&res, &A_row, &s);

    // print res
    printf(" [ Test ] : mat x vec (res = A * s) over R_3329, primitive polynimial is x^256 + 1 \n\n");
    uint32_t decimalvalue;
    for (int i=0; i<3; i++) {
        decimalvalue = 0;
        printf("[ ");
        for (int j=0; j<4; j++) { // coef index
            printf(" %3d*x^%d ", res.vec[i].coef[j], j);
            decimalvalue += res.vec[i].coef[j]<<j; // 출력을 위한 value, 실제로는 빅넘 범위임
        }
        printf(" ( == %d) ] over R_q \n", decimalvalue);
    }
}

// main function
int main(void) {

    /* poly row * col test */
    // mul_row_with_vec_test(); // polyvec('row of A') * polyvec(col; vector 's')

    /* matrix 'A' * polyvector s */
    mul_mat_with_vec_test();


    return 0;
}

// EOF