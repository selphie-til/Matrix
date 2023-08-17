#include "matrix.hpp"
#include "gtest/gtest.h"

TEST(test_matrix, equal) {
    Matrix A(1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix B(1000, 1000, 30);
    B = A;
    EXPECT_TRUE( A==B );
}

TEST(test_matrix, plus) {
    Matrix A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix C(1000,1000,30);
    C = A + B;
    Matrix D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] + B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix, minus) {
    Matrix A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix C(1000,1000,30);
    C = A - B;
    Matrix D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] - B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix, elm) {
    Matrix A( 1000, 1000, 30);
    A.gen_rnd_elm();

    EXPECT_TRUE( *(A.elm( 20, 20, 5, 5)) == A(605, 605));
}

TEST(test_matrix, zero){
    Matrix A(1000, 1000);
    A.zero();
    for(uint32_t i=0; i<1000; ++i){
        for(uint32_t j=0; j<1000; ++j){
            EXPECT_TRUE( A(i,j) == 0 );
        }
    }
}