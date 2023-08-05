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
    for(int i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] + B[i];
    }

    EXPECT_TRUE( C==D );
}