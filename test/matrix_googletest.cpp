#include <complex>
#include "matrix.hpp"
#include "gtest/gtest.h"

TEST(test_matrix_float, equal) {
    Matrix<float> A(1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<float> B(1000, 1000, 30);
    B = A;
    EXPECT_TRUE( A==B );
}

TEST(test_matrix_float, plus) {
    Matrix<float> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<float> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<float> C(1000,1000,30);
    C = A + B;
    Matrix<float> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] + B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_float, minus) {
    Matrix<float> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<float> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<float> C(1000,1000,30);
    C = A - B;
    Matrix<float> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] - B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_float, elm) {
    Matrix<float> A( 1000, 1000, 30);
    A.gen_rnd_elm();

    EXPECT_TRUE( *(A.elm( 20, 20, 5, 5)) == A(605, 605));
}

TEST(test_matrix_float, zero){
    Matrix<float> A(1000, 1000);
    A.zero();
    for(uint32_t i=0; i<1000; ++i){
        for(uint32_t j=0; j<1000; ++j){
            EXPECT_TRUE( A(i,j) == 0 );
        }
    }
}

TEST(test_matrix_double, equal) {
    Matrix<double> A(1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<double> B(1000, 1000, 30);
    B = A;
    EXPECT_TRUE( A==B );
}

TEST(test_matrix_double, plus) {
    Matrix<double> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<double> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<double> C(1000,1000,30);
    C = A + B;
    Matrix<double> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] + B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_double, minus) {
    Matrix<double> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<double> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<double> C(1000,1000,30);
    C = A - B;
    Matrix<double> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] - B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_double, elm) {
    Matrix<double> A( 1000, 1000, 30);
    A.gen_rnd_elm();

    EXPECT_TRUE( *(A.elm( 20, 20, 5, 5)) == A(605, 605));
}

TEST(test_matrix_double, zero){
    Matrix<double> A(1000, 1000);
    A.zero();
    for(uint32_t i=0; i<1000; ++i){
        for(uint32_t j=0; j<1000; ++j){
            EXPECT_TRUE( A(i,j) == 0 );
        }
    }
}
