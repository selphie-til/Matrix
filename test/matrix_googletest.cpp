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

TEST(test_matrix_complex_float, equal) {
    Matrix<std::complex<float>> A(1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<std::complex<float>> B(1000, 1000, 30);
    B = A;
    EXPECT_TRUE( A==B );
}

TEST(test_matrix_complex_float, plus) {
    Matrix<std::complex<float>> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<std::complex<float>> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<std::complex<float>> C(1000,1000,30);
    C = A + B;
    Matrix<std::complex<float>> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] + B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_complex_float, minus) {
    Matrix<std::complex<float>> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<std::complex<float>> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<std::complex<float>> C(1000,1000,30);
    C = A - B;
    Matrix<std::complex<float>> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] - B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_complex_float, elm) {
    Matrix<std::complex<float>> A( 1000, 1000, 30);
    A.gen_rnd_elm();

    EXPECT_TRUE( *(A.elm( 20, 20, 5, 5)) == A(605, 605));
}

TEST(test_matrix_complex_float, zero){
    Matrix<std::complex<float>> A(1000, 1000);
    std::complex<float> c{0.0,0.0};
    A.zero();
    for(uint32_t i=0; i<1000; ++i){
        for(uint32_t j=0; j<1000; ++j){
            EXPECT_TRUE( A(i,j) == c );
        }
    }
}

TEST(test_matrix_complex_double, equal) {
    Matrix<std::complex<double>> A(1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<std::complex<double>> B(1000, 1000, 30);
    B = A;
    EXPECT_TRUE( A==B );
}

TEST(test_matrix_complex_double, plus) {
    Matrix<std::complex<double>> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<std::complex<double>> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<std::complex<double>> C(1000,1000,30);
    C = A + B;
    Matrix<std::complex<double>> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] + B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_complex_double, minus) {
    Matrix<std::complex<double>> A( 1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix<std::complex<double>> B( 1000, 1000, 30);
    B.gen_rnd_elm();

    Matrix<std::complex<double>> C(1000,1000,30);
    C = A - B;
    Matrix<std::complex<double>> D(1000, 1000, 30);
    for(uint32_t i=0; i< D.m()*D.n(); ++i){
        D[i] = A[i] - B[i];
    }

    EXPECT_TRUE( C==D );
}

TEST(test_matrix_complex_double, elm) {
    Matrix<std::complex<double>> A( 1000, 1000, 30);
    A.gen_rnd_elm();

    EXPECT_TRUE( *(A.elm( 20, 20, 5, 5)) == A(605, 605));
}

TEST(test_matrix_complex_double, zero){
    Matrix<std::complex<double>> A(1000, 1000);
    std::complex<double> z {0.0, 0.0};
    A.zero();
    for(uint32_t i=0; i<1000; ++i){
        for(uint32_t j=0; j<1000; ++j){
            EXPECT_TRUE( A(i,j) == z );
        }
    }
}