#include "matrix.hpp"
#include "gtest/gtest.h"

TEST(test_matrix, equal) {
    Matrix A(1000, 1000, 30);
    A.gen_rnd_elm();
    Matrix B(1000, 1000, 30);
    B = A;
    EXPECT_TRUE( A==B );
}
