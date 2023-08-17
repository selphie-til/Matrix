//  Matrix.cpp

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <random>

#include "matrix.hpp"

Matrix::~Matrix() {
    
    delete[] top_;
    top_ = nullptr;
}

void Matrix::gen_rnd_elm() {
    
    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<> dist(0.0, 1.0);
    
    for (uint64_t i = 0; i < m_*n_; ++i) {
        top_[i] = dist(engine);
    }
}

double* Matrix::elm(const uint32_t &ti, const uint32_t &tj) const {
    assert(ti >= 0);
    assert(ti < p_);
    assert(tj >= 0);
    assert(tj < q_);
    
    uint32_t pos = 0;
    
    pos += ti * (mb_ * n_);
    pos += (ti == p_ - 1) ? tj * (m_ % mb_) * nb_ : tj * mb_ * nb_;
    
    return &top_[pos];
}

double* Matrix::elm(const uint32_t &ti, const uint32_t &tj,
                    const uint32_t &i, const uint32_t &j) const {
    assert(i >= 0);
    assert(i < (ti == (p_ - 1) ? m_ % mb_ : mb_));
    assert(j >= 0);
    assert(j < (tj == (q_ - 1) ? n_ % nb_ : nb_));
    
    assert(ti >= 0);
    assert(ti < p_);
    assert(tj >= 0);
    assert(tj < q_);
    
    uint64_t pos = 0;
    
    pos += ti * (mb_ * n_);
    pos += (ti == p_ - 1) ? tj * (m_ % mb_) * nb_ : tj * mb_ * nb_;
    pos += i;
    pos += (ti == p_ - 1) ? j * (m_ % mb_) : j * mb_;
    
    return &top_[pos];
}

void Matrix::file_out(const char *fname) {
    
    std::ofstream matf(fname);
    if (!matf) {
        std::cerr << "Unable to open " << fname << std::endl;
        exit(1);
    }
    
    matf.precision(5);
    
    for (uint32_t i = 0; i < m_; i++) {
        for (uint32_t j = 0; j < n_; j++) {
            uint64_t p = 0;

            // (i / mb_) ti
            if ((i / mb_) != p_ - 1) {
                p += (i / mb_) * mb_ * n_;  // tiの場所
                p += (j / nb_) * mb_ * nb_; // tjの場所

                p += (j % nb_) * mb_ + (i % mb_); // i,jの場所
            } else {
                // ここで最終行の時と差別化
                p += (p_ - 1) * mb_ * n_; // tiの場所
                p += (m_ % mb_ == 0) ? (j / nb_) * mb_ * nb_
                                     : (j / nb_) * (m_ % mb_) * nb_; // tjの場所
                p += (m_ % mb_ == 0) ? (j % nb_) * mb_ + (i % mb_)
                                     : (j % nb_) * (m_ % mb_) + (i % mb_);
            }

            matf << top_[p] << " ";
        }
        matf << std::endl;
    }
    matf.close();
}

Matrix& Matrix::operator=(const Matrix &T){
    assert(m_ == T.m_);
    assert(n_ == T.n_);

    m_ = T.m_;
    n_ = T.n_;
    mb_ = T.mb_;
    nb_ = T.nb_;
    p_ = T.p_;
    q_ = T.q_;
    
    for (uint64_t i = 0; i < m_*n_; i++)
        top_[i] = T.top_[i];
    
    return *this;
}

Matrix Matrix::operator+(const Matrix &T) const {
    assert(m_ == T.m_);
    assert(n_ == T.n_);

    Matrix M(m_, n_, mb_, nb_);
    for (uint64_t i=0; i<m_*n_; ++i){
        M[i] = (*this)[i] + T[i];
    }

    return M;
}

Matrix Matrix::operator-(const Matrix &T) const {
    assert(m_ == T.m_);
    assert(n_ == T.n_);

    Matrix M(m_, n_, mb_, nb_);
    for (int i=0; i<m_*n_; ++i){
        M[i] = (*this)[i] - T[i];
    }

    return M;
}

bool Matrix::operator==(const Matrix &T) {
    if(m_ != T.m_ || n_ != T.n_) {
        return false;
    }
    for (uint64_t i = 0; i < m_*n_; i++) {
        if (top_[i] != T.top_[i]){
            return false;
        }
    }
    return true;
}

double &Matrix::operator[](const uint64_t &i) const {
    assert(i >= 0);
    assert(i < m_*n_);
    
    return top_[i];
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */

double &Matrix::operator()(const uint32_t &i, const uint32_t &j)
const {
    assert(i >= 0);
    assert(i < m_);
    assert(j >= 0);
    assert(j < n_);


    const uint32_t ti = i / mb_;
    const uint32_t tj = j / nb_;
    const uint32_t tp = i % mb_;
    const uint32_t tq = j % nb_;

    //間違っている ti, tj の位置を考慮しないといけない
    return top_[ (mb_*n_)*ti + tp + (mb_*nb_)*tj + (mb_*tq)];
}
/*
  template<typename CharT, typename Traits>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os, const Matrix& ma){
    
    int m_ = ma.m();
    int n_ = ma.n();
    int mb_ = ma.mb();
    int nb_ = ma.nb();
    int p_ = ma.p();
    int q_ = ma.q();
    
    for (int i=0; i<m_; i++) {
    for (int j=0; j<n_; j++) {
    int p=0;
    
    // (i / mb_) ti
    if( (i / mb_) != p_-1){
    p += (i / mb_)*mb_*n_; //tiの場所
    p += (j / nb_)*mb_*nb_;//tjの場所
    
    
    p += (j % nb_)*mb_ + (i % mb_);//i,jの場所
    }
    else{
    //ここで最終行の時と差別化
    p += (p_-1)*mb_*n_;//tiの場所
    p += (m_%mb_ == 0) ? (j/nb_)*mb_*nb_ :
    (j/nb_)*(m_%mb_)*nb_;//tjの場所 p += (m_%mb_ == 0) ? (j%nb_)*mb_ + (i%mb_) :
    (j%nb_)*(m_%mb_) + (i%mb_);
    }
    
    os << ma[p] << " ";
    }
    os << std::endl;
    }
    
    return os;
    }
*/

std::ostream &operator<<(std::ostream &os, const Matrix &ma) {
    
    uint32_t m_ = ma.m();
    uint32_t n_ = ma.n();
    uint32_t mb_ = ma.mb();
    uint32_t nb_ = ma.nb();
    uint32_t p_ = ma.p();
    uint32_t q_ = ma.q();
    
    for (uint32_t i = 0; i < m_; i++) {
        for (uint32_t j = 0; j < n_; j++) {
            uint32_t p = 0;

            // (i / mb_) ti
            if ((i / mb_) != p_ - 1) {
                p += (i / mb_) * mb_ * n_;  // tiの場所
                p += (j / nb_) * mb_ * nb_; // tjの場所

                p += (j % nb_) * mb_ + (i % mb_); // i,jの場所
            } else {
                // ここで最終行の時と差別化
                p += (p_ - 1) * mb_ * n_; // tiの場所
                p += (m_ % mb_ == 0) ? (j / nb_) * mb_ * nb_
                                     : (j / nb_) * (m_ % mb_) * nb_; // tjの場所
                p += (m_ % mb_ == 0) ? (j % nb_) * mb_ + (i % mb_)
                                     : (j % nb_) * (m_ % mb_) + (i % mb_);
            }

            os << ma[p] << " ";
        }
        os << std::endl;
    }

    return os;
}

void Matrix::zero() {
    for (uint64_t i = 0; i < m_*n_; ++i)
        top_[i] = 0.0;
}
