//  Matrix.cpp

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <random>

#include "matrix.hpp"

Matrix::Matrix() {
    m_ = 0;
    n_ = 0;

    mb_ = 0;
    nb_ = 0;

    p_ = 0;
    q_ = 0;

    top_ = nullptr;
}

Matrix::Matrix(const unsigned int &m, const unsigned int &n) {
    
    m_ = m;
    n_ = n;
    
    mb_ = m;
    nb_ = n;
    
    p_ = 1;
    q_ = 1;
    
    top_ = new double[m_ * n_];
};

Matrix::Matrix(const unsigned int &m, const unsigned int &n,
               const unsigned int &ts) {
    
    m_ = m;
    n_ = n;
    
    mb_ = ts;
    nb_ = ts;
    
    p_ = (m_ % mb_ == 0) ? m_ / mb_ : m_ / mb_ + 1;
    q_ = (n_ % nb_ == 0) ? n_ / nb_ : n_ / nb_ + 1;
    
    top_ = new double[m_ * n_];
};

Matrix::Matrix(const unsigned int &m, const unsigned int &n,
               const unsigned int &mb, const unsigned int &nb) {
    
    m_ = m;
    n_ = n;
    
    mb_ = mb;
    nb_ = nb;
    
    p_ = (m_ % mb_ == 0) ? m_ / mb_ : m_ / mb_ + 1;
    q_ = (n_ % nb_ == 0) ? n_ / nb_ : n_ / nb_ + 1;
    
    top_ = new double[m_ * n_];
};

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
    
    for (int i = 0; i < m_ * n_; i++)
        top_[i] = dist(engine);
}

double* Matrix::elm(const unsigned int &ti, const unsigned int &tj) const {
    assert(ti >= 0);
    assert(ti < p_);
    assert(tj >= 0);
    assert(tj < q_);
    
    unsigned int pos = 0;
    
    pos += ti * (mb_ * n_);
    pos += (ti == p_ - 1) ? tj * (m_ % mb_) * nb_ : tj * mb_ * nb_;
    
    return &top_[pos];
}

double* Matrix::elm(const unsigned int &ti, const unsigned int &tj,
                    const unsigned int &i, const unsigned int &j) const {
    assert(i >= 0);
    assert(i < (ti == (p_ - 1) ? m_ % mb_ : mb_));
    assert(j >= 0);
    assert(j < (tj == (q_ - 1) ? n_ % nb_ : nb_));
    
    assert(ti >= 0);
    assert(ti < p_);
    assert(tj >= 0);
    assert(tj < q_);
    
    unsigned int pos = 0;
    
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
    
    for (unsigned int i = 0; i < m_; i++) {
        for (unsigned int j = 0; j < n_; j++) {
            unsigned int p = 0;

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
    
    for (unsigned int i = 0; i < m_ * n_; i++)
        top_[i] = T.top_[i];
    
    return *this;
}

Matrix Matrix::operator+(const Matrix &T) const {
    assert(m_ == T.m_);
    assert(n_ == T.n_);

    Matrix M(m_, n_, mb_, nb_);
    for(int i=0; i<m_*n_; ++i){
        M[i] = (*this)[i] + T[i];
    }

    return M;
}

Matrix Matrix::operator-(const Matrix &T) const {
    assert(m_ == T.m_);
    assert(n_ == T.n_);

    Matrix M(m_, n_, mb_, nb_);
    for(int i=0; i<m_*n_; ++i){
        M[i] = (*this)[i] - T[i];
    }

    return M;
}

bool Matrix::operator==(const Matrix &T) {
    if (m_ != T.m_ || n_ != T.n_) {
        return false;
    }
    for (int i = 0; i < m_ * n_; i++) {
        if (top_[i] != T.top_[i]){
            return false;
        }
    }
    return true;
}

double &Matrix::operator[](const unsigned int &i) const {
    assert(i >= 0);
    assert(i < m_ * n_);
    
    return top_[i];
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */

double &Matrix::operator()(const unsigned int &i, const unsigned int &j)
const {
    assert(i >= 0);
    assert(i < m_);
    assert(j >= 0);
    assert(j < n_);


    const unsigned int ti = i / mb_;
    const unsigned int tj = j / nb_;
    const unsigned int tp = i % mb_;
    const unsigned int tq = j % nb_;

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
    
    unsigned int m_ = ma.m();
    unsigned int n_ = ma.n();
    unsigned int mb_ = ma.mb();
    unsigned int nb_ = ma.nb();
    unsigned int p_ = ma.p();
    unsigned int q_ = ma.q();
    
    for (unsigned int i = 0; i < m_; i++) {
        for (unsigned int j = 0; j < n_; j++) {
            unsigned int p = 0;

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
