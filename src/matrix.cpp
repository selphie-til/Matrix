//  Matrix.cpp

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <random>
#include <cstring>
#include <complex>

#include "matrix.hpp"

template<typename T>
Matrix<T>::~Matrix() {
    
    delete[] top_;
    top_ = nullptr;
}

template<>
void Matrix<float>::gen_rnd_elm() {

    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    for (uint64_t i = 0; i < m_*n_; ++i) {
        top_[i] = dist(engine);
    }
}

template<>
void Matrix<double>::gen_rnd_elm() {
    
    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (uint64_t i = 0; i < m_*n_; ++i) {
        top_[i] = dist(engine);
    }
}

template<>
void Matrix<std::complex<float>>::gen_rnd_elm() {

    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    for (uint64_t i = 0; i < m_*n_; ++i) {
        top_[i] = std::complex<float> {dist(engine), dist(engine)};
    }
}

template<>
void Matrix<std::complex<double>>::gen_rnd_elm() {

    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (uint64_t i = 0; i < m_*n_; ++i) {
        top_[i] = std::complex<double> {dist(engine), dist(engine)};
    }
}

template<typename T>
T* Matrix<T>::elm(const uint32_t &ti, const uint32_t &tj) const {
    assert(ti >= 0);
    assert(ti < p_);
    assert(tj >= 0);
    assert(tj < q_);
    
    uint32_t pos = 0;
    
    pos += ti * (mb_ * n_);
    pos += (ti == p_ - 1) ? tj * (m_ % mb_) * nb_ : tj * mb_ * nb_;
    
    return &top_[pos];
}

template<typename T>
T* Matrix<T>::elm(const uint32_t &ti, const uint32_t &tj,
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

template<typename T>
void Matrix<T>::file_out(const char *fname) {
    
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

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &M){
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    m_ = M.m_;
    n_ = M.n_;
    mb_ = M.mb_;
    nb_ = M.nb_;
    p_ = M.p_;
    q_ = M.q_;
    
    for (uint64_t i = 0; i < m_*n_; i++)
        top_[i] = M.top_[i];
    
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &M) const {
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    Matrix N(m_, n_, mb_, nb_);
    for (uint64_t i=0; i<m_*n_; ++i){
        N[i] = (*this)[i] + M[i];
    }

    return N;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &M) const {
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    Matrix N(m_, n_, mb_, nb_);
    for (uint64_t i=0; i<m_*n_; ++i){
        N[i] = (*this)[i] - M[i];
    }

    return N;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &M) {
    if(m_ != M.m_ || n_ != M.n_) {
        return false;
    }
    for (uint64_t i = 0; i < m_*n_; i++) {
        if (this->top_[i] != M.top_[i]){
            return false;
        }
    }
    return true;
}

template<typename T>
T &Matrix<T>::operator[](const uint64_t &i) const {
    assert(i >= 0);
    assert(i < m_*n_);
    
    return top_[i];
}

template<typename T>
T &Matrix<T>::operator()(const uint32_t &i, const uint32_t &j) const {
    assert(i >= 0);
    assert(i < m_);
    assert(j >= 0);
    assert(j < n_);

    const uint32_t ti = i / mb_;
    const uint32_t tj = j / nb_;
    const uint32_t tp = i % mb_;
    const uint32_t tq = j % nb_;

    return top_[ (mb_*n_)*ti + tp + (mb_*nb_)*tj + (mb_*tq)];
}

/*
template<typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &ma) {
    
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
*/

template<typename T>
void Matrix<T>::zero() {

    std::memset( top_, 0, sizeof(T)*m_*n_);

}

template class Matrix<float>;
template class Matrix<double>;
template class Matrix<std::complex<float>>;
template class Matrix<std::complex<double>>;