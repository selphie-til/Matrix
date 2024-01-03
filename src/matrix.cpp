//  Matrix.cpp

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <random>
#include <complex>

#include "matrix.hpp"

template<typename T>
Matrix<T>::~Matrix() {
    this->top_ = nullptr;
}

template<typename T>
void Matrix<T>::gen_rnd_elm() {
    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<T> dist(0.0, 1.0);

    for (uint64_t i = 0; i < m_*n_; ++i) {
        ((this->top_).get())[i] = dist(engine);
    }
}

template<typename T>
T* Matrix<T>::elm(const uint32_t &ti, const uint32_t &tj) const {
    assert(ti >= 0);
    assert(ti < (this->p_));
    assert(tj >= 0);
    assert(tj < (this->q_));

    uint32_t pos = 0;

    pos += ti * ((this->mb_) * (this->n_));
    pos += (ti == (this->p_ - 1)) ? tj * (this->m_ % this->mb_) * (this->nb_) : tj * this->mb_ * this->nb_;

    return &(top_.get()[pos]);
}

template<typename T>
T* Matrix<T>::elm(const uint32_t &ti, const uint32_t &tj,
                  const uint32_t &i, const uint32_t &j) const {
    assert(i >= 0);
    assert(i < (ti == (this->p_ - 1) ? (this->m_) % (this->mb_) : this->mb_));
    assert(j >= 0);
    assert(j < (tj == (this->q_ - 1) ? (this->n_) % (this->nb_) : this->nb_));

    assert(ti >= 0);
    assert(ti < this->p_);
    assert(tj >= 0);
    assert(tj < this->q_);

    uint64_t pos = 0;

    pos += ti * (this->mb_ * this->n_);
    pos += (ti == (this->p_ - 1)) ? tj * (this->m_ % this->mb_) * this->nb_ : tj * this->mb_ * this->nb_;
    pos += i;
    pos += (ti == (this->p_ - 1)) ? j * (this->m_ % this->mb_) : j * this->mb_;

    return &(top_.get()[pos]);
}

template<typename T>
void Matrix<T>::file_out(const char *fname) {

    std::ofstream matf(fname);
    if (!matf) {
        std::cerr << "Unable to open " << fname << std::endl;
        exit(1);
    }

    matf.precision(5);

    for (uint32_t i = 0; i < this->m_; i++) {
        for (uint32_t j = 0; j < this->n_; j++) {
            uint64_t p = 0;

            // (i / mb_) ti
            if ((i / this->mb_) != this->p_ - 1) {
                p += (i / this->mb_) * this->mb_ * this->n_;  // tiの場所
                p += (j / this->nb_) * this->mb_ * this->nb_; // tjの場所

                p += (j % this->nb_) * this->mb_ + (i % this->mb_); // i,jの場所
            } else {
                // ここで最終行の時と差別化
                p += (this->p_ - 1) * this->mb_ * this->n_; // tiの場所
                p += ((this->m_ % this->mb_) == 0) ? (j / this->nb_) * this->mb_ * this->nb_
                                     : (j / this->nb_) * (this->m_ % this->mb_) * this->nb_; // tjの場所
                p += ((this->m_ % this->mb_) == 0) ? (j % this->nb_) * this->mb_ + (i % this->mb_)
                                     : (j % this->nb_) * (this->m_ % this->mb_) + (i % this->mb_);
            }

            matf << (this->top_)[p] << " ";
        }
        matf << std::endl;
    }
    matf.close();
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &M){
    assert(this->m_ == M.m_);
    assert(this->n_ == M.n_);

    this->m_ = M.m_;
    this->n_ = M.n_;
    this->mb_ = M.mb_;
    this->nb_ = M.nb_;
    this->p_ = M.p_;
    this->q_ = M.q_;

    for (uint64_t i = 0; i < (this->m_)*(this->n_); i++)
        (this->top_).get()[i] = (M.top_).get()[i];

    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::addElements(const Matrix<T> &M) const {
    assert(this->m_ == M.m_);
    assert(this->n_ == M.n_);

    Matrix N(this->m_, this->n_, this->mb_, this->nb_);
    uint64_t totalElements = (this->m_) * (this->n_);
    for( uint64_t i=0; i<totalElements; ++i){
        N[i] = (*this)[i] + M[i];
    }

    return N;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &M) const {
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    return this->addElements(M);
}

template<typename T>
Matrix<T> Matrix<T>::minusElements(const Matrix<T> &M) const{
    assert(this->m_ == M.m_);
    assert(this->n_ == M.n_);

    Matrix N( this->m_, this->n_, this->mb_, this->nb_);
    uint64_t totalElements = (this->m_) * (this->n_);
    for( uint64_t i=0; i<totalElements; ++i){
        N[i] = (*this)[i] - M[i];
    }

    return N;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &M) const {
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    return this->minusElements(M);
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
    assert(i < (this->m_)*(this->n_));

    return top_[i];
}

template<typename T>
T &Matrix<T>::operator()(const uint32_t &i, const uint32_t &j) const {
    assert(i >= 0);
    assert(i < (this->m_));
    assert(j >= 0);
    assert(j < (this->n_));

    const uint32_t ti = i / this->mb_;
    const uint32_t tj = j / this->nb_;
    const uint32_t tp = i % this->mb_;
    const uint32_t tq = j % this->nb_;

    return top_[ (this->mb_*this->n_)*ti + tp + (this->mb_*this->nb_)*tj + (this->mb_*tq)];
}

template<typename T>
std::ostream &operator<< (std::ostream &os, const Matrix<T> &ma) {

    uint32_t m_ = ma.m();
    uint32_t n_ = ma.n();
    uint32_t mb_ = ma.mb();
    uint32_t nb_ = ma.nb();
    uint32_t p_ = ma.p();

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

            os << ma[p] << " ";
        }
        os << std::endl;
    }

    return os;
}

template<typename T>
void Matrix<T>::zero() {

    std::fill(top_.get(), top_.get() + (this->m_) + (this->n_), 0);

}

template class Matrix<float>;
template class Matrix<double>;