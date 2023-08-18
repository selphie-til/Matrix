//  Matrix.hpp
#pragma once

#include <iostream>
#include <cstdint>

template <typename T>
class Matrix {
    
protected:
    T *top_;     // pointer to the matrix
    uint32_t m_;  // number of lows of the matrix or lda
    uint32_t n_;  // number of columns of the matrix
    uint32_t mb_{}; // number of lows of the tile
    uint32_t nb_{}; // number of columns of the tile
    uint32_t p_{};  // number of low tiles
    uint32_t q_{};  // number of column tiles
    
public:
    // Default constructor
    Matrix():m_{0}, n_{0},
             mb_{0}, nb_{0},
             p_{0}, q_{0}, top_{nullptr} {};
    // Constructor
    Matrix(const uint32_t &m, const uint32_t &n):m_{m}, n_{n},
                                                 mb_{m}, nb_{n},
                                                 p_{1}, q_{1}, top_{new T[m*n]} {};
    Matrix(const uint32_t &m, const uint32_t &n, const uint32_t &ts):m_{m}, n_{n},
                                                                     mb_{ts}, nb_{ts},
                                                                     p_{(m % ts == 0) ? m / ts : m / ts + 1},
                                                                     q_{(n % ts == 0) ? n / ts : n / ts + 1},
                                                                     top_{new T[m*n]} {};
    Matrix(const uint32_t &m, const uint32_t &n,
           const uint32_t &mb, const uint32_t &nb):m_{m}, n_{n},
                                                   mb_{mb}, nb_{nb},
                                                   p_{(m % mb == 0) ? m / mb : m / mb + 1},
                                                   q_{(n % nb == 0) ? n / nb : n / nb + 1},
                                                   top_{new T[m*n]} {};

    // Destructor
    ~Matrix();
    
    // Getters
    T *top() { return top_; }
    uint32_t m() const { return m_; }
    uint32_t n() const { return n_; }
    uint32_t mb() const { return mb_; }
    uint32_t mb(const int ti, const int tj) const {
        return (ti == p_ - 1 ? m_ % mb_ : mb_);
    }
    uint32_t nb() const { return nb_; }
    uint32_t nb(const int ti, const int tj) const {
        return (tj == q_ - 1 ? n_ % nb_ : nb_);
    }
    uint32_t p() const { return p_; }
    uint32_t q() const { return q_; }
    
    // Assign random numbers to matrix elements
    void gen_rnd_elm();
    // the pointer of the (i,j) element
    T *elm(const uint32_t &ti, const uint32_t &tj) const;
    // the pointer of the (i,j) element (ti,tj)
    T *elm(const uint32_t &ti, const uint32_t &tj,
           const uint32_t &i, const uint32_t &j) const;
    
    void file_out(const char *fname);

    void zero();
    
    // Operator overload
    Matrix &operator=(const Matrix &M);
    Matrix operator+(const Matrix &M) const;
    Matrix operator-(const Matrix &M) const;
    bool operator==(const Matrix &M);
    T &operator[](const uint64_t &i) const;
    T &operator()(const uint32_t &i, const uint32_t &j) const;

    /*
    friend std::ostream &operator<<(std::ostream &os, const Matrix<T> &ma);
     */
};
