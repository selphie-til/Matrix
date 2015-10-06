//  Matrix.cpp


#include <iostream>
#include <cstdlib>
#include <cassert>
#include <random>

#include "Matrix.hpp"

using namespace std;

Matrix::Matrix(){
  m_ = 0;
  n_ = 0;
  top_ = nullptr;

}

Matrix::Matrix( const unsigned int m, const unsigned int n){

  m_ = m;
  n_ = n;

  top_ = new double[ m_ * n_ ];

};

Matrix::~Matrix(){

  delete[] top_;
  top_ = nullptr;

}

void Matrix::gen_rnd_elm(){

  //乱数生成器
  random_device rdev;
  //ランダムなシードの設定
  mt19937 engine(rdev());
  uniform_real_distribution<> dist(0.0,1.0);

  for(int i=0; i< m_*n_; i++)
    top_[i] = dist(engine);

}

double* Matrix::elm(const int i, const int j)
{
  assert( i >= 0 ); assert( i < m_);
  assert( j >= 0 ); assert( j < n_);

  return &top_[ i + j * m_];
}

Matrix &Matrix::operator=( const Matrix &T )
{
	assert( m_ == T.m_ );
	assert( n_ == T.n_ );

	for (unsigned int i = 0; i < m_ * n_; i++)
		top_[i] = T.top_[i];

	return *this;
}

double &Matrix::operator[]( const unsigned int i ) const
{
	assert( i >= 0 );
	assert( i < m_ * n_ );

	return top_[ i ];
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */

double &Matrix::operator()( const unsigned int i, const unsigned int j ) const
{
	assert( i >= 0 );  assert( i < m_ );
	assert( j >= 0 );  assert( j < n_ );

	return top_[ i + j * m_ ];
}
