//  Matrix.cpp

#include <iostream>
#include <fstream>
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

  mb_ = m;
  nb_ = n;

  p_ = 1;
  q_ = 1;

  top_ = new double[ m_ * n_ ];

};

Matrix::Matrix( const unsigned int m, const unsigned int n, const unsigned int ts)
{

  m_ = m;
  n_ = n;

  mb_ = ts;
  nb_ = ts;

  p_ = ( m_%mb_ == 0) ? m_ / mb_ : m_ / mb_ + 1;
  q_ = ( n_%nb_ == 0) ? n_ / nb_ : n_ / nb_ + 1;



  top_ = new double[ m_ * n_ ];

};

Matrix::Matrix( const unsigned int m, const unsigned int n, 
		const unsigned int mb, const unsigned int nb)
{

  m_ = m;
  n_ = n;

  mb_ = mb;
  nb_ = nb;

  p_ = ( m_%mb_ == 0) ? m_ / mb_ : m_ / mb_ + 1;
  q_ = ( n_%nb_ == 0) ? n_ / nb_ : n_ / nb_ + 1;

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

double* Matrix::elm(const int ti, const int tj)
{
  assert( ti >= 0 ); assert( ti < p_);
  assert( tj >= 0 ); assert( tj < q_);

  int pos = 0;

  pos += ti*(mb_*n_);
  pos += (ti == p_-1) ? tj*(m_%mb_)*nb_ : tj*mb_*nb_;

  return &top_[ pos ];

}

double* Matrix::elm(const int ti, const int tj, const int i, const int j)
{
  assert( i >= 0 ); assert( i < (ti==(p_ -1) ? m_%mb_ : mb_) );
  assert( j >= 0 ); assert( j < (tj==(q_ -1) ? n_%nb_ : nb_) );

  assert( ti >= 0 ); assert( ti < p_);
  assert( tj >= 0 ); assert( tj < q_);

  int pos = 0;

  pos += ti*(mb_*n_);
  pos += (ti == p_-1) ? tj*(m_%mb_)*nb_ : tj*mb_*nb_;
  pos += i;
  pos += (ti == p_-1) ? j*(m_%mb_) : j*mb_;

  return &top_[ pos ];
}

void Matrix::file_out( char* fname )
{
  
  ofstream matf(fname);
  if (!matf) {
    cerr << "Unable to open " << fname << endl;
    exit(1);
  }

  matf.precision(5);

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
	p += (m_%mb_ == 0) ? (j/nb_)*mb_*nb_ : (j/nb_)*(m_%mb_)*nb_;//tjの場所
	p += (m_%mb_ == 0) ? (j%nb_)*mb_ + (i%mb_) : (j%nb_)*(m_%mb_) + (i%mb_);
      }

      matf << top_[p] << " ";
    }
    matf << endl;
  }
  matf.close();
  return;
}


Matrix &Matrix::operator=( const Matrix &T )
{
	assert( m_ == T.m_ );
	assert( n_ == T.n_ );

	mb_ = T.mb_;
	nb_ = T.nb_;
	p_ = T.p_;
	q_ = T.q_;
	
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
	p += (m_%mb_ == 0) ? (j/nb_)*mb_*nb_ : (j/nb_)*(m_%mb_)*nb_;//tjの場所
	p += (m_%mb_ == 0) ? (j%nb_)*mb_ + (i%mb_) : (j%nb_)*(m_%mb_) + (i%mb_);
      }

      os << ma[p] << " ";
    }
    os << std::endl;
  }
  
  return os;
}
*/


std::ostream& operator<<(std::ostream& os, const Matrix& ma){
  
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
	p += (m_%mb_ == 0) ? (j/nb_)*mb_*nb_ : (j/nb_)*(m_%mb_)*nb_;//tjの場所
	p += (m_%mb_ == 0) ? (j%nb_)*mb_ + (i%mb_) : (j%nb_)*(m_%mb_) + (i%mb_);
      }

      os << ma[p] << " ";
    }
    os << std::endl;
  }
  
  return os;
}
