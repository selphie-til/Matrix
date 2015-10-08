#include <iostream>
#include <cstdlib>
#include <cassert>
#include <iomanip>

#include "Matrix.hpp"

using namespace std;

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
    os << endl;
  }
  
  return os;
}

int main(const int argc, const char* argv[])
{

  if( argc < 3 ){
    cerr << "Usage: ./Matrix_test [M] [N] [ts]\n";
    return -1;
  }

  const int M = atoi(argv[1]);
  const int N = atoi(argv[2]);
  const int ts = atoi(argv[3]);

  Matrix *A = new Matrix(M,N,ts);


  A->gen_rnd_elm();

  cout << setprecision(5);

  cout << *A << endl;

  for(int i=0; i<M*N; i++){
      cout << (*A)[i] << endl;;
  }
  char filename[] = "test.dat";
  A->file_out(filename);

  return 0;

}
