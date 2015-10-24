#include <iostream>
#include <cstdlib>
#include <cassert>
#include <iomanip>

#include "Matrix.hpp"

using namespace std;


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
