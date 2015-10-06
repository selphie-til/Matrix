#include <iostream>
#include <cstdlib>
#include <cassert>
#include <iomanip>

#include "Matrix.hpp"

using namespace std;

int main(const int argc, const char* argv[])
{

  if( argc < 3 ){
    cerr << "Usage: ./Matrix_test [M] [N]\n";
    return -1;
  }

  const int M = atoi(argv[1]);
  const int N = atoi(argv[2]);

  Matrix *A = new Matrix(M,N);

  A->gen_rnd_elm();

  cout << setprecision(5);

  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++)
      cout << *(A->elm(i,j)) << ' ';
    cout << endl;
  }

  return 0;

}
