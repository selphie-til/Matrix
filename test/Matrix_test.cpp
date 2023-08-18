#include <iostream>
#include <cstdlib>
#include <cassert>
#include <iomanip>

#include "matrix.hpp"

using namespace std;


int main(const int argc, const char* argv[])
{
    
    if( argc < 3 ){
	cerr << "Usage: ./Matrix_test [M] [N] [ts]\n";
	return -1;
    }
    
    const uint32_t M = atoi(argv[1]);
    const uint32_t N = atoi(argv[2]);
    const uint32_t ts = atoi(argv[3]);
    
    auto *A = new Matrix<double>(M,N,ts);
    auto *B = new Matrix<double>(M,N,ts);
    
    A->gen_rnd_elm();
    *B = *A;

    if( *A == *B ){
	cout << "True" << endl;
    }else{
	cout << "False" << endl;
    }
    
    cout << setprecision(5);
    
    cout << *A << endl;
    
    for(int i=0; i<M*N; i++){
        cout << (*A)[i] << endl;;
    }
    char filename[] = "test.dat";
    A->file_out(filename);
    
    return 0;
}
