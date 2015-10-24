//  Matrix.hpp

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

class Matrix {

protected:
  double* top_;		// pointer to the matrix
  unsigned int m_;	// number of lows of the matrix or lda
  unsigned int n_;	// number of columns of the matrix
  unsigned int mb_;	// number of lows of the tile
  unsigned int nb_;	// number of columns of the tile
  unsigned int p_;      // number of low tiles
  unsigned int q_;      // number of column tiles


public:
  // Default constructor
  Matrix();
  // Constructor
  Matrix( const unsigned int m, const unsigned int n );
  Matrix( const unsigned int m, const unsigned int n, const unsigned int ts );
  Matrix( const unsigned int m, const unsigned int n,
	  const unsigned int mb, const unsigned int nb);

  // Destructor
  ~Matrix();

  // Getters
  double* top() { return top_; }
  unsigned int m() const { return m_; }
  unsigned int n() const { return n_; }
  unsigned int mb() const { return mb_; }
  unsigned int mb(const int ti,const int tj){ return (ti == p_-1 ? m_%mb_: mb_ ); }
  unsigned int nb() const { return nb_; }
  unsigned int nb(const int ti,const int tj){ return (tj == q_-1 ? n_%nb_: nb_ ); }
  unsigned int p() const { return p_; }
  unsigned int q() const { return q_; }
  
  // Assign randam numbers to matrix elements
  void gen_rnd_elm();
  //the pointer of the (i,j) element
  double* elm(const int, const int);
  //the pointer of the (i,j) element (ti,tj)
  double* elm(const int, const int, const int, const int);

  void file_out( char* fname );
  
  // Operator overload
  Matrix &operator=( const Matrix& T );
  double &operator[]( const unsigned int i ) const;
  double &operator()( const unsigned int i, const unsigned int j ) const;

  /*
  template<typename CharT, typename Traits>
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, const Matrix& ma);
  */
  friend std::ostream& operator<<(std::ostream& os, const Matrix& ma);
};


#endif /* MATRIX_HPP_ */

