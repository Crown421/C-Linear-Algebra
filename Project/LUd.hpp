#ifndef LUDDEF
#define LUDDEF

#include <algorithm>
#include "matrix.hpp"

class Matrix;

class LUd
{
private:
  int d;
  Matrix lu;
  int* indx;

public:
  ~LUd();
  LUd(const Matrix& A);
  void solve(Matrix& b) ;// const;
  Matrix L();
  Matrix U();
  Matrix P();
  void det(double &deter) const;
};

double det(const Matrix & A);

#endif
