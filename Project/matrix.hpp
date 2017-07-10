#ifndef MATRIXDEF
#define MATRIXDEF

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>  // common maths functions
#include <initializer_list>
#include <algorithm>
#include "Exception.hpp"//  This class throws errors using the class "error"

//class LUd;

class Matrix
{
  //friend class LUd;
private:
   // member variables
   double** mData;   // data stored in Matrix
   int mNRow;       // number of rows of Matrix
   int mNCol;       // number of columns of Matrix
   int mRowRes;
   int mColRes;
   bool mCopy;      // states whether the matrix refers to anothers elements, relevant for destructor
   Matrix(int nrRow, int nrCol, double**& pt_to_cont);
   Matrix(int nrRow, int nrCol, char empty);

public:
///// constructor
  // No default constructor
  // overridden copy constructor
  Matrix(const Matrix& m1);
  Matrix(int nrRow, int nrCol, double val = 0);
  Matrix(int n, double val = 0);
  Matrix(int nrRow, int nrCol,std::initializer_list<double> l);
  Matrix(std::string filename);
  Matrix(int nrRow, int nrCol, bool name, int resnrRow, int resnrCol);
  // named constructor
  static Matrix MatrixWithRes(int nrRow, int nrCol, int resnrRow, int resnrCol);
  // destructor
  ~Matrix();

/////indexing
  double& operator()(int i, int j);
  const double& operator()(int i, int j) const;
  double& operator()(int i);
  const double& operator()(int i) const;
  Matrix operator()(int i, char c);
  Matrix operator()(char c, int i);
  Matrix operator()(std::initializer_list<int> range1, std::initializer_list<int> range2);

/////assignment
  Matrix& operator=(const Matrix& m);
  Matrix& operator+=(const Matrix& m);
  Matrix& operator-=(const Matrix& m);
  Matrix& operator*=(Matrix m);
  Matrix& operator*=(const double& a);
  //Matrix& operator/=(Matrix b);
  Matrix& operator/=(const double& a);

///// storage structure altering members
  // transpose
  void T();

  // append
  void append(const Matrix& m);

  // swap rows
  void swaprow(int i, int j);

  // givens rotation
  double givens(int i, int j);

  void givens(int i, int j, double theta);

  void gaussj(Matrix& b);

  // All "friend" external operators and functions are declared as friend inside the class (here)
  // but their actual prototype definitions occur outside the class.
////// Binary operators
  //friend Matrix operator+(const Matrix& m1, const Matrix& m2);
  //friend Matrix operator-(const Matrix& m1, const Matrix& m2);
  //friend Matrix operator*(Matrix m1, const Matrix& m2);
  //friend Matrix operator*(const Matrix& m, const double& a);
  //friend Matrix operator*(const double& a, const Matrix& m);
  //friend Matrix operator/(const Matrix& m, const double& a);
  //friend Matrix operator/(const Matrix& m1, const Matrix& m1);
  // Unary operator
  //TODO friend Matrix operator-(const Matrix& v);
  // TODO write .T method for transpose

   //other operators


/////output
  //friend std::ostream& operator<<(std::ostream& output, const Matrix& m);

//// member methods


   //double T(i,j);
   //norm (as a member method)
   double norm(int p=2) const;
   // functions that are friends
   //TODO friend
   //friend void transpose(Matrix& A);
   friend int size(const Matrix& m, int dir);
   friend int reservedsize(const Matrix& m, int dir);
  //  friend void gaussj(Matrix & m, Matrix& b);
   friend bool iscopy(const Matrix& m);
  //  friend void badmult(Matrix )
   //friend Matrix inv(const Matrix & m);
   //friend Matrix eye(int rows, int cols);
   //friend Matrix eye(int dim);
   //friend Matrix diag( Matrix& m, int offs);
};

std::ostream& operator<<(std::ostream& output, const Matrix& m);
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m, const double& a);
Matrix operator*(const double& a, const Matrix& m);
Matrix operator/(const Matrix& m, const double& a);
Matrix operator/(const double& a, const Matrix& m);
Matrix operator/(const Matrix& A, Matrix b);
double norm(Matrix& v, int p=2);
//Matrix operator/(const Matrix& v, const double& a);
//Matrix operator/(const Matrix& m1, const Matrix& m1);

Matrix transpose(const Matrix& A);
int size(const Matrix& m, int dir);
int reservedsize(const Matrix& m, int dir);
void gaussj(Matrix & m, Matrix& b);
Matrix inv(Matrix & m);
Matrix eye(int rows, int cols);
Matrix eye(int dim);
Matrix diag(const Matrix& m, int offs = 0);

Matrix GMRES(const Matrix & A, const Matrix& b,const  Matrix& x0, double Tol, std::ostream &file = std::cout);

Matrix GMRES2(const Matrix & A, const Matrix& b,const  Matrix& x0, double Tol, std::ostream &file = std::cout);

double mat2d(const Matrix& M);

double scaprod(const Matrix& m1, const Matrix& m2);

bool iscopy(const Matrix& m);

#endif
