#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "matrix.hpp"
#include "LUd.hpp"

int main()
{
  // { // test SAFETY flag
  //   Matrix m(2,2);
  //   double a = m(3,1);
  //   std::cout << "print SAFETY test" << a << '\n';
  // }
  {
  Matrix m(2,2); //test zero initialization
  m(1,2) = 1.0;   // test single element assignment
  std::cout << "first element = " << m(1,1) << '\n';
  m(2,1) = m(3);   // test single element assignment
  std::cout << "first element = " << m(2,1) << '\n';
  std::cout << "Init matrix m = " << m << '\n'; // test overloaded output operator
  }

  {
  Matrix m(2); //test zero initialization
  std::cout << "Init matrix m with one index " << m << '\n'; // test overloaded output operator
  }


  Matrix m(2,1.0); //test zero initialization
  std::cout << "Init matrix m with one index and 1 everywhere (basically ones) " << m << '\n'; // test overloaded output operator


  Matrix m2(m); // test copy constructor

  std::cout << "Matrix copy = " << m2 << '\n';

  Matrix m4(3,3,
    {0.0, 1.0 , 2.0,
     4.0, 5.0, 6.0,
     7.0, 8.0}
    );

  std::cout << "List initialized, one element to few, so last one zero,  m4 = " << m4 << '\n';

  std::cout << "m4(3) = " << m4(3) << '\n'; // test single index access

std::cout << "/////////////////////////////////////////////////////////////" << '\n';
  Matrix m5(m4(1,':')); // copy constructor with overloaded () for row access, giving a matrix class

  std::cout << "First row of m4, m5= " << m5 << '\n';

  std::cout << "m4(1,:) = \n" << m4(1,':') << '\n';

  std::cout << "m4(:,1) = \n" << m4(':',1) << '\n';

  std::cout << "m4(1,:) = " << m4(1,':') << '\n';

  Matrix m6=m4(1,':'); // initialization via conversion constr, essentially copy constr

  std::cout << "First row of m4, m6= " << m5 << '\n';

  Matrix m7(1,3,
        {9.0, 9.0, 9.0}
  );


  std::cout << "Replacement row, m7= " << m7 << '\n';

  Matrix m8(3,3);
  m8=m4;
  std::cout << "Copy of m4, m8= " << m8 << '\n';

  m8(1,':')=m7;

  std::cout << "Replaced first row of m8 with m7, m8= " << m8 << '\n';

  std::cout << "Trying to assign row vector to column m8, changes result " << '\n';
  try
  {
    m8(':',1)=m7;
    std::cout << "Replacement column, m7= " << m8 << '\n';
  }
  catch (Exception &ex)
  {
      ex.DebugPrint();
  }

  Matrix m9(3,1,
        {1.0, 1.0, 1.0}
  );

  std::cout << "Replacement column, m9= " << m9 << '\n';

  m8(':',1)=m9;
  std::cout << "Replaced first column of m8 with m9, m8= " << m8 << '\n';

  m8.append(m9);
  std::cout << "Appended m9 to m8, m8= " << m8 << '\n';

  Matrix m10(4,1,
        {1.0, 1.0, 1.0,1.0}
  );

  std::cout << "Longer column column, m10= " << m10 << '\n';

  m8.append(m10);
  std::cout << "Appended m10 to m8, m8= " << m8 << '\n';

  m8.append(m9);
  std::cout << "Appended m9 to m8 again, m8= " << m8 << '\n';

  Matrix m11 = m4+m4;
  std::cout << "m4+m4= m11 via init =" << m11 << '\n';

  Matrix m12(3,3);
  m12 = m4-m4;
  std::cout << "m4-m4= m12 via pre init =" << m12 << '\n';

  Matrix m13 = m4*m4;
  std::cout << "m4*m4= m11 via init =" << m13 << '\n';

  //Matrix m13(1,1,{1.0});
  //m12(1,2)=m13;

  Matrix m14(3,3,
  { 2.0, 7.5, 0.0,
    0.0, 3.0, 0.0,
    0.0, 0.5, 5.0});
  Matrix m14inv = m14;
  Matrix b14 = m9;

  std::cout << "m14, size (" << size(m14,1) << ", " << size(m14,2) << ") and test (size(.,3))" << size(m14,3) << " gives " << m14 << '\n';
  m14inv.gaussj(b14);

  std::cout << "test gaussj output" << '\n';

  std::cout << "m4 inverse m14inv" << m14inv << '\n';
  std::cout << "solution b14" << b14 << '\n';

  std::cout << "test solution m14* b14 " << m14* b14 << '\n';
  std::cout << "test solution again" << m14* b14 << '\n';
  std::cout << "test inverse" << m14inv* m14 << '\n';

  Matrix m14inv2 = inv(m14);
  std::cout << "m4 inverse via separate fun" << m14inv2 << '\n';

  try
    {
      std::cout << "unit matrix (3,4) (two args)" << eye(3,4) << '\n';
    }
  catch (Exception &ex)
  {
    ex.DebugPrint();
  }

  std::cout << "unit matrix (3,3) (one arg)" << eye(3) << '\n';

  std::cout << "m4 repeat" << m4 << '\n';
  LUd lud(m4);
  std::cout << "LU decomp done" << '\n';
  Matrix l=lud.L();
  std::cout << "L of m4" << l << '\n';
  Matrix u=lud.U();
  std::cout << "U of m4" << u << '\n';
  Matrix p=lud.P();
  std::cout << "P of m4" << p << '\n';

  Matrix blu(3,4,
      {1.0, 2.0, 3.0, 4.0,
       1.0, 2.0, 3.0, 4.0,
       1.0, 2.0, 3.0, 4.0});

  Matrix blubs = blu;
  std::cout << "test out, solve for rhs " << blu << '\n';
  lud.solve(blu);

  std::cout << "solution blu" << blu << '\n';
  std::cout << "test solution m4* blu " << m4* blu << '\n';

  //Matrix tmp = l*u;
  std::cout << "Test decomp just l*u" << l*u << '\n';
  //Matrix tmp2=p*tmp;
  std::cout << "Test decomp total p*l*u" << p*l*u << '\n';

  std::cout << "Solution via backslash" << m4/blubs << '\n';

  std::cout << "Determinant of m4" << det(m4) << '\n';

  Matrix diagm4 = diag(m4);
  std::cout << "Diag of m4" << diagm4 << '\n';

  Matrix diagm4m1 = diag(m4,-1);
  std::cout << "-1 Diag of m4" << diagm4m1 << '\n';

  Matrix diagm41 = diag(m4,1);
  std::cout << "1 Diag of m4" << diagm41 << '\n';

  //Matrix diagm44 = diag(m4,4);
  //std::cout << "4 Diag of m4" << diagm44 << '\n';

  Matrix Dm7 = diag(m7);
  std::cout << "Diagonal matrix" << Dm7 << '\n';

  Matrix Dm72 = diag(m7,1);
  std::cout << "Diagonal matrix with offset +1" << Dm72 << '\n';

  Matrix Dm73 = diag(m7,-1);
  std::cout << "Diagonal matrix with offset -1" << Dm73 << '\n';

  {
    Matrix A(2,4,
    {1.0,2.0,3.0,4.0,
     5.0,6.0,7.0,8.0});
    std::cout << "To be transposed, A=" << A << '\n';
    Matrix AT2= transpose(A);
    std::cout << "Transposed via transposed, A=" << AT2 << '\n';
    A.T();
    std::cout << "Transposed via member, A=" << A << '\n';
  }

  Matrix msubm(m4({1,2},{1,2})); // copy constructor with overloaded () for row access, giving a matrix class

  std::cout << "First row of m4, m5= " << msubm << '\n';

  std::cout << "m4(1,:) = \n" << m4({1,2},{1,2}) << '\n';

  Matrix m16=m4;
  m16({1,2},{1,2})=eye(2);
  std::cout << "Take copy of m4=" << m4 <<"\n and insert 2,2 identity matrix at (1,1)" << m16 << '\n';

  Matrix m17(3,3,
          { 1.0, 2.0 , 3.0,
            4.0, 5.0, 6.0,
            0.0, 7.0,8.0 });
  std::cout << "copy via =" << '\n';
  Matrix m17c = m17;
  std::cout << "copy via ()" << '\n';
  Matrix m17c2(m17);
  Matrix QT = eye(3);
  double tmpcosth;
  std::cout << "Test matrix for gives" << m17 << '\n';
  tmpcosth=m17.givens(1,2);
  QT.givens(1,2,tmpcosth);
  std::cout << "Givens (1,2) applied" << m17 << '\n';
  tmpcosth=m17.givens(2,3);
  QT.givens(2,3,tmpcosth);
  std::cout << "Givens (2,3) applied" << m17 << '\n';

  std::cout << "QT" << QT << '\n';

  std::cout << "QT * m17c" << QT* m17c << '\n';

  std::cout << "Q * R (m17)" << transpose(QT) * m17 << '\n';

  std::cout << "m10" << m10 << '\n';
  std::cout << "norm of m10" << norm(m10) << '\n';


  Matrix M(4,4,
          { 1.0, 2.0 , 3.0, 5.0,
            4.0, 5.0, 6.0, 6.0,
            0.0, 7.0,8.0, 9.0,
            10.0, 4.0, 6.0, 3.0});
  Matrix b(4,1,{1.0,1.0, 1.0, 1.0});
  Matrix x0(4,1);
  // Matrix M(5,5,
  //         { 0.0, 0.0 , 0.0, 0.0, 1.0 ,
  //           1.0, 0.0, 0.0, 0.0, 0.0,
  //           0.0, 1.0, 0.0, 0.0, 0.0,
  //           0.0, 0.0, 1.0, 0.0, 0.0,
  //           0.0, 0.0, 0.0, 1.0, 0.0});
  // Matrix b(5,1,{1.0,0.0, 0.0, 0.0, 0.0});
  // Matrix x0(5,1);

  // std::cout << "//////////////////////////////////////////////////////////" << '\n';
  // std::cout << " 4*b = " << 4*M({1,2},{1,2}) << " b*4 = " << M({1,2},{1,2}) +  M({1,2},{1,2}) << " b = " << M << '\n';
  // std::cout << "copy sum \n" << M({3,4},{1,2}) +  M({3,4},{1,2})<< '\n';
  // std::cout << "noncopy sum\n" << M+M << '\n';
  // std::cout << "copy sub \n" << M({1,2},{1,2}) -  M({1,2},{1,2})<< '\n';
  // std::cout << "noncopy sub\n" << M-M << '\n';
  // std::cout << "control " << M << '\n';
  // std::cout << "//////////////////////////////////////////////////////////" << '\n';
  std::cout << "start GMRES" << '\n';
  Matrix r =   M * b;
  // std::cout << "test r outside of GMRES " << r << '\n';
  Matrix x = GMRES(M,b,x0,1e-10);

  std::cout << "reference: b= " << b << " solution" << x << "A*x" << M*x << '\n';

  std::cout << "start GMRES2" << '\n';
  Matrix xG2 = GMRES2(M,b,x0,1e-10);

  std::cout << "reference: b= " << b << " solution xG2" << xG2 << "A*x" << M*xG2 << '\n';

  Matrix x2 = M/b;

  std::cout <<  "solution: x2= " << x2 << "A*x" << M*x2 << '\n';

  exit(0);
}
