#include <iostream>
#include "LUd.hpp"
#include "matrix.hpp"

LUd::LUd(const Matrix& A)
    : d(1), lu(A)
//LUd::LUd(const double b)
{
  const double TINY=1.0e-40;
  int imax;
  double cmp, temp;
  if ( size(A,1) != size(A,2))
  {
    std::cerr << "Error, Matrix is not square" << '\n';
    exit(1);
  }
  indx = new int[size(A,2)];
  //int d=1;
  int n;
  n = size(A,2);
  //
  for (int k=1; k<= n; k++)
  {
    cmp=0.0;
    for (int i=k; i <= n ; i++)
    {
      if ( fabs(lu(i,k)) > cmp )
      {
        cmp = fabs(lu(i,k));
        imax = i;
      }
    }
    if (k!= imax)
    {
      //std::swap(lu.mData[k],lu.mData[imax]);
      lu.swaprow(k,imax);
      d=-d;
    }

    indx[k-1]=imax;
    if ( lu(k,k) == 0.0)
    {
      lu(k,k) = TINY;
      //std::cerr << "Error, Matrix is singular" << '\n';
      //exit(1);
    }
    for (int i=k+1; i<=n; i++)
    {
      temp=lu(i,k) /= lu(k,k);
      for (int j=k+1 ; j <= n ; j++)
      {
        lu(i,j) -= temp * lu(k,j);
      }
    }
  }
}

LUd::~LUd()
{
  delete[] indx;
}

void LUd::solve(Matrix& b) //const
{
  if (size(b,1) != size(lu,1))
  {
    std::cerr << "Error: LUd, rhs dimension mismatch" << '\n';
  }

  int nrowb=size(b,1);
  int ncolb=size(b,2);
  int ip;
  int startsum=0;
  double tmpsum;

  for (int i= 1; i <=nrowb; i++)
  {
    for (int k=1; k<= ncolb; k++ )
    {
      ip = indx[i-1];
      tmpsum= b(ip,k);
      b(ip,k) = b(i,k);
      if (startsum !=0)
      {
        for (int j = startsum-1; j < i; j++)
        {
          tmpsum -= lu(i,j) * b(j,k);
        }
      }
      else
      {
        if (tmpsum != 0.0)
        {
          startsum = i+1;
        }
      }
      b(i, k) = tmpsum;
    }
  }

  for (int i=nrowb; i>0; i--)
  {
    for (int k= 1; k<=ncolb; k++)
    {
      tmpsum=b(i,k);
      for (int j= i+1; j<= nrowb; j++)
      {
        tmpsum -= lu(i,j)*b(j,k);
      }
      b(i,k)=tmpsum/lu(i,i);
    }
  }
}

Matrix LUd::L()
{
  Matrix L(size(lu,1),size(lu,2));

  for (int i = 1; i<=size(lu,1); i++)
  {
    for (int j = 1; j<= i; j++)
    {
      L(i,j)=lu(i,j);
    }
    L(i,i)=1.0;
    // for (int j = i+1; j<size(lu,2); j++)
    // {
    //   L(i,j) = 0.0;
    // }
  }
  return L;
}

Matrix LUd::U()
{
  Matrix U(size(lu,1),size(lu,2));
  for (int i = 1; i<=size(lu,1); i++)
  {
    // for (int j = 0; j< i; j++)
    // {
    //   U(i,j)= 0.0;
    // }
    for (int j = i; j<=size(lu,2); j++)
    {
      U(i,j) = lu(i,j);
    }
  }
  return U;
}

Matrix LUd::P()
{
  Matrix P = eye(size(lu,1));
  for (int i=size(lu,1); i>0; i--)
  {
    P.swaprow(i,indx[i-1]);
    //std::swap(P.mData[i],P.mData[indx[i]]);
  }
  return P;
}

void LUd::det(double &deter) const
{
  deter=d;
  int n=size(lu,1);
  for (int i = 1; i<= n; i++)
  {
    deter*= lu(i,i);
  }
}

double det(const Matrix& A)
{
  LUd lud(A);
  double deter;
  lud.det(deter);
  return deter;
}
