Matrix& Matrix::operator=(const Matrix& m)
{
  if (m.mNRow > mNRow || m.mNCol > mNCol)
    {
      throw Exception("length mismatch",
		  "matrix assignment operator - matrix to be copied is too big");
    }
  else
    {
      for (int i=0; i<m.mNRow; i++)
	     {
         for (int j=0; j<m.mNCol; j++)
         {
           mData[i][j] = m.mData[i][j];
         }
	     }
       if (m.mNRow < mNRow)
         {
           for (int i=m.mNRow; i<mNRow; i++)
     	     {
             for (int j=0; j<mNCol; j++)
             {
               mData[i][j] = 0.0;
             }
     	     }
           std::cerr << "matrix assignment operator - copied matrix had too few rows";
           std::cerr << " and has been extended with zero rows\n";
         }
       if (m.mNCol < mNCol)
         {
           for (int i=0; i<mNRow; i++)
           {
             for (int j=m.mNCol; j<mNCol; j++)
             {
               mData[i][j] = 0.0;
             }
           }
           std::cerr << "matrix assignment operator - copied matrix had too few columns";
           std::cerr << " and has been extended with zero columns\n";
         }
    }

    return *this;
}





/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// scraped due to row/ column assignment
Matrix& Matrix::operator=(const Matrix& m)
{
  if (m.mNRow != mNRow || m.mNCol != mNCol)
    {
      std::cerr << "assignment changed matrix size" << '\n';
      double** tmpData;
      if ( (m.mNCol == mNCol) && (m.mNRow != mNRow) )
      {
        tmpData=new double* [m.mNRow];
        //minrow=std::min(m.mNRow, mNRow);
        for (int i = 0; i< m.mNRow; i++)
        {
          tmpData[i]=m.mData[i];
        }
        for (int i = m.mNRow; i< mNRow ; i++ )
        {
          delete[] mData[i];
          std::cout << "Deleted sth in assignment" << '\n';
        }
      }
      else
      {
        tmpData=new double* [m.mNRow];
        //minrow=std::min(m.mNRow, mNRow);
        for (int i = 0; i< m.mNRow; i++)
        {
          tmpData[i] = new double [m.mNCol];
          for (int j = 0; j< m.mNCol; j++)
          {
            tmpData[i][j]=m.mData[i][j];
          }
        }
        for (int i = m.mNRow; i< mNRow ; i++ )
        {
          delete[] mData[i];
        }
      }
      mData=tmpData;
    }
  else
    {
      for (int i=0; i<m.mNRow; i++)
	     {
         for (int j=0; j<m.mNCol; j++)
         {
           mData[i][j] = m.mData[i][j];
         }
	     }

    }

    return *this;
}

Matrix operator+(const Matrix& m1, const Matrix& m2)
{
  if (m1.mNRow != m2.mNRow || m1.mNCol != m2.mNCol)
  {
    throw Exception("dimension mismatch",
    "addition is not defined for matrices of different size");
  }

  Matrix M(m1.mNRow,m1.mNCol,'e');

  for (int i=0; i< m1.mNRow; i++)
  {
    for (int j=0; j<m1.mNCol; j++)
    {
      M.mData[i][j] = m1.mData[i][j] + m2.mData[i][j];
    }
  }
  return M;
}



/// multiplication
Matrix operator*(const Matrix& m1, const Matrix& m2)
{
  if ( m1.mNCol != m2.mNRow )
  {
    throw Exception("dimension mismatch",
    "multiplication is not defined for matrices of mismatching row/columns");
  }

  Matrix M(m1.mNRow,m2.mNCol,'e');

  Matrix m2t=m2;
  m2t.T();

  for (int i=0; i< m1.mNRow; i++)
  {
    for (int j=0; j<m2.mNCol; j++)
    {
      //M.mData[i][j] = m1.mData[i][0] * m2.mData[0][j];
      M.mData[i][j] = m1.mData[i][0] * m2t.mData[j][0];
      for (int k = 1; k<m1.mNCol ;k++)
      {
          //M.mData[i][j] = M.mData[i][j] + m1.mData[i][k] * m2.mData[k][j];
          M.mData[i][j] = M.mData[i][j] + m1.mData[i][k] * m2t.mData[j][k];
      }
    }
  }
  return M;
}


Matrix y = H({1,l},{1,l})/(res[0]*QT({1,l},{1,1}));


std::cout << "lu decom \n Hhat" << H({1,l},{1,l}) << " res* QT" << res[0]*QT({1,l},{1,1}) << "Solution " <<  y << '\n';

Matrix Hhat= H({1,l},{1,l});
Matrix b2 = (res[0]*QT({1,l},{1,1}));
std::cout << "Gaus jordan \n" << "Hhat" << Hhat << " res* QT" << b2 << '\n';
gaussj(Hhat,b2);
std::cout << "Solution " << b2 << " solution test " << Hhat*b2 << "inverse test" << Hhat * H({1,l},{1,l}) << '\n';
