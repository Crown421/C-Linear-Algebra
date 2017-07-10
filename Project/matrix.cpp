#include "matrix.hpp"
#include "LUd.hpp"

//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// constructors
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///// copy constructor TODO add shared pointer magic
Matrix::Matrix(const Matrix& m1)
{
  mNRow = m1.mNRow;
  mNCol = m1.mNCol;
  mRowRes = mNRow;
  mColRes = mNCol;
  // std::cout << "deb mRow= "<< mNRow << " mCol= " << mNCol << " resRow= " << mRowRes << " rescol= " << mColRes << '\n';
  mCopy = false;
  mData = new double*[mNRow];
  for (int i=0; i<mNRow; i++)
  {
    mData[i] =  new double[mNCol];
    for (int j=0; j<mNCol; j++){
      mData[i][j] = m1.mData[i][j];
    }
  }
  // std::cout << "Constructor" << '\n';
}

///// constructor to initialize given size with zeros
Matrix::Matrix(int nrRow, int nrCol, double val)
{
  #ifndef SAFETYOFF
  if (nrRow < 1 || nrCol < 1)
  {
    std::cerr << "Constructor error, dimension must be positive" << '\n';
    exit(1);
  }
  #endif
  mNRow = mRowRes = nrRow;
  mNCol = mColRes = nrCol;
  mCopy = false;
  mData = new double*[mNRow];
  for (int i=0; i<mNRow; i++)
  {
    mData[i] =  new double[mNCol];
    std::fill_n(mData[i],mNCol,val);
  }
  // std::cout << "Constructor" << '\n';
}

Matrix::Matrix(int n, double val)
  : Matrix(n,n,val)
{
  // std::cout << "Constructor" << '\n';
}

// constructor to INITIALIZE with a list of values
Matrix::Matrix(int nrRow, int nrCol, std::initializer_list<double> l)
{
  #ifndef SAFETYOFF
  if (nrRow < 1 || nrCol < 1)
  {
    std::cerr << "Constructor error, dimension must be positive" << '\n';
    exit(1);
  }
  #endif
  mNRow = mRowRes = nrRow;
  mNCol = mColRes = nrCol;
  mCopy = false;
  mData = new double*[mNRow];
  int k = 0;
  int lenlist = l.end()-l.begin();

  #ifndef SAFETYOFF
  if (lenlist>mNRow*mNCol)
  {
    std::cerr<<"Too many initialization values\n";
    std::cerr<<"extra entries are ignored\n";
  }
  else if (lenlist<mNRow*mNCol)
  {
    std::cerr << "Too few initialization values" << '\n';
    std::cerr << "Remaining elements initialized as 0.0" << '\n';
  }
  #endif

  for (int i=0; i<mNRow; i++)
  {
    mData[i] =  new double[mNCol];

    for (int j=0; j<mNCol; j++)
    {
      if (k<lenlist)
      {
        mData[i][j] = l.begin()[k];
      }
      else
      {
        mData[i][j]=0.0;
      }
      k+=1;
    }
  }
  // std::cout << "Constructor" << '\n';
}

// made privat to ensure use only in member functions, where one can call properly. no fail safes against mismatching nrRow/nrCol and elements in array, dangerous as pointers are manipulated
Matrix::Matrix(int nrRow, int nrCol, double**& pt_to_cont)
{
  #ifndef SAFETYOFF
  if (nrRow < 1 || nrCol < 1)
  {
    std::cerr << "Constructor error, dimension must be positive" << '\n';
    exit(1);
  }
  #endif

  mNRow = mRowRes = nrRow;
  mNCol = mColRes = nrCol;
  mCopy = true;

  if (nrRow==1)
  {
        //mData= pt_to_cont;
    mData = new double*[1];
    mData[0]= pt_to_cont[0];
    // for (int i=0;i<nrCol; i++)
    // {
    //   std::cout << "pt[i] = " << pt_to_cont[i] << ' ';
    //   std::cout << "test+i = " << pt_to_cont[0]+i << ' ';
    //   std::cout << "mData[] = " << &mData[0][i] << ' ';
    //   std::cout << "mData+1 = " << mData[0]+i << '\n';
    // }
    //mData = &pt_to_cont[0];
  }
  else //if (nrCol == 1)
  {
    mData = new double*[mNRow];
    for (int i=0; i<mNRow; i++)
    {
      mData[i]=pt_to_cont[i];
    }
  }
  // else
  // {
  //   throw Exception("access error", "Can only be used for single rows/columns");
  //   mData = new double*[mNRow];
  //   for (int i=0)
  // }
  // std::cout << "Constructor" << '\n';
}

// privat constructor for arithmetik, no unneccessary assignment that will be overwritten later, also private
Matrix::Matrix(int nrRow, int nrCol, char empty)
{
  #ifndef SAFETYOFF
  if (nrRow < 1 || nrCol < 1)
  {
    std::cerr << "Constructor error, dimension must be positive" << '\n';
    exit(1);
  }
  if (empty!='e')
  {
    std::cerr << "Warning, this should be 'e'" << '\n';
  }
  #endif
  mNRow = mRowRes = nrRow;
  mNCol = mColRes = nrCol;
  mCopy = false;
  mData = new double*[mNRow];
  for (int i=0; i<mNRow; i++)
  {
    mData[i] =  new double[mNCol];
  }
  // std::cout << "Constructor" << '\n';
}




Matrix::Matrix(std::string filename)
{
  std::ifstream in(filename);

  if (!in) {
    std::cout << "Cannot open file.\n";
    exit(1);
  }
  in >> mNRow;
  mRowRes = mNRow;
  in >> mNCol;
  mColRes = mNCol;
  mCopy = false;

  //std::cout << "Opened file and got " << mNRow << " " << mNCol << '\n';
  mData =  new double*[mNRow];

  for (int i=0; i<mNRow; i++)
  {
    mData[i] =  new double[mNCol];
    for (int j=0; j<mNCol; j++)
    {
        in >> mData[i][j];
    }
  }
  in.close();
}


Matrix::Matrix(int nrRow, int nrCol, bool name, int resnrRow, int resnrCol)
{
  #ifndef SAFETYOFF
  // if (nrRow < 1 || nrCol < 1)
  // {
  //   std::cerr << "Constructor error, dimension must be positive" << '\n';
  //   exit(1);
  // }
  if (nrRow>resnrRow || nrCol>resnrCol)
  {
    std::cerr << "number of reserved rows/cols has to be larger than the actual number" << '\n';
  }
  #endif

  mNRow = nrRow;
  mRowRes = resnrRow;
  mNCol = nrCol;
  mColRes = resnrCol;
  mCopy = false;
  mData = new double*[mRowRes];
  for (int i=0; i<mRowRes; i++)
  {
    mData[i] =  new double[mColRes];
    std::fill_n(mData[i],mColRes,0);
  }
  // std::cout << "Constructor" << '\n';
}

static Matrix MatrixWithRes(int nrRow, int nrCol, int resnrRow, int resnrCol)
{
  return Matrix(nrRow, nrCol, 1, resnrRow, resnrCol);
}


///// destructor, special behaviour for copies
Matrix::~Matrix()
{
  // std::cout << "deb destructor\n mRow= "<< mNRow << " mCol= " << mNCol << " resRow= " << mRowRes << " rescol= " << mColRes << " and matrix " << *this << '\n';

  if ( !mCopy )
  {
    for (int i=0; i<mRowRes; i++)
    {
      delete[] mData[i];
    }
  }
  delete[] mData;
  // std::cout << "destructor" << '\n';
}

//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// indexing, get/set
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// mData, mNRow, mNCol  fields of () preceding matrix, we return reference to mData[i-1][j-1], allowing us to change it
const double& Matrix::operator()(int i, int j) const
{
  #ifndef SAFETYOFF
  if (i<1 || j<1)
  {
    std::cerr << "exception in M("<<i <<","<< j <<")" << '\n';
    throw Exception("out of range",
    "accessing matrix through () - index too small/negative");
  }
  else if (i > mNRow || j>mNCol)
  {
    std::cerr << "exception in M("<<i <<","<< j <<")" << '\n';
    throw Exception("length mismatch",
    "accessing matrix through () - index too high");
  }
  #endif

  return mData[i-1][j-1];
}

double& Matrix::operator()(int i, int j)
{
  return const_cast<double&>(static_cast<const Matrix &>(*this)(i, j));
}




// get single element
// can use double&, gives (temporary) reference address to object that is staying around
const double& Matrix::operator()(int i) const
{
  #ifndef SAFETYOFF
  if (i> mNRow*mNCol)
  {
    std::cerr << "exception in M("<< i << ")" << '\n';
    throw Exception("length mismatch",
    " accessing matrix through () - index too high");
  }
  #endif

  int col= (i-1)/mNRow;
  int row= (i-1) % mNRow;

  return mData[row][col];
}




double& Matrix::operator()(int i)
{
  return const_cast<double&>(static_cast<const Matrix &>(*this)(i));
}




// get row of matrix
// return temporary object, matrix defined in here
Matrix Matrix::operator()(int i, char c)
{
  #ifndef SAFETYOFF
  if (c!=':')
  {
    throw Exception("wrong use",
    " char must be :");
  }
  if (i>mNRow)
  {
    std::cerr << "exception in M("<<i <<","<< ":)" << '\n';
    throw Exception("length mismatch",
    " accessing matrix through () - index too high");
  }
  #endif

  double** tmp_pt;
  tmp_pt = new double* [1];

  tmp_pt[0]=&mData[i-1][0];

  Matrix M(1,mNCol,tmp_pt);

  return M;
}




// get column of matrix
Matrix Matrix::operator()(char c, int i)
{
  #ifndef SAFETYOFF
  if (c!=':')
  {
    throw Exception("wrong use",
    " char must be :");
  }
  if (i>mNRow)
  {
    std::cerr << "exception in M(:,"<<i << ")" << '\n';
    throw Exception("length mismatch",
    " accessing matrix through () - index too high");
  }
  #endif

  //Matrix M(mNRow, 1);
  double** tmp_pt;
  tmp_pt = new double* [mNRow];

  for (int j=0; j<mNRow; j++)
  {
    tmp_pt[j]=&mData[j][i-1];
  }

  Matrix M(mNRow,1,tmp_pt);

  return M;
}


Matrix Matrix::operator()(std::initializer_list<int> r1, std::initializer_list<int> r2)
{
  #ifndef SAFETYOFF
  if (r1.end()-r1.begin() != 2 || r2.end()-r2.begin() != 2)
  {
    std::cerr << "get range error: only ranges of 2 are allowed" << '\n';
  }
  #endif

  int row1=r1.begin()[0];
  int row2=r1.begin()[1];
  int col1=r2.begin()[0];
  int col2=r2.begin()[1];
  int nrrow=row2-row1+1;
  int nrcol=col2-col1+1;
  // std::cout << nrrow << " " << row1 << " "<< row2 << " "<< nrcol << " " << '\n';

  #ifndef SAFETYOFF
  if (nrrow<0 || row1<=0 || row2>mNRow || nrcol<0 || col1<=0 || col2>mNCol )
  {
    std::cerr << "Access error, invalid range (increasing lists of two elements in range)" << '\n';
  }
  #endif

  double** tmp_pt;
  tmp_pt = new double* [nrrow];

  for (int j=0; j<nrrow; j++)
  {
    tmp_pt[j]=&mData[row1-1+j][col1-1];
  }

  Matrix M(nrrow,nrcol,tmp_pt);

  return M;
}


//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// storage structure altering members
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// transpose
void Matrix::T()
{
  double** target;
  target = new double*[mColRes];
  for (int i = 0; i < mColRes; i++)
  {
    target[i] = new double[mRowRes];
    for (int j = 0; j < mRowRes; j++)
    {
      target[i][j]= mData[j][i];
    }
  }
  for (int i = 0; i< mRowRes; i++)
  {
    delete[] mData[i];
  }
  delete[] mData;

  mData=target;
  std::swap(mNRow,mNCol);
  std::swap(mRowRes,mColRes);
  // return *this;
}


// append
void Matrix::append(const Matrix& m)
{
  // if appended vector has more rows, pad original with zeros from below
  // std::cout << " deb appending mRow= "<< m.mNRow+mNRow << " mCol= " << mNCol+m.mNCol << " resRow= " << mRowRes << " rescol= " << mColRes << '\n';
  int tarrow = std::max(m.mNRow,mNRow);
  if ( (mRowRes >= tarrow) && (mColRes >= mNCol+m.mNCol) )
  {
    // std::cout << "enough reserved" << '\n';
    for (int i=0; i< m.mNRow;  i++)
    {
      for (int j=0; j< m.mNCol; j++)
      {
        mData[i][j+mNCol]= m.mData[i][j];
      }
    }
    mNCol = mNCol+m.mNCol;
    mNRow = tarrow;
  }
  else
  {
    // std::cout << "not enough reserved space" << '\n';
    int minrow = std::min(m.mNRow,mNRow);
    int tarcol = mNCol+m.mNCol;
    double** tmpData;
    tmpData = new double* [tarrow];

    for (int i=0; i<minrow; i++)
    {
      tmpData[i] = new double[tarcol];
      std::copy(mData[i], mData[i]+mNCol, tmpData[i]);
      std::copy(m.mData[i], m.mData[i]+m.mNCol, tmpData[i]+mNCol);
      delete[] mData[i];
    }

    if ( mNRow > minrow)
    {
      Matrix Z(1 ,m.mNCol);
      for (int i=minrow; i<mNRow; i++)
      {
        tmpData[i] = new double[tarcol];
        std::copy(mData[i], mData[i]+mNCol, tmpData[i]);
        std::copy(Z.mData[0], Z.mData[0]+m.mNCol, tmpData[i]+mNCol);
        delete[] mData[i];
      }
      std::cerr<<"matrix += - matrices have different nRow\n";
      std::cerr<<"appended matrix padded with 0.0\n";
    }
    else if ( m.mNRow > minrow)
    {
      Matrix Z(1 ,mNCol);
      for (int i = mNRow; i< m.mNRow; i++)
      {
        tmpData[i] = new double[mNCol+m.mNCol];
        std::copy(Z.mData[0], Z.mData[0]+mNCol, tmpData[i]);
        std::copy(m.mData[i], m.mData[i]+m.mNCol, tmpData[i]+mNCol);
      }
      mNRow=m.mNRow;
      // std::cerr<<"matrix += - matrices have different nRow\n";
      // std::cerr<<"matrix padded with 0.0\n";
    }

    delete[] mData;
    mData = tmpData;
    mNCol = tarcol;
  }
  // return *this;
}




// swaprow
void Matrix::swaprow(int i, int j)
{
  #ifndef SAFETYOFF
  if ( i< 1 || j < 1 || i > mNRow || j>mNRow)
  {
    std::cerr << "Error: row swap index out of bounds" << '\n';
    exit(1);
  }
  #endif

  std::swap(mData[i-1],mData[j-1]);
  // return *this;
}



// givens rotation
double Matrix::givens(int i, int j)
{
  //std::cout << "i="<< i << "j=" << j << '\n';
  double costheta = 1;
  if (mData[j-1][i-1]!=0)
  {
    costheta = mData[i-1][i-1] / sqrt(pow(mData[i-1][i-1],2.0) + pow(mData[j-1][i-1],2.0));
  }
  //std::cout << "costheta" << costheta << '\n';
  givens(i,j,costheta);
  return costheta;
}

void Matrix::givens(int i, int j, double costheta)
{
    double c = costheta;
    double s = sqrt(1-pow(costheta,2.0));
    i-=1; // switch to zero indexing
    j-=1;
    double tempi;
    double tempj;
    for (int k = 0; k < mNCol; k++)
    {
      tempi = c * mData[i][k] + s * mData[j][k];
      tempj = -s * mData[i][k] + c * mData[j][k];
      mData[i][k] = tempi;
      mData[j][k] = tempj;
    }
}


//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// assignment
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// assignment operator
Matrix& Matrix::operator=(const Matrix& m)
{
  // std::cout << " deb assignment mRow= "<< mNRow << " mCol= " << mNCol << " resRow= " << mRowRes << " rescol= " << mColRes << '\n';

  if (this == &m)      // Same object?
  {
    return *this;
  }

  #ifndef SAFETYOFF
  if (m.mNRow != mNRow || m.mNCol != mNCol)
  {
    std::cerr << "Exception assignment dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
	  "matrix assignment operator - lhs and rhs have different dimension");
  }
  #endif

  for (int i=0; i<m.mNRow; i++)
   {
     for (int j=0; j<m.mNCol; j++)
     {
       mData[i][j] = m.mData[i][j];
     }
   }

   mRowRes = mNRow;
   mColRes = mNCol;

  return *this;
}


// define appending of columns, comparably expensive
Matrix& Matrix::operator+=(const Matrix& m)
{
  // std::cout << "deb sum" << '\n';
  #ifndef SAFETYOFF
  if (mNRow != m.mNRow || mNCol != m.mNCol)
  {
    std::cerr << "Exception: addition dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
    "addition is not defined for matrices of different size");
  }
  #endif

  for (int i=0; i< mNRow; i++)
  {
    for (int j=0; j<mNCol; j++)
    {
      mData[i][j] = mData[i][j] + m.mData[i][j];
    }
  }
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& m)
{
  // std::cout << "deb subtraction" << '\n';
  #ifndef SAFETYOFF
  if (mNRow != m.mNRow || mNCol != m.mNCol)
  {
    std::cerr << "Exception: subtraction dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
    "subtraction is not defined for matrices of different size");
  }
  #endif

  for (int i=0; i< mNRow; i++)
  {
    for (int j=0; j<mNCol; j++)
    {
      mData[i][j] = mData[i][j] - m.mData[i][j];
    }
  }
  return *this;
}

//Matrix& Matrix::operator*=(const Matrix& m)
Matrix& Matrix::operator*=(Matrix m)
{
  // std::cout << "deb multiplication" << '\n';
  #ifndef SAFETYOFF
  if ( mNCol != m.mNRow )
  {
    std::cerr << "Exception: multiplication dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
    "multiplication is not defined for matrices of mismatching row/columns");
  }
  #endif
  // A m,n; n,p
  //int mr = mNRow; int nrc = mNCol;
  int pc = m.mNCol;
  double** tmpData;
  m.T();

  //std::cout << "test2" << m << '\n';
  // std::cout << "deb mRow= "<< mNRow << " mCol= " << mNCol << " resRow= " << mRowRes << " rescol= " << mColRes << '\n';
  // std::cout << "deb LHS= "<< *this << "RHS= " << m << '\n';
  tmpData = new double* [mNRow];
  for (int i=0; i< mNRow; i++)
  {
    // std::cout << "deb first loop" << '\n';
    tmpData[i]= new double[pc];
    for (int j=0; j<pc; j++)
    {
      // std::cout << "deb test " << i << j << "0" << '\n';
      // tmpData[i][j] = mData[i][0] * m.mData[0][j];
      tmpData[i][j] = mData[i][0] * m.mData[j][0];
      for (int k = 1; k<mNCol ;k++)
      {
          //std::cout << "test " << i << j << k << '\n';
          // tmpData[i][j] = tmpData[i][j] + mData[i][k] * m.mData[k][j];
          tmpData[i][j] = tmpData[i][j] + mData[i][k] * m.mData[j][k];
      }

    }
    delete[] mData[i];
  }
  // std::cout << "deb delete mData" << '\n';
  delete[] mData;
  mData=tmpData;
  mNCol=pc;
  // std::cout << "deb return multiplication" << '\n';

  return *this;
}

Matrix& Matrix::operator*=(const double& a)
{
  for (int i=0; i< mNRow; i++)
  {
    for (int j=0; j<mNCol; j++)
    {
      mData[i][j] = a * mData[i][j];
    }
  }
  return *this;
}




Matrix& Matrix::operator/=(const double& a)
{
  for (int i=0; i< mNRow; i++)
  {
    for (int j=0; j<mNCol; j++)
    {
      mData[i][j] = mData[i][j] / a;
    }
  }
  return *this;
}


//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// member functions
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int size(const Matrix& m, int dir)

{
int out;
  switch (dir)
  {
    case 1: out = m.mNRow; break;
    case 2: out = m.mNCol; break;
    default:
    {
      out = m.mNRow;
      std::cerr << "size error, undefined dimension given, default gives nr of rows" << '\n';

    }
  }
  return out;
}

int reservedsize(const Matrix& m, int dir)

{
int out;
  switch (dir)
  {
    case 1: out = m.mRowRes; break;
    case 2: out = m.mColRes; break;
    default:
    {
      out = m.mNRow;
      std::cerr << "size error, undefined dimension given, default gives nr of rows" << '\n';

    }
  }
  return out;
}


double Matrix::norm(int p) const
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i=0; i<mNRow; i++)
  {
    for (int j=0; j< mNCol; j++)
    {
      temp = fabs(mData[i][j]);
      norm_val += pow(temp, p);
    }

  }

  return pow(norm_val, 1.0/((double) (p)));
}




//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// arithmetik operators free functions
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// plus
Matrix operator+(const Matrix& m1, const Matrix& m2)
{
  #ifndef SAFETYOFF
  if (size(m1,1) != size(m2,1) || size(m1,2) != size(m2,2))
  {
    std::cerr << "Exception: addition dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
    "addition is not defined for matrices of different size");
  }
  #endif

  Matrix m1c (m1);
  m1c += m2;
  return m1c;
}




 // subtraction of matrices
Matrix operator-(const Matrix& m1, const Matrix& m2)
{
  #ifndef SAFETYOFF
  if (size(m1,1) != size(m2,1) || size(m1,2) != size(m2,2))
  {
    std::cerr << "Exception: subtraction dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
    "subtraction is not defined for matrices of different size");
  }
  #endif
  Matrix m1c (m1);
  m1c -= m2;
  return m1c;
}




// matrix multiplication
Matrix operator*(const Matrix& m1, const Matrix& m2)
{
  #ifndef SAFETYOFF
  if ( size(m1,2) != size(m2,1) )
  {
    std::cerr << "Exception: multiplication dimension mismatch" << '\n';
    throw Exception("dimension mismatch",
    "multiplication is not defined for matrices of mismatching row/columns");
  }
  #endif

  Matrix m1c (m1);
  m1c *= m2;
  return m1c;
}




Matrix operator*(const double& a, const Matrix& m)
{
  Matrix mc (m);
  mc *= a;
  return mc;
}

Matrix operator*(const Matrix& m, const double& a)
{
  Matrix mc (m);
  mc *= a;
  return mc;
}


Matrix operator/(const Matrix& A, Matrix b)
{
  LUd lud(A);
  lud.solve(b);
  return b;
}

Matrix operator/(const Matrix& m, const double& a)
{
  Matrix mc (m);
  mc /= a;
  return mc;
}

Matrix operator/(const double& a, const Matrix& m)
{
  Matrix mc (m);
  mc /= a;
  return mc;
}



///// output/ print a matrix
std::ostream& operator<<(std::ostream& output, const Matrix& m)
{
  int nrrow=size(m,1);
  int nrcol=size(m,2);
  output << "\n(";
  for (int i=1; i<=nrrow; i++)
    {
      for (int j=1; j<=nrcol; j++)
      {
        output << std::fixed << std::setprecision(3) << m(i,j);
        if (j != nrcol)
  	       output  << ", ";
      }
      if (i != nrrow)
         output  << ";\n ";
      else
          output  << ")\n";
    }
  return output;  // for multiple << operators.
}



Matrix transpose(const Matrix & A)
{
  Matrix M(A);
  M.T();
  return M;
}

double norm(Matrix& v, int p)
{
  return v.norm(p);
}



void Matrix::gaussj(Matrix& b)
{
  #ifndef SAFETYOFF
  if ( mNRow != mNCol || mNRow!=b.mNRow)
  {
    std::cerr << "Error, Matrix is not square or A,b have different row number" << '\n';
    exit(1);
  }
  #endif

  int indxr [mNCol];
  int indxc [mNCol];
  int ipiv [mNCol];
  std::fill_n(ipiv,mNCol,0);
  double cmp, dum, pivinv;
  int icol, irow;

  for (int i = 0; i<mNCol; i++)
  {
    cmp=0.0;
    for (int j = 0; j<mNRow; j++)
    {
      if (ipiv[j] != 1)
      {
        for (int k=0; k<mNCol; k++)
        {
          if (ipiv[k]==0)
          {
            if (fabs(mData[i][j])>=cmp)
            {
              cmp=mData[j][k];
              irow=j;
              icol=k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    //std::cout << "Pivot = " << pivrow[j] << '\n';
    if (irow != icol){
      std::swap(mData[irow], mData[icol]);
      std::swap(b.mData[irow], b.mData[icol]);
      // std::cout << "Swapped" << '\n';
      //std::cout << "Swapped b = " << b << '\n';
    }
    indxr[i]=irow;
    indxc[i]=icol;

    if ( mData[icol][icol] == 0.0)
    {
      std::cerr << "Error, Matrix is singular" << '\n';
      exit(1);
    }

    pivinv = 1.0/ mData[icol][icol];
    mData[icol][icol] = 1.0; //??
    for (int l = 0; l < mNCol; l++)
    {
      mData[icol][l] *= pivinv;
    }
    for (int l = 0; l < b.mNCol; l++)
    {
      b.mData[icol][l] *= pivinv;
    }
    //std::cout << "Normalized = " << m << '\n';

    for (int ll=0; ll<mNRow; ll++ )
    {
      if (ll != icol)
      {
        dum=mData[ll][icol];
        mData[ll][icol]=0.0;
        for (int l = 0; l < mNCol; l++)
        {
          mData[ll][l] -= mData[icol][l]*dum;
        }
        for (int l = 0; l < b.mNCol; l++)
        {
          b.mData[ll][l] -= b.mData[icol][l]*dum;
        }
        //std::cout << "Remove = " << m << '\n';
      }
    }
  }
  for (int l =mNCol-1; l>=0 ; l--)
  {
    if (indxr[l] != indxc[l]){
      for (int k=0 ; k< mNRow; k++)
      std::swap(mData[k][indxr[l]], mData[k][indxc[l]]);
    }
  }
}



Matrix inv(Matrix& m)
{
  Matrix b(size(m,1),1);
  Matrix minv=m;
  minv.gaussj(b);
  return minv;
}




Matrix eye(int rows, int cols)
{
  Matrix I(rows,cols);
  int dim=std::min(rows,cols);
  for (int i = 1; i<= dim; i++)
  {
    I(i,i) = 1.0;
  }
  return I;
}




Matrix eye(int dim)
{
  //Matrix I(3,3); //= eye(dim,dim);
  //return I;
  return eye(dim,dim);
}




// if input is vector creates appropriate diagonal matrix,
// if input is matrix, return row vector with elements from diagonal
Matrix diag(const Matrix& m, int offs)
{
  if (size(m,1)==1 || size(m,2)==1)
  {
    int dim = std::max(size(m,2),size(m,1));
    Matrix D(dim+abs(offs));
    int coloff = std::max(0,offs);
    int rowoff = std::max(0,-offs);
    for (int i=1; i<= dim; i++)
    {
      D(i+rowoff,i+coloff) = m(i);
      //D(1,1) = 1;//m(1);
    }
    return D;
  }
  else
  {
    #ifndef SAFETYOFF
    if ( !( (offs < size(m,2)) && (offs > -size(m,1))))
    {
      //std::cout << "logic" << (offs < size(m,2))  << (offs > -size(m,1)) << ( (offs < size(m,2)) && (offs > -size(m,1))) <<'\n';
      std::cerr << "Error: diag, offset out of bounds" << '\n';
      exit(1);
    }
    #endif

    int dim=std::min(size(m,1),size(m,2))- abs(offs);
    Matrix D(1,dim);
    int coloff = std::max(0,offs);
    int rowoff = std::max(0,-offs);
    for (int i = 1; i<=dim; i++)
    {
      D(i)=m(i+rowoff,i+coloff);
    }
    return D;
  }
}


Matrix GMRES(const Matrix & A, const Matrix& b,const  Matrix& x0, double Tol, std::ostream &file)
{
  if (size(A,2) != size(A,1))
  {
    std::cerr << "Error, GMRES requires square matrix" << '\n';
    exit(1);
  }
  // std::cout << "deb inside" << '\n';

  int n = size(A,2);
  int l;
  double res[n];
  double hl1l, s;
  double costheta[n];

  // std::cout << "deb instatiated" << '\n';

  Matrix r = b - A * x0;
  // std::cout << "deb computed r" << '\n';
  res[0]  = norm(r);
  if (file)
  {
    file << std::setprecision(15) << res[0] << ", ";
  }
  Matrix V = r/norm(r);
  // for (int i = 1; i<= size(V,1); i++ )
  // {
  //   file << std::setprecision(15) << V(i) << ", ";
  // }
  // file << "\n";
  Matrix w(n,1);
  Matrix H(2,1);

  // std::cout << "deb start loop" << '\n';
  for  ( l=1;  l <= n ; l++)
  {
    // std::cout << "l= " << l << " n= " << n << '\n';
    w = A*V(':',l);
    {
      Matrix hl(l+1,1);
      for (int j = 1; j<=l; j++)
      {
        hl(j)=scaprod(V(':',j) , w);
        w -= (hl(j) * V(':',j));
      }
      hl(l+1) = hl1l = norm(w);

      for (int j = 1; j<l; j++)
      {
        hl.givens(j,j+1, costheta[j]);
      }

      // std::cout << "deb append H" << '\n';

      if (l==1)
      {
        H=hl;
      }
      else
      {
        H.append(hl);
      }
      costheta[l] = H.givens(l,l+1);
      s = sqrt(1-pow(costheta[l],2.0));
      res[l] = s * res[l-1];

      if (file)
      {
        file << std::setprecision(15) << res[l] << ", ";
      }
      // std::cout << "residue" << res[l] << '\n';

      // std::cout << std::scientific << "hl1l = " << fabs(hl1l) << " ref = " << 1e-17 << '\n';

      if ((res[l] < Tol) || fabs(hl1l) <= 1e-14 || l==n  )
      {
        break;
      }
      else
      {
        // for (int i = 1; i<= size(w,1); i++ )
        // {
        //   file << std::setprecision(15) << w(i)/hl(l+1) << ", ";
        // }
        // file << "\n";
        // std::cout << "deb append V" << '\n';
        V.append(w/hl(l+1));
      }
    }
  }
  if (file)
  {
    file << "\n";// << l << ", ";
  }
  // std::cout << "V=" << V << "H=" << H << '\n';

  Matrix QT=eye(l+1,1);
  for (int k = 1;  k <= l ; k++)
  {
    QT.givens(k,k+1,costheta[k]);
  }

/// Back substitution to solve
  Matrix y = res[0]*QT({1,l},{1,1});

  for (int i = l; i > 0; i--)
  {
    double sum = y(i);
    for (int j= i+1; j<= l; j++)
    {
      sum-= H(i,j) * y(j);
    }
    y(i) = sum/H(i,i);

  }

  // std::cout << "lu decom \n Hhat" << H({1,l},{1,l}) << " res* QT" << res[0]*QT({1,l},{1,1}) << "Solution " <<  y << "solution test" << H({1,l},{1,l})*y << '\n';

  Matrix x = (x0 + V({1,n},{1,l}) * y) ;

  // std::cout << "Solution" << x << '\n';

  return x;
}



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/// GMRES 2
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


Matrix GMRES2(const Matrix & A, const Matrix& b,const  Matrix& x0, double Tol, std::ostream &file)
{
  if (size(A,2) != size(A,1))
  {
    std::cerr << "Error, GMRES requires square matrix" << '\n';
    exit(1);
  }
  int n = size(A,2);
  int l;
  double res[n];
  double hl1l, s;
  double costheta[n];


  Matrix r= b - A * x0;
  res[0]  = norm(r);
  if (file)
  {
    file << std::setprecision(15) << res[0] << ", ";
  }
  Matrix V = MatrixWithRes(n,1,n,n);
  V(':',1) = r/norm(r);
  // for (int i = 1; i<= size(V,1); i++ )
  // {
  //   file << std::setprecision(15) << V(i) << ", ";
  // }
  // file << "\n";
  Matrix w(n,1);
  // Matrix H(2,1);
  // std::cout << "build H(1,0)" << '\n';
  Matrix H = MatrixWithRes(1,0,n+1,n);


  for  ( l=1;  l <= n ; l++)
  {
    // std::cout << "l= " << l << " n= " << n << '\n';
    w = A*V(':',l);
    //H.resResize(1,1);
    {
      Matrix hl(l+1,1);
      for (int j = 1; j<=l; j++)
      {
        hl(j)=scaprod(V(':',j) , w);
        w -= (hl(j) * V(':',j));
      }
      hl(l+1) = hl1l = norm(w);

      for (int j = 1; j<l; j++)
      {
        hl.givens(j,j+1, costheta[j]);
      }
      // std::cout << "old H" << H << '\n';

      // std::cout << "deb resRow= " << reservedsize(H,1) << " rows= " << size(H,1) << " resCol= " << reservedsize(H,2) << " cols= " << size(H,2) <<  '\n';
      H.append(hl);
      // std::cout << "new H" << H << '\n';
      costheta[l] = H.givens(l,l+1);
      s = sqrt(1-pow(costheta[l],2.0));
      res[l] = s * res[l-1];
      if (file)
      {
        file << std::setprecision(15) << res[l] << ", ";
      }
      // std::cout << "residue" << res[l] << '\n';

      // std::cout << std::scientific << "hl1l = " << fabs(hl1l) << " ref = " << 1e-17 << '\n';

      if ((res[l] < Tol) || fabs(hl1l) <= 1e-14 || l==n  )
      {
        break;
      }
      else
      {
        // for (int i = 1; i<= size(w,1); i++ )
        // {
        //   file << std::setprecision(15) << w(i)/hl(l+1) << ", ";
        // }
        // file << "\n";
        // std::cout << "deb append to V" << '\n';
        // std::cout << "deb resRow= " << reservedsize(V,1) << " rows= " << size(V,1) << " resCol= " << reservedsize(V,2) << " cols= " << size(V,2) <<  '\n';
        // V.resResize(0,1);
        V.append(w/hl(l+1));
        // std::cout << "deb V" << V << '\n';
      }
    }
  }
  if (file)
  {
    file << "\n";// << l << ", ";
  }
  // std::cout << "V=" << V << "H=" << H << '\n';

  Matrix QT=eye(l+1,1);
  for (int k = 1;  k <= l ; k++)
  {
    QT.givens(k,k+1,costheta[k]);
  }

/// Back substitution to solve
  Matrix y = res[0]*QT({1,l},{1,1});

  for (int i = l; i > 0; i--)
  {
    double sum = y(i);
    for (int j= i+1; j<= l; j++)
    {
      sum-= H(i,j) * y(j);
    }
    y(i) = sum/H(i,i);

  }

  // std::cout << "lu decom \n Hhat" << H({1,l},{1,l}) << " res* QT" << res[0]*QT({1,l},{1,1}) << "Solution " <<  y << "solution test" << H({1,l},{1,l})*y << '\n';

  Matrix x = (x0 + V({1,n},{1,l}) * y) ;

  // std::cout << "Solution" << x << '\n';

  return x;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//// random extra functions
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double mat2d(const Matrix& M)
{
  return M(1,1);
}


double scaprod(const Matrix& m1, const Matrix& m2)
{
  #ifndef SAFETYOFF
  if (size(m1,1) != size(m2,1) || size(m1,2) != 1 || size(m2,2)!=1)
  {
    std::cerr << "Error: scalar product needs two column vector of equal size" << '\n';
  }
  #endif

  double tmp =0;
  int n= size(m1,1);

  for (int i = 1; i <=n; i++)
  {
    tmp += m1(i)*m2(i);
  }
  return tmp;
}


bool iscopy(const Matrix& m)
{
  return m.mCopy;
}
