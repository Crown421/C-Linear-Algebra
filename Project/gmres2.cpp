#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <string>
//#include <fstream>
#include <cstdio>
#include <ctime>
#include <chrono>
#include "Exception.hpp"
#include "matrix.hpp"
#include "LUd.hpp"



int main(int argc, char* argv[])
{
  std::string filenamearg;
  if (argc < 2)
  {
    std::cerr << "No file name for matrix given\n Please provide name now" << '\n';
    std::cin >> filenamearg;
  }
  else
  {
    filenamearg = argv[1];
  }

  ///////////////////////// TODO change here
  std::string outfilename = "../../../../../Dropbox/Software/MatLab_Files/Oxford/C++ Supplements/results/" + filenamearg + "_gmres2.txt";

  std::string flag;
  if (argc < 3)
  {
    std::cerr << "No flag given, run base version" << '\n';
  }
  else
  {
    flag = argv[2];
    outfilename = "../../../../../Dropbox/Software/MatLab_Files/Oxford/C++ Supplements/results/" + filenamearg + flag +"_gmres2.txt";
  }

  std::string filename = "../Data/" + filenamearg + ".txt";

  std::cout << "filename" << filename << '\n';

  Matrix M(filename);

  std::cout << "Matrix is read, outfile is" << outfilename << '\n';

  Matrix b(size(M,1),1,1.0);

  Matrix x0(size(M,1),1);

  std::ofstream outfile(outfilename);

  if (!outfile)
  {
    std::cerr << "Did not open file" << '\n';
    exit(1);
  }

  auto wcts = std::chrono::system_clock::now();
  std::clock_t startcputime = std::clock();

  //////////////////////////////////////////////////////////////

  Matrix x = GMRES2(M,b,x0,1e-10,outfile);

  //////////////////////////////////////////////////////////////
  double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);

  std::cout << "Finished in " << cpu_duration << " seconds [CPU Clock] " << wctduration.count() << " seconds [Wall Clock]" <<  std::endl;

  // outfile << std::setprecision(13) << cpu_duration << "\n";
  outfile << std::setprecision(13) << wctduration.count() << "\n";
  for (int i = 1; i<= size(x,1); i++ )
  {
    outfile << std::setprecision(15) << x(i) << ", ";
  }

  Matrix tempb=b-M*x;
  double res = norm(tempb);
  outfile << std::setprecision(15) << "\n" << res;

  outfile.close();
  std::cout << "res = " << res << '\n';


  return(0);
}
