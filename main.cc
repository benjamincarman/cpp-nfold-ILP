/*******************************************************************************
 *
 * File: main.cc
 *
 * Author: Benjamin Carman
 *
 * Date: June 4, 2019
 *
 * Description:
 *
 ******************************************************************************/

 #include "gurobi_c++.h"
 #include "nfold.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  GRBEnv env = GRBEnv();
  env.set("LogFile", "info.log");
  env.start();

  for (int i = 1; i < argc; i++)
  {
    NFold nfInstance(&env);
    string filename = argv[i];

    cout << "Beginning New Solve from \"" << filename << "\":" << endl;
    ifstream ins;
    ins.open(filename);
    nfInstance.instantiate(ins);
    ins.close();

    nfInstance.setGraverComplexity(5);
    nfInstance.solve();
    nfInstance.outputState(cout);
  }
  return EXIT_SUCCESS;
}
