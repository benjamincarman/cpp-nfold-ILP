/*******************************************************************************
 *
 * File: main.cc
 *
 * Author: Benjamin Carman
 *
 * Date: June 4, 2019
 *
 * Description: A testing harness for NFold class.
 *
 ******************************************************************************/

#include "gurobi_c++.h"
#include "nfold.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace std;

int main(int argc, char *argv[])
{
  //Start Gurobi environment for each solve iteration
  GRBEnv env = GRBEnv();
  env.set("LogFile", "info.log");
  env.start();

  //Output files for objective value data and time data
  fstream fObjVals;
  fObjVals.open("objVals.txt", fstream::app);
  fstream fTimes;
  fTimes.open("times.txt", fstream::app);

  //Run each input file from command line arguments
  for (int i = 1; i < argc; i++)
  {
    string filename = argv[i];
    ifstream ins;
    ins.open(filename);

    if (!ins.fail())
    {
      //Measure time at start of solve iteration
      auto start = chrono::high_resolution_clock::now();

      //Create nfold instance
      NFold nfInstance(&env);
      nfInstance.inputState(ins);
      ins.close();
      cout << "-------------------------------------------------------------------" << endl;
      cout << "Beginning New Solve from \"" << filename << "\":" << endl;
      cout << "-------------------------------------------------------------------" << endl << endl;

      //Use gc of 5
      nfInstance.setGraverComplexity(5);
      bool solved = nfInstance.solve();
      if (!solved)
      {
        cout << "Error: Instance not solved. Continuing to next input file." << endl;
        continue;
      }
      nfInstance.outputState(cout);

      //Measure, output, and record time taken
      auto finish = chrono::high_resolution_clock::now();
      chrono::duration<double> elapsed = finish - start;
      cout << "Elapsed time: " << elapsed.count() << " s" << endl << endl;
      fTimes << elapsed.count() << endl;

      //Record final objective value
      int objValue = nfInstance.getOptimizedObjectiveValue();
      fObjVals << objValue << endl;

    }
    else
    {
      cout << "Failed to open file: " << filename << ". Skipping solve." << endl;
    }
  }

  //Close output files
  fObjVals.close();
  fTimes.close();

  return EXIT_SUCCESS;
}
