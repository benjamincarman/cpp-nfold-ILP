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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "nfold.h"

using namespace std;

int main()
{
  cout << "Hello n-fold ILP!" << endl;

  NFold testInput;

  ifstream ins;
  ins.open("tests/QCmax_m_30_lengths_2_3_7_13_weigths_6_5_2_1_smallest_100_largest_500_slack_r_0.90_formatted.txt");
  testInput.instantiate(ins);
  ins.close();
  testInput.setGraverComplexity(5);
  testInput.outputState(cout);

  vector<int> step = testInput.findGoodStep(1);

  cout << "Good Step: " << endl;
  for (size_t i = 0; i < step.size(); i++)
  {
    cout << step[i] << ' ';
  }
  cout << endl;
  return EXIT_SUCCESS;
}
