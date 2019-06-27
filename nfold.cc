/*******************************************************************************
 *
 * File: nfold.cc
 *
 * Author: Benjamin Carman
 *
 * Date: June 4, 2019
 *
 * Description:
 *
 ******************************************************************************/

#include "nfold.h"
#include "gurobi_c++.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

NFold::NFold()
{

}

NFold::NFold(unsigned int n1, std::vector<int> objective1, std::vector<int> lowerBound1,
             std::vector<int> upperBound1, std::vector<std::vector<int> > topMatrix1,
             std::vector<std::vector<int> > diagMatrix1, std::vector<int> b1,
             std::optional<std::vector<int> > initial)
{
  n = n1;
  r = topMatrix.size();
  s = diagMatrix.size();
  t = topMatrix[0].size(); //TODO: error check
  objective = objective1;
  lowerBound = lowerBound1;
  upperBound = upperBound1;
  topMatrix = topMatrix1;
  diagMatrix = diagMatrix1;

  buildConstraintMatrix();

  b = b1;
  currentSolution.resize(n * t);

  if (initial == std::nullopt)
  {
    //find inital feasible solution
  }
  else
  {
    currentSolution = initial.value();
  }
}

NFold::NFold(unsigned int n1, unsigned int r1, unsigned int s1, unsigned int t1,
             std::vector<int> objective1, std::vector<int> lowerBound1,
             std::vector<int> upperBound1, std::vector<std::vector<double> > constraintMatrix1,
             std::vector<int> b1, std::optional<std::vector<int> > initial)
{
  n = n1;
  r = r1;
  s = s1;
  t = t1;
  r = topMatrix.size();
  s = diagMatrix.size();
  t = topMatrix[0].size(); //TODO: error check
  objective = objective1;
  lowerBound = lowerBound1;
  upperBound = upperBound1;
  constraintMatrix = constraintMatrix1;

  //TODO Maybe set top and diagonal matrices here. Maybe.

  b = b1;
  currentSolution.resize(n * t);

  if (initial == std::nullopt)
  {
    //find inital feasible solution
  }
  else
  {
    currentSolution = initial.value();
  }
}

void NFold::instantiate(istream &ins)
{
  //Input matrix dimension values
  ins >> n;
  ins >> r;
  ins >> s;
  ins >> t;

  //Size all vectors/matrices accordingly
  objective.resize(n * t);
  lowerBound.resize(n * t);
  upperBound.resize(n * t);
  topMatrix.resize(r);
  for (size_t i = 0; i < r; i++)
  {
    topMatrix[i].resize(t);
  }
  diagMatrix.resize(s);
  for (size_t i = 0; i < s; i++)
  {
    diagMatrix[i].resize(t);
  }

  b.resize(r + (n * s));
  currentSolution.resize(n * t);

  //Input vector/matrix values
  for (size_t i = 0; i < n * t; i++)
  {
    ins >> objective[i];
  }
  for (size_t i = 0; i < n * t; i++)
  {
    ins >> lowerBound[i];
  }
  for (size_t i = 0; i < n * t; i++)
  {
    ins >> upperBound[i];
  }
  for (size_t i = 0; i < r; i++)
  {
    for (size_t j = 0; j < t; j++)
    {
      ins >> topMatrix[i][j];
    }
  }
  for (size_t i = 0; i < s; i++)
  {
    for (size_t j = 0; j < t; j++)
    {
      ins >> diagMatrix[i][j];
    }
  }

  buildConstraintMatrix();

  for (size_t i = 0; i < r + (n * s); i++)
  {
    ins >> b[i];
  }
  for (size_t i = 0; i < n * t; i++)
  {
    ins >> currentSolution[i];
  }
  //ADD SOME ERROR HANDLING
}

void NFold::buildConstraintMatrix()
{
  constraintMatrix.resize(r + n * s);

  for (size_t i = 0; i < r; i++)
  {
    constraintMatrix[i].resize(n * t);
    for (size_t j = 0; j < n; j++)
    {
      for (size_t k = 0; k < t; k++)
      {
        constraintMatrix[i][j * t + k] = 1.0 * topMatrix[i][k];
      }
    }
  }

  for (size_t i = r; i < r + n * s; i++)
  {
    constraintMatrix[i].resize(n * t);
    for (size_t j = 0; j < n; j++)
    {
      if ((i - r) / s == j)
      {
        for (size_t k = 0; k < t; k++)
        {
          constraintMatrix[i][j * t + k] = 1.0 * diagMatrix[(i - r) % s][k];
        }
      }
      else
      {
        for (size_t k = 0; k < t; k++)
        {
          constraintMatrix[i][j * t + k] = 0.0;
        }
      }
    }
  }
}

void NFold::outputState(std::ostream &outs)
{
  outs << "n: " << n << endl
       << "r: " << r << endl
       << "s: " << s << endl
       << "t: " << t << endl;

  outs << endl << "Objective:" << endl << "<";
  for (size_t i = 0; i < n * t - 1; i++)
  {
    outs << objective[i] << ", ";
  }
  outs << objective[n * t - 1] << ">" << endl;

  outs << endl << "Lower Bound:" << endl << "<";
  for (size_t i = 0; i < n * t - 1; i++)
  {
    outs << lowerBound[i] << ", ";
  }
  outs << lowerBound[n * t - 1] << ">" << endl;

  outs << endl << "Upper Bound:" << endl << "<";
  for (size_t i = 0; i < n * t - 1; i++)
  {
    outs << upperBound[i] << ", ";
  }
  outs << upperBound[n * t - 1] << ">" << endl;

  outs << endl << "A Matrix:" << endl << "[ ";
  for (size_t i = 0; i < r; i++)
  {
    if (i != 0) outs << "  ";
    outs << "[";
    for (size_t j = 0; j < t - 1; j++)
    {
      outs << topMatrix[i][j] << ", ";
    }
    outs << topMatrix[i][t - 1];
    if (i != r - 1) {outs << "]," << endl;}
    else {outs << "] ]" << endl;}
  }

  outs << endl << "D Matrix:" << endl << "[ ";
  for (size_t i = 0; i < s; i++)
  {
    if (i != 0) outs << "  ";
    outs << "[";
    for (size_t j = 0; j < t - 1; j++)
    {
      outs << diagMatrix[i][j] << ", ";
    }
    outs << diagMatrix[i][t - 1];
    if (i != s - 1) {outs << "]," << endl;}
    else {outs << "] ]" << endl;}
  }

  outs << endl << "An Matrix:" << endl << "[ ";
  for (size_t i = 0; i < r + n * s; i++)
  {
    if (i != 0) outs << "  ";
    outs << "[";
    for (size_t j = 0; j < n * t - 1; j++)
    {
      outs << constraintMatrix[i][j] << ", ";
    }
    outs << constraintMatrix[i][n * t - 1];
    if (i != r + n * s - 1) {outs << "]," << endl;}
    else {outs << "] ]" << endl;}
  }

  outs << endl << "RHS (b): " << endl << "<";
  for (size_t i = 0; i < r + n * s - 1; i++)
  {
    outs << b[i] << ", ";
  }
  outs << b[r + n * s - 1] << ">" << endl;

  outs << endl << "Solution: " << endl << "<";
  for (size_t i = 0; i < n * t - 1; i++)
  {
    outs << currentSolution[i] << ", ";
  }
  outs << currentSolution[n * t - 1] << ">" << endl;
}

void NFold::setGraverComplexity(unsigned int gc)
{
  graverComplexity = gc; //Maybe check a bound?
}

void NFold::solve()
{
  //NEED DOUBLE CHECK IF INITIAL SOLUTION WAS FOUND?

  bool done = false;

  //Keep obtaining and adding best steps until step of 0 found
  while (!done)
  {
    vector<int> graverBestStep = findGraverBestStep();

    bool done = false;

    //If returns a zero step or a step in wrong direction
    if (innerProduct(objective, graverBestStep) > -1)
    {
      done = true;
    }
    else //Good step. Add to current solution and continue
    {
      for (size_t i = 0; i < n * t; i++)
      {
        currentSolution[i] += graverBestStep[i];
      }
    }
  }
}

vector<int> NFold::findGraverBestStep()
{
  //Find infinite norm of u - l
  int infNorm = 0;
  for (size_t i = 0; i < (n * t); i++)
  {
    int absDifference = abs(upperBound[i] - lowerBound[i]);
    if (absDifference > infNorm)
    {
      infNorm = absDifference;
    }
  }

  int upperLimit = int(floor(log2(infNorm)));
  int lambda = 1;

  vector<int> bestStep(n * t, 0);
  //int bestObjective
  for (size_t i = 0; i < upperLimit; i++)
  {
    vector<int> goodStep = findGoodStep(lambda);

    if (innerProduct(objective, goodStep) > -1)
    {
      break;
    }

    int exhaustedLambda = INT_MAX;



    //Exhaust the lamda and see if it is better step than ones seen so far
    //May have a break if w dot good_step >= 0?
    lambda *= 2; //TODO: Maybe generalize this to other c-apx values
  }

  return bestStep;
}

vector<int> NFold::findGoodStep(int lambda)
{
  try {
    //Create set Gamma and solve ILP for each
    GRBEnv env = GRBEnv();
    env.set("LogFile", "info.log");
    env.start();
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_LogToConsole, 0);

    //Create bounds for solution. Use double type to match addVars() requirements
    //graverComplexity * lambda from Altmanova implementation. Need to address theory
    double *lb = new double [n * t];
    for (size_t i = 0; i < (n * t); i++)
    {
      lb[i] = ceil(max(-1.0 * graverComplexity * lambda, 1.0 * (lowerBound[i] - currentSolution[i]) / lambda));
    }

    double *ub = new double [n * t];
    for (size_t i = 0; i < (n * t); i++)
    {
      ub[i] = floor(min(1.0 * graverComplexity * lambda, 1.0 * (upperBound[i] - currentSolution[i]) / lambda));
    }

    double *obj = new double [n * t];
    for (size_t i = 0; i < (n * t); i++)
    {
      obj[i] = 1.0 * objective[i];
    }

    char *vTypesInt = new char [n * t];
    for (size_t i = 0; i < (n * t); i++)
    {
      vTypesInt[i] = GRB_INTEGER;
    }

    GRBVar *h = model.addVars(lb, ub, obj, vTypesInt, NULL, (n * t));

    //Add constraints from constraint matrix
    for (size_t i = 0; i < r + n * s; i++)
    {
      GRBLinExpr constraint;
      constraint.addTerms(&constraintMatrix[i][0], h, n * t);
      model.addConstr(constraint, GRB_EQUAL, 0);//1.0 * b[i]);
    }

    //Add constraints for l1 norm less than graver complexity and minimize
    double *lbZero = new double [n * t];
    double *ubInfinity = new double [n * t];
    double *arrayOnes = new double [n * t]; //For graverConstraint below
    for (size_t i = 0; i < n * t; i++)
    {
      lbZero[i] = 0;
      ubInfinity[i] = GRB_INFINITY;
      arrayOnes[i] = 1;
    }

    GRBVar *positive = model.addVars(lbZero, ubInfinity, 0, vTypesInt, NULL, (n * t));
    GRBVar *negative = model.addVars(lbZero, ubInfinity, 0, vTypesInt, NULL, (n * t));

    for (size_t i = 0; i < n * t; i++)
    {
      model.addConstr(positive[i] - negative[i], GRB_EQUAL, h[i]);
    }

    GRBLinExpr graverConstraint;
    graverConstraint.addTerms(arrayOnes, positive, n * t);
    graverConstraint.addTerms(arrayOnes, negative, n * t);
    model.addConstr(graverConstraint, GRB_LESS_EQUAL, graverComplexity);

    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.optimize();

    vector<int> zeros(n * t, 0);
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
      double objectiveValue = model.get(GRB_DoubleAttr_ObjVal);
      cout << "Objective value: " << objectiveValue << endl;

      if (objectiveValue > -1)
      {
        return zeros; //Step found is non-improving.
      }
      else
      {
        vector<int> goodStep;
        goodStep.resize(n * t);
        for (size_t i = 0; i < n * t; i++)
        {
          goodStep[i] = int(h[i].get(GRB_DoubleAttr_X));
        }

        return goodStep;
      }
    }

    //Infeasible
    return zeros;
  }
  catch(GRBException e)
  {
    cout << "Error in findGoodStep() with error code: " << e.getErrorCode() << endl;
    cout << "Error Message: " << e.getMessage() << endl;
    cout << "Exiting program..." << endl;
    exit(-1);
  }
  catch(...)
  {
    cout << "Exception during optimization in findGoodStep()" << endl;
    cout << "Exiting program..." << endl;
    exit(-1);
  }
}

int NFold::innerProduct(const std::vector<int> &v1, const std::vector<int> &v2)
{
  if (v1.size() != v2.size())
  {
    cout << "Error in call to innerProduct. Vector sizes do not match." << endl;
    exit(-1);
  }

  int prod = 0;
  for (size_t i = 0; i < v1.size(); i++)
  {
    prod += v1[i] * v2[i];
  }

  return prod;
}
