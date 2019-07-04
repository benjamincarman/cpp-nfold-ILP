/*******************************************************************************
 *
 * File: nfold.cc
 *
 * Author: Benjamin Carman
 *
 * Date: June 4, 2019
 *
 * Description: An implementation file of the NFold class which allows for the
 *              the creation of an ILP instance with the n-fold form. This class
 *              verifies the ILP's structure and can solve it by finding
 *              augmenting steps using Gurobi to solve auxillary programs until
 *              an optimum solution is reached. See documentation, README, and
 *              paper by Altmanová, Knop, and Koutecký for more information.
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

NFold::NFold(GRBEnv *e, unsigned int gc)
{
  env = e;
  solved = false;
  safelyInstantiated = false;
  initial = std::nullopt;
  graverComplexity = gc;
}

NFold::NFold(GRBEnv *e, unsigned int n1, std::vector<int> objective1, std::vector<int> lowerBound1,
             std::vector<int> upperBound1, std::vector<std::vector<int> > topMatrix1,
             std::vector<std::vector<int> > diagMatrix1, std::vector<int> b1,
             unsigned int gc, std::optional<std::vector<int> > initial1)
{
  env = e;
  n = n1;
  r = topMatrix1.size();
  s = diagMatrix1.size();
  t = topMatrix1[0].size();
  objective = objective1;
  lowerBound = lowerBound1;
  upperBound = upperBound1;
  topMatrix = topMatrix1;
  diagMatrix = diagMatrix1;

  b = b1;
  currentSolution.resize(n * t);

  string error = checkDataValidity();
  if (error != "")
  {
    cout << error << endl;
    exit(-2);
  }

  buildConstraintMatrix();

  graverComplexity = gc;

  if (initial1 == std::nullopt)
  {
    cout << "No initial solution provided." << endl;

    //Commence finding initial solution here
    findInitialFeasibleSolution();
    if (initial == std::nullopt)
    {
      cout << "No feasible solution for this instance exists." << endl;
      exit(-2);
    }
  }
  else
  {
    initial = initial1.value();
  }

  //Check validity of initial solution
  error = checkInitialSolution();
  if (error != "")
  {
    cout << error << endl;
    exit(-2);
  }

  solved = false;

  //Error checks pass, object is safely instantiated
  safelyInstantiated = true;
}

void NFold::setGraverComplexity(unsigned int gc)
{
  if (gc >= 2)
  {
    solved = false;
    graverComplexity = gc;
  }
  else
  {
    cout << "Error: Invalid graver complexity." << endl;
    exit (-2);
  }
}

int NFold::getOptimizedObjectiveValue() const
{
  if (solved)
  {
    return innerProduct(objective, currentSolution);
  }
  else
  {
    cout << "Cannnot get optimized objective value. Instance not solved." << endl;
    exit(-2);
  }
}

vector<int> NFold::getOptimizedSolution() const
{
  if (solved)
  {
    return currentSolution;
  }
  else
  {
    cout << "Cannot get optimized solution. Instance not solved." << endl;
    exit(-2);
  }
}

void NFold::inputState(istream &ins)
{
  solved = false;

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

  for (size_t i = 0; i < r + (n * s); i++)
  {
    ins >> b[i];
  }

  if (ins.eof())
  {
    cout << "Error: Input file does not contain enough input data. NFold not instantiated." << endl;
    exit(-2);
  }

  //Call error check function
  string error = checkDataValidity();
  if (error != "")
  {
    cout << error << endl;
    exit(-2);
  }

  buildConstraintMatrix();

  vector<int> initialInput;
  initialInput.resize(n * t);
  for (size_t i = 0; i < n * t; i++)
  {
    ins >> initialInput[i];
  }

  if (!ins.eof())
  {
    initial = initialInput;
  }
  else
  {
    cout << "No initial solution provided." << endl;

    //Commence finding initial solution here
    findInitialFeasibleSolution();
    if (initial == std::nullopt)
    {
      cout << "No feasible solution for this instance exists." << endl;
      exit(-2);
    }
  }

  //Check validity of initial solution
  error = checkInitialSolution();
  if (error != "")
  {
    cout << error << endl;
    exit(-2);
  }

  //At this point, all is safely instantiated
  safelyInstantiated = true;
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

void NFold::outputState(std::ostream &outs) const
{
  outs << "n: " << n << endl
       << "r: " << r << endl
       << "s: " << s << endl
       << "t: " << t << endl;

  outs << endl << "Objective:" << endl;
  writeVector(objective, outs);

  outs << endl << "Lower Bound:" << endl;
  writeVector(lowerBound, outs);

  outs << endl << "Upper Bound:" << endl;
  writeVector(upperBound, outs);

  outs << endl << "Top Matrix:" << endl;
  writeMatrix(topMatrix, outs);

  outs << endl << "Diagonal Matrix:" << endl;
  writeMatrix(diagMatrix, outs);

  outs << endl << "An Matrix:" << endl;
  writeMatrix(constraintMatrix, outs);

  outs << endl << "RHS (b): " << endl;
  writeVector(b, outs);

  outs << endl << "Solution: " << endl;
  writeVector(currentSolution, outs);

  outs << endl << "Objective Value: " << innerProduct(objective, currentSolution) << endl;
}

string NFold::checkDataValidity() const
{
  //Check sizes of each vector and matrix
  if (objective.size() != n * t)
  {
    return "Instantiation Error: Objective vector is not of size (n * t)";
  }
  if (lowerBound.size() != n * t)
  {
    return "Instantiation Error: Lower Bound vector is not of size (n * t)";
  }
  if (upperBound.size() != n * t)
  {
    return "Instantiation Error: Upper Bound vector is not of size (n * t)";
  }
  if (topMatrix.size() != r)
  {
    return "Instantiation Error: Top matrix does not have r rows.";
  }
  for (size_t i = 0; i < r; i++)
  {
    if (topMatrix[i].size() != t)
    {
      return "Instantiation Error: Top matrix does not have t columns.";
    }
  }
  if (diagMatrix.size() != s)
  {
    return "Instantiation Error: Diagonal matrix does not have s rows.";
  }
  for (size_t i = 0; i < s; i++)
  {
    if (diagMatrix[i].size() != t)
    {
      return "Instantiation Error: Diagonal matrix does not have t columns.";
    }
  }
  if (b.size() != r + n * s)
  {
    return "Instantiation Error: RHS vector is not of size (r + n * s)";
  }

  //Check that lower <= upper
  for (size_t i = 0; i < n * t; i++)
  {
    if (lowerBound[i] > upperBound[i])
    {
      return "Instantiation Error: Lower bounds must be less than or equal to upper bounds.";
    }
  }

  return "";
}

string NFold::checkInitialSolution() const
{
  if (initial == std::nullopt)
  {
    return "Initial Solution Instantiation Error: No initial solution provided.";
  }

  vector<int> inSol = initial.value();

  if (inSol.size() != n * t)
  {
    return "Initial Solution Instantiation Error: Initial solution is not of size (n * t).";
  }

  for (size_t i = 0; i < n * t; i++)
  {
    if (inSol[i] < lowerBound[i] || inSol[i] > upperBound[i])
    {
      return "Initial Solution Instantiation Error: Initial solution is not within upper and lower bounds.";
    }
  }

  //Check feasibility of initial solution
  vector<int> crossProduct;
  crossProduct.resize(r + n * s);

  //Create inSol vector of doubles to work with constraint matrix
  vector<double> inSolDouble;
  inSolDouble.resize(n * t);
  for (size_t i = 0; i < n * t; i++)
  {
    inSolDouble[i] = 1.0 * inSol[i];
  }

  //Calculate cross product
  for (size_t i = 0; i < r + n * s; i++)
  {
    int prod = int(innerProduct(constraintMatrix[i], inSolDouble));
    crossProduct[i] = prod;
  }

  if (crossProduct != b)
  {
    return "Initial Solution Instantiation Error: Provided initial solution is infeasible.";
  }

  //No error found, return empty string
  return "";
}


bool NFold::solve()
{
  if (initial == std::nullopt)
  {
    cout << "No initial solution. Unable to solve." << endl;
    solved = false;
    return false;
  }

  if (solved == true)
  {
    cout << "Instance already solved." << endl;
    return true;
  }

  currentSolution = initial.value();

  cout << "Begin solve with GC: " << graverComplexity << endl;
  bool done = false;
  //Keep obtaining and adding best steps until step of 0 found
  while (!done)
  {
    cout << "Current solution is: " << endl << "<";
    for (size_t i = 0; i < n * t - 1; i++)
    {
      cout << currentSolution[i] << ", ";
    }
    cout << currentSolution[n * t - 1] << ">" << endl;

    vector<int> graverBestStep = findGraverBestStep();

    //If returns a zero step or a step in wrong direction
    if (innerProduct(objective, graverBestStep) > -1)
    {
      done = true;
    }
    else //Good step. Add to current solution and continue
    {
      cout << "Graver best step found:" << endl << "<";
      for (size_t i = 0; i < n * t - 1; i++)
      {
        cout << graverBestStep[i] << ", ";
      }
      cout << graverBestStep[n * t - 1] << ">" << endl;

      for (size_t i = 0; i < n * t; i++)
      {
        currentSolution[i] += graverBestStep[i];
      }
    }
  }

  solved = true;
  return true;
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
  int bestObjectiveValue = 0;
  for (size_t i = 0; i < upperLimit; i++)
  {
    vector<int> goodStep = findGoodStep(lambda);
    cout << "Good Step Found: " << endl << "<";
    for (size_t i = 0; i < n * t - 1; i++)
    {
      cout << goodStep[i] << ", ";
    }
    cout << goodStep[n * t - 1] << ">" << endl;

    if (innerProduct(objective, goodStep) > -1)
    {
      cout << "No good step found with lambda: " << lambda << endl;
      break;
    }

    int exhaustedLambda = INT_MAX;
    for (size_t i = 0; i < n * t; i++)
    {
      if (goodStep[i] != 0)
      {
        double lowerFeasibleStepLength = (1.0 *(lowerBound[i] - currentSolution[i]))/goodStep[i];
        double upperFeasibleStepLength = (1.0 *(upperBound[i] - currentSolution[i]))/goodStep[i];

        int bestFeasibleStepLength = int(floor(max(lowerFeasibleStepLength, upperFeasibleStepLength)));

        if (exhaustedLambda > bestFeasibleStepLength)
        {
          exhaustedLambda = bestFeasibleStepLength;
        }
      }
    }

    //Just in case
    if (exhaustedLambda == INT_MAX || exhaustedLambda < lambda)
    {
      exhaustedLambda = lambda;
    }

    for (size_t i = 0; i < n * t; i++)
    {
      goodStep[i] *= exhaustedLambda;
    }
    int currentObjectiveValue = innerProduct(objective, goodStep);

    if (exhaustedLambda != lambda)
    {
      cout << "Lambda changed! Exhausted Lambda: " << exhaustedLambda << endl;
    }
    cout << "Objective value for lambda = " << exhaustedLambda << ", is: " << currentObjectiveValue << endl;
    if (currentObjectiveValue < bestObjectiveValue)
    {
      bestObjectiveValue = currentObjectiveValue;
      bestStep = goodStep;
    }

    //Exhaust the lamda and see if it is better step than ones seen so far
    lambda *= 2; //TODO: Maybe generalize this to other c-apx values
  }

  return bestStep;
}

vector<int> NFold::findGoodStep(int lambda)
{
  cout << "Finding good step with lambda: " << lambda << endl;
  try {
    //Create set Gamma and solve ILP for each
    GRBModel model = GRBModel(*env);
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
      //cout << "Objective value: " << objectiveValue << endl;

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

void NFold::findInitialFeasibleSolution()
{
  cout << "Trying to find initial feasible solution." << endl;

  vector<int> iObjective;
  vector<int> iLowerBound;
  vector<int> iUpperBound;
  vector<vector<int> > iTopMatrix(r, vector<int>(2 * (r + s) + t, 0));
  vector<vector<int> > iDiagMatrix(s, vector<int>(2 * (r + s) + t, 0));
  vector<int> iInitial;

  //Set sizes of vectors
  iObjective.resize(n * (2 * (r + s) + t));
  iLowerBound.resize(n * (2 * (r + s) + t));
  iUpperBound.resize(n * (2 * (r + s) + t));
  iInitial.resize(n * (2 * (r + s) + t));

  size_t newBrickSize = t + 2 * (r + s);

  //Initialize objective
  for (size_t brick = 0; brick < n; brick++)
  {
    for (size_t i = 0; i < t; i++)
    {
      iObjective[brick * newBrickSize + i] = 0;
    }
    for (size_t i = t; i < newBrickSize; i++)
    {
      iObjective[brick * newBrickSize + i] = 1;
    }
  }

  //Initialize lower bound
  for (size_t brick = 0; brick < n; brick++)
  {
    for (size_t i = 0; i < t; i++)
    {
      iLowerBound[brick * newBrickSize + i] = lowerBound[i];
    }
    for (size_t i = t; i < newBrickSize; i++)
    {
      iLowerBound[brick * newBrickSize + i] = 0;
    }
  }

  //Get bMax value for upper bound
  int bMax = 0;
  for (size_t i = 0; i < r + n * s; i++)
  {
      if (abs(b[i]) > bMax)
      {
        bMax = abs(b[i]);
      }
  }

  //Initialize upper bound
  for (size_t brick = 0; brick < n; brick++)
  {
    for (size_t i = 0; i < t; i++)
    {
      iUpperBound[brick * newBrickSize + i] = upperBound[i];
    }
    for (size_t i = t; i < newBrickSize; i++)
    {
      iUpperBound[brick * newBrickSize + i] = bMax;
    }
  }

  //Initialize top matrix
  for (size_t i = 0; i < r; i++)
  {
    for (size_t j = 0; j < t; j++)
    {
      iTopMatrix[i][j] = topMatrix[i][j];
    }

    iTopMatrix[i][i + t] = 1;
    iTopMatrix[i][i + t + r] = -1;
  }

  //Initialize the diagonal matrix
  for (size_t i = 0; i < s; i++)
  {
    for (size_t j = 0; j < t; j++)
    {
      iDiagMatrix[i][j] = diagMatrix[i][j];
    }

    iDiagMatrix[i][i + t + r + r] = 1;
    iDiagMatrix[i][i + t + r + r + s] = -1;
  }

  //Iniitialize the initial solution for this subproblem
  size_t index = 0;

  for (size_t i = 0; i < t; i++)
  {
    if (lowerBound[i] != INT_MIN)
    {
      iInitial[i] = lowerBound[i];
    }
    else if (upperBound[i] != INT_MIN)
    {
      iInitial[i] = upperBound[i];
    }
    else
    {
      iInitial[i] = 0;
    }
  }
  index = t;
  for (size_t i = 0; i < r; i++)
  {
    if (b[i] >= 0)
    {
      iInitial[i + index] = b[i];
      iInitial[i + index + r] = 0;
    }
    else
    {
      iInitial[i + index] = 0;
      iInitial[i + index + r] = -1 * b[i];
    }
  }
  index = t + 2 * r;
  for (size_t i = 0; i < s; i++)
  {
    if (b[i + r] >= 0)
    {
      iInitial[i + index] = b[i + r];
      iInitial[i + index + s] = 0;
    }
    else
    {
      iInitial[i + index] = 0;
      iInitial[i + index + s] = -1 * b[i + r];
    }
  }
  index = newBrickSize;

  //Fill in the rest
  for (size_t sSection = 1; sSection < n; sSection++)
  {
    for (size_t i = 0; i < 2 * r + t; i++)
    {
      iInitial[i + index] = 0;
    }
    index += 2 * r + t;
    for (size_t i = 0; i < s; i++)
    {
      if (b[r + sSection * s] >= 0)
      {
        iInitial[i + index] = b[r + sSection * s];
        iInitial[i + index + s] = 0;
      }
      else
      {
        iInitial[i + index] = 0;
        iInitial[i + index + s] = -1 * b[r + sSection * s];
      }
    }

    index += 2 * s;
  }

  //Create instance
  NFold auxILP(env, n, iObjective, iLowerBound, iUpperBound,
               iTopMatrix, iDiagMatrix, b, graverComplexity, iInitial);
  auxILP.solve();
  auxILP.outputState(cout);

  //Check solution of auxILP
  vector<int> auxSolution = auxILP.getOptimizedSolution();

  bool valid = true;
  for (size_t brick = 0; brick < n; brick++)
  {
    for (size_t i = t; i < newBrickSize; i++)
    {
      if (auxSolution[brick * newBrickSize + i] != 0)
      {
        valid = false;
        break;
      }
    }
    if (!valid)
    {
      break;
    }
  }

  if (!valid)
  {
    initial = std::nullopt;
  }
  else
  {
    vector<int> finalSolution;
    finalSolution.resize(n * t);

    for (size_t brick = 0; brick < n; brick++)
    {
      for (size_t i = 0; i < t; i++)
      {
        finalSolution[brick * t + i] = auxSolution[brick * newBrickSize + i];
      }
    }
    cout << "Initial Feasible Solution Found: " << endl;
    auxILP.writeVector(finalSolution, cout);
    initial = finalSolution;
  }
}

template <class T>
T NFold::innerProduct(const std::vector<T> &v1, const std::vector<T> &v2) const
{
  if (v1.size() != v2.size())
  {
    cout << "Error in call to innerProduct. Vector sizes do not match." << endl;
    exit(-1);
  }

  T prod = 0;
  for (size_t i = 0; i < v1.size(); i++)
  {
    prod += v1[i] * v2[i];
  }

  return prod;
}

template <class T>
void NFold::writeVector(const std::vector<T> &v, std::ostream &outs) const
{
  outs << "<";
  for (size_t i = 0; i < v.size() - 1; i++)
  {
    outs << v[i] << ", ";
  }
  outs << v[v.size() - 1] << ">" << endl;
}

template <class T>
void NFold::writeMatrix(const std::vector<std::vector<T> > &m, std::ostream &outs) const
{
  outs << "[ ";
  for (size_t i = 0; i < m.size(); i++)
  {
    if (i != 0) outs << "  ";
    outs << "[";
    for (size_t j = 0; j < m[i].size() - 1; j++)
    {
      outs << m[i][j] << ", ";
    }
    outs << m[i][m[i].size() - 1];
    if (i != m.size() - 1) {outs << "]," << endl;}
    else {outs << "] ]" << endl;}
  }
}
