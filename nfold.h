/*******************************************************************************
 *
 * File: nfold.h
 *
 * Author: Benjamin Carman
 *
 * Date: June 4, 2019
 *
 * Description: An header file of the NFold class which allows for the
 *              the creation of an ILP instance with the n-fold form. This class
 *              verifies the ILP's structure and can solve it by finding
 *              augmenting steps using Gurobi to solve auxillary programs until
 *              an optimum solution is reached. See documentation, README, and
 *              paper by Altmanová, Knop, and Koutecký for more information.
 *
 ******************************************************************************/

#ifndef NFOLD_H
#define NFOLD_H

#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <optional>

class NFold {
public:
  //Constructors
  NFold(GRBEnv *e, unsigned int gc = 2);
  NFold(GRBEnv *e, unsigned int n1, std::vector<int> objective1, std::vector<int> lowerBound1,
        std::vector<int> upperBound1, std::vector<std::vector<int> > topMatrix1,
        std::vector<std::vector<int> > diagMatrix1, std::vector<int> b1,
        unsigned int gc = 2, std::optional<std::vector<int> > initial1 = std::nullopt);

  //Setter
  void setGraverComplexity(unsigned int gc);

  //Getter
  int getOptimizedObjectiveValue() const;
  std::vector<int> getOptimizedSolution() const;

  //I/O
  void inputState(std::istream &ins);
  void buildConstraintMatrix();
  void outputState(std::ostream &outs) const;

  //Error checking
  std::string checkDataValidity() const;
  std::string checkInitialSolution() const;

  //Solver
  bool solve();
  std::vector<int> findGraverBestStep();
  std::vector<int> findGoodStep(int lambda);
  void findInitialFeasibleSolution();

  //Helper
  template <class T>
  T innerProduct(const std::vector<T> &v1, const std::vector<T> &v2) const;
  template <class T>
  void writeVector(const std::vector<T> &v, std::ostream &outs) const;
  template <class T>
  void writeMatrix(const std::vector<std::vector<T> > &m, std::ostream &outs) const;


private:
  GRBEnv *env;
  unsigned int n;
  unsigned int r;
  unsigned int s;
  unsigned int t;
  std::vector<int> objective;
  std::vector<int> lowerBound;
  std::vector<int> upperBound;
  std::vector<std::vector<int> > topMatrix; //Matrix along top of constraintMatrix
  std::vector<std::vector<int> > diagMatrix; //Matrix along diagonal of constraintMatrix
  std::vector<std::vector<double> > constraintMatrix;
  std::vector<int> b;
  std::vector<int> currentSolution;
  unsigned int graverComplexity;
  std::optional<std::vector<int> > initial;
  bool solved;
  bool safelyInstantiated;
};

#endif
