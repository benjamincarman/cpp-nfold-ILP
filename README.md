# C++ Implementation for n-fold ILP Solver

#### **Code by:** Benjamin Carman
#### **Algorithm proposed by:** Altmanová, Knop, and Koutecký
#### *Implemented as part of a DAAD Rise-funded summer internship at Universität Bonn*

## Description

### Overview

This program implements a C++ version of the algorithm proposed by Altmanová, Knop, and Koutecký. Their equivalent implementation using SageMath can be found at https://github.com/katealtmanova/nfoldexperiment.

This algorithm is a solver for Integer Linear Programs (ILPs) of a special "*n*-fold" structure. This structure comes from the constraint matrix being of a special form where the only nonzero entries exist in blocks along the top and along the main diagonal of the matrix. For example, as in the following:


          [ [T   T   T   T   T   T],
            [D   0   0   0   0   0],
            [0   D   0   0   0   0],
      A =   [0   0   D   0   0   0],
            [0   0   0   D   0   0],
            [0   0   0   0   D   0],
            [0   0   0   0   0   D] ]

Where T is a submatrix with *r* rows and *t* columns and D is a submatrix with *s* rows and *t* columns. A complete description of this special structure can be found in Altmanová, et al.'s paper: https://arxiv.org/abs/1802.09007.

This program utilizes this special structure to find augmenting steps from an initial feasible solution until an optimum solution is reached.

## Usage

### Prerequisites

* A C++17 or newer compiler
* A Gurobi installation with a valid and current license
  * Academic licenses available. More information at https://www.gurobi.com/academia/academic-program-and-licenses/

### Compiling and Running NFoldSolve

1. If you are new to Gurobi, follow their documentation to get Gurobi installed and working with a C++ API: https://www.gurobi.com/documentation
2. Ensure that the environment variable GRB_LICENSE_FILE is set to point to your gurobi.lic file by adding it to your .bashrc file.
3. Ensure that the macros INC and CPPLIB in the Makefile correctly point to the include/ and lib/ directories of your Gurobi installation.
4. Finally to run all test examples, simply run:
        make
        ./NFoldSolve tests/*.txt

### Running Your Own Test Cases

Running your own test cases should be no problem. You can write your own main program to construct the NFold class with parameters of your choosing, or you can simply write a test case in a text file and run it through the existing main program. In order to do this, ensure the text file has the following form:

        n             <-- Where n is the dimension of the ILP
        r             <-- Where r is number of rows of the top submatrix
        s             <-- Where s is the number of rows of the diagonal matrix
        t             <-- Where t is the number of columns of the submatrices

        obj[1]        <-- Where obj is the objective vector
        obj[2]
        .
        .
        .
        obj[n*t]

        l[1]          <-- Where l is the lower bound vector
        l[2]
        .
        .
        .
        l[n*t]

        u[1]
        u[2]
        .
        .
        .
        u[n*t]          <-- Where u is the upper bound vector

        top[1,1] top[1,2] ... top[1,t]        <-- Where top is the top submatrix
        top[2,1] top[2,2] ... top[2,t]
        .             ...   .
        .             ...   .
        .             ...   .
        top[r,1] top[r,2] ... top[r,t]

        diag[1,1] diag[1,2] ... diag[1,t]     <-- Where diag is the diagonal submatrix
        diag[2,1] diag[2,2] ... diag[2,t]
        .             ...   .
        .             ...   .
        .             ...   .
        diag[s,1] diag[s,2] ... diag[s,t]

        b[1]          <-- Where b is the right hand side vector
        b[2]
        .
        .
        .
        b[r+n*s]

        i[1]          <-- Where i is an initial feasible feasible solution
        i[2]              * NOTE: This is optional and can be omitted. *
        .
        .
        .
        i[n*t]

Where each parameter is in the above order and whitespace delimited. Ensure that the test case parameters are all of valid size and form or else the test case will not be accepted.
