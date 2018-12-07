#include <stdio.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include "quicksort.hpp"
#include "testFiles/TestQuickSort.hpp"
#include "testSolver.hpp"
#include "testFiles/testWeakly.hpp"
#include "experiments.hpp"

#include "subprocedureTests/greedy2Test.hpp"

void initialTest()
{
  // Setting up test
  int numVar = 3;
  Data myD;

  // Setting up for mosek
  Solution sol(numVar);
  sol = solveMosek(myD);

  // Print mosek solution
  std::cout << "Mosek: " << std::endl;
  sol.printSol();
  std::cout << (myD.isFeasible(sol) ? "Feasible" : "Infeasible") << std::endl;

  // Setting up for custom
  Solution customSol(numVar);
  vector<double> v;
  for ( int i = 0 ; i < myD.N ; ++i )
    v.push_back(i+1);
  
  cout << "B2: " << myD.B2 << endl;
  customSort(myD, v, 0, numVar-1, customSol);

  // Print custom solution
  std::cout << "Custom: " << std::endl;
  customSol.printSol();
  std::cout << (myD.isFeasible(customSol) ? "Feasible" : "Infeasible") << std::endl;


  return;
}

int main()
{
  cout << endl << endl;
  int numTests = 5;
  int numRuns = 5;
  int inputInitialSize = 0;
  int inputGrowth = 5;

  quicktest();
  
  //runWeaklyPolyTests( numTests, numRuns, inputInitialSize, inputGrowth);

  //ios::pos_type position = 0;
  //testOnSpecificData(position);

  //cout << "Test BEGUN" << endl;
  //testGreedy();
    return 0;
}
