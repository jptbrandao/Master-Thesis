#include <iostream>
#include "greedyTest.hpp"
#include <cmath>
#include <iomanip>

double evalGreedySol( Data &myD, Solution &sol )
{
  double total = 0;

  for ( int i = 0 ; i < myD.N ; ++i )
  {
    total += sol[i] *myD.c[i];
  }

  return total;
}


bool checkSolutionIsValid(Data &myD, Solution &sol);

double euclideanDist( Solution &x, Solution &y, int n)
{
  double total = 0;

  for ( int i = 0 ; i < n ; ++i )
  {
    double aux = x[i] - y[i];
    total += aux*aux;
  } 

  return sqrt(total);
}

double maxDist( Solution &x, Solution &y, int n)
{
  double max = 0;

  for ( int i = 0 ; i < n ; ++i )
  {
    double aux = fabs(x[i] - y[i]);
    if ( max < aux ) max = aux;
  }

  return max;
}

double avgDist( Solution &x, Solution &y, int n )
{
  double avg = 0;
  
  for ( int i = 0 ; i < n ; ++i )
  {
    double aux = fabs(x[i] - y[i]);
    avg += aux;
  }

  avg = avg/n;
  return avg;
}

void printSolDiff( Solution &x, Solution &y, int n )
{
  for ( int i = 0 ; i < n ; ++i )
  {
    cout << i << " - " << x[i] - y[i] << endl;
  }
}

double tolerance = 0.001;
bool testGreedy() 
{
  cout << "Testing Simple Greedy Method " << endl;
  // Test Parameters 
  int numRuns = 100;         // Number of Runs per test
  int numTests = 100;        // Number of tests
  int numVarIncrease = 5;   // Increase in number of variables between one test and the next
  int count = 0;

  // Test Setup
  unsigned int seed = time(NULL);
  int numVar = 0; 
  
  // Test Metrics 
  int ownEqualMosek = 0;    // Number of runs where solutions were equal
  int ownBetterMosek = 0;   // Number of runs where own was better than mosek
  int ownWorseMosek = 0;    // Number of runs where own was worse than mosek
  int ownInvalid = 0;       // Number of runs where own's solution was invalid
  int mosekInvalid = 0;     // Number of runs where mosek's solution was invalid
  int maxGrt1Perc = 0;      // Number of runs where max difference in solution differs more than 1%
  int avgGrt1Perc = 0;      // Number of runs where average difference in solution differs by more than 1% from own solution
  int eucGrt1Perc = 0;      // Number of runs where euclidean difference in solution differs by more than 1% from own solution

  // Test Execution
  for ( int test = 0 ; test < numTests ; ++test )
  {
    cout << endl << " -------- TEST: " << test+1 << "--------- " << endl;
    numVar += numVarIncrease;
    for ( int run = 0 ; run < numRuns ; ++run )
    {
      //cout << endl << " RUN: " << run+1 << endl;
      Data myD(numVar, seed);

      // Run own solution
      //cout << "Running own Greedy Solution" << endl;
      Solution ownSol(myD.N);
      simpleGreedy(myD, ownSol);
      //cout << "Finished running own Greedy Solution" << endl;

      // Run Mosek's
      //cout << "Running Mosek Solution " << endl;
      Solution mosekSol(myD.N);
      mosekSol = solveGreedyWithMosek(myD);
      //cout << "Finished running Mosek Solution " << endl;

      //cout << "Evaluating and comparing both solutions " << endl;
      // Compare Solutions

      // check if solutions are valid
      bool isOwnValid = true;
      if ( checkSolutionIsValid(myD, ownSol) == false )
      {
        //cout << "Own Invalid" << endl;
        ownInvalid++;
        isOwnValid = false;
      }

      bool isMosekValid = true;
      if ( checkSolutionIsValid(myD, mosekSol) == false )
      {
        //cout << "Mosek Invalid" << endl;
        mosekInvalid++;
        isMosekValid = false;

        double euclD = euclideanDist(ownSol, mosekSol, myD.N);
        double maxD = maxDist(ownSol, mosekSol, myD.N);
        double avgD = avgDist(ownSol, mosekSol, myD.N);
        cout << endl << "Count: " << count << " || Eucl: " << euclD << " || Max: " << maxD << " || Avg: " << avgD << endl;
        //cout << "OWN - MOSEK" << endl;
        //printSolDiff( ownSol, mosekSol, myD.N);
        //cout << "TotalDiff: " << fabs(ownTotal-mosekTotal) << " || PercDif: " << percentDiff << " || Own - Mosek " << ownTotal << " - " << mosekTotal <<  endl;
      }

      /*if ( isOwnValid == false )
      {
        //cout << endl << "Own Invalid -- Test x Run : " << test << " x " << run << endl;
      }

      if ( isMosekValid == false ) 
      {
        //cout << endl << "Mosek Invalid -- Test x Run : " << test << " x " << run << endl;
      }*/

      double ownTotal = evalGreedySol(myD, ownSol);
      double mosekTotal = evalGreedySol(myD, mosekSol);
      
      double percentDiff = fabs((ownTotal - mosekTotal)/mosekTotal); 
      double absDiff = fabs(ownTotal - mosekTotal);
      
      if ( percentDiff < tolerance || absDiff < tolerance )
      {
        //cout << "Same answer" << endl;
        ownEqualMosek++;
      }
      else 
      {
        if ( ownTotal < mosekTotal ) 
        {
          //cout << "Own better" << endl;
          ownBetterMosek++;
        }
        else
        {
          //cout << "Mosek better" << endl;
          ownWorseMosek++;
        }
      }

      // Comment following if running new tests
      if ( count == -1 )
      {
        cout << endl << "OWN: " << ownTotal << endl;
        ownSol.printSol();

        cout << endl << "MOSEK: " << mosekTotal  << endl;
        mosekSol.printSol();
      }
      count++;

      //cout << "Finished evaluating and comparing both solutions " << endl;
    }
  }
  
  cout << "Number of Same Answer: " << ownEqualMosek << endl;
  cout << "Number of Own Better: " << ownBetterMosek << endl;
  cout << "Number of Own Worse: " << ownWorseMosek << endl;

  cout << endl << "Number of Own Invalid: " << ownInvalid << endl;
  cout << "Number of Mosek Invalid: " << mosekInvalid << endl;

  cout << "Total Tests: " << numTests << endl;
  cout << "Total Runs: " << numRuns << endl;
  cout << "Total solved: " << numTests*numRuns << endl;
  cout << "SEED: " << seed << endl;
}

bool checkSolutionIsValid(Data &myD, Solution &sol)
{
  double total = 0;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    // non negativity
    if ( sol[i] < 0 )
    {
      cout << "Violation: Nonnegativity, index - " << i << " || value - " << sol[i] << endl;
      return false;
    }

    total += sol[i];
    
    // nested constraints;
    double upperCstTolerance = tolerance * myD.upperCst[i];
    if ( myD.nestCst[i] == true && total > myD.upperCst[i] + upperCstTolerance ) 
    {
      cout << "Violation: Upper nested, index - " << i << " || original v result vs PercDiff " << myD.upperCst[i] << " vs " << total << " vs " << (fabs(myD.upperCst[i] - total)/myD.upperCst[i]) << endl; 
      return false;
    }
  }

  // resource constraint
  if ( fabs( myD.B1 - total )/myD.B1 > tolerance ) 
  {
    cout << "Violation: B1 resource - original vs result vs PercDiff- " << myD.B1 << " vs " << total << " vs " << (fabs(myD.B1 - total)/myD.B1) << endl;
    return false;
  }
  

  return true;

}
