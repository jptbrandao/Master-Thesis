#include <iostream>
#include "greedy2Test.hpp"
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <math.h>

double evalGreedySol( Data &myD, Solution &sol )
{
  double total = 0;

  for ( int i = 0 ; i < myD.N ; ++i )
  {
    total += sol[i] *myD.c[i];
  }

  return total;
}

double evalGreedySolWithLagr( Data &myD, Solution &sol, double lambda )
{
  double total = 0;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    total += sol[i]*(myD.c[i] + lambda*myD.p[i]);
  }

  return total;
}

void printSolWithLagr( Solution &x, double lambda, Data &myD)
{
  int j = 0;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    cout << "Index: " << i << " - amount: " << x[i] << " - cost: " << myD.c[i] + lambda*myD.p[i] << " - UBound: " << myD.upperCst[j] << endl;

    if ( myD.nestCst[i] == true ) 
      j++;
  }
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

void printVector ( vector<double> &x )
{
  cout << endl;
  for ( int i = 0 ; i < x.size() ; ++i )
  {
    cout << i << " - " << x[i] << endl;
  }
}

vector<double> getLambdaPoints( Data &myD ) 
{
  // lambda = c_i - c_j / p_j - p_i

  vector<double> lambdas;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    for ( int j = i+1 ; j < myD.N ; ++j )
    {
      lambdas.push_back( (myD.c[i] - myD.c[j])/(myD.p[j]-myD.p[i]) );
    }
  }

  sort(lambdas.begin(), lambdas.end());

  for ( int i = 0 ; i < lambdas.size() - 1 ; ++i )
  {
    if ( lambdas[i] == lambdas[i+1] )
    {
      lambdas.erase( lambdas.begin()+i+1 );
      i--;
    }
  }

  for ( int i = 0 ; i < lambdas.size(); ++i )
  {
    if ( isinf(lambdas[i]) ) 
    {
      lambdas.erase(lambdas.begin()+i);
      i--;
      continue;
    }

    if ( isnan(lambdas[i]) )
    {
      lambdas.erase(lambdas.begin()+i);
      i--;
      continue;
    }
  } 

  return lambdas;
}

void solveGreedy( Data &myD, Solution &ownSol, double lambda)
{
  Solution left(myD.N), right(myD.N);
  greedyStrong( myD, lambda, left, right);
  convexSum(myD, left, right, ownSol);
}

double tolerance = 0.01;
bool testGreedy() 
{
  cout << "Testing complicated Greedy Method " << endl;
  // Test Parameters 
  int numRuns = 10;         // Number of Runs per test
  int numTests = 10;        // Number of tests
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

  int equalObjDiffSol = 0;
  int mosekObjGrterDiffSol = 0;
  int ownObjgGrterDiffSol = 0;

  // Test Execution
  for ( int test = 0 ; test < numTests ; ++test )
  {
    cout << endl << " -------- TEST: " << test+1 << "--------- " << endl;
    numVar += numVarIncrease;
    for ( int run = 0 ; run < numRuns ; ++run )
    {
      cout << endl << " RUN: " << run+1 << endl;
      Data myD(numVar, seed);

      // get lambda points
      //cout << "Getting lambdas... " << endl;
      vector<double> lambdas = getLambdaPoints(myD);
      //cout << "Finished getting lambdas" << endl;

      //cout << "Testing on all lambdas" << endl;
      for ( int L = 0 ; L < lambdas.size() ; ++L )
      {
        double lambda = lambdas[L];
      
        // Run own solution
        //cout << "Running own Greedy Solution" << endl;
        Solution ownSol(myD.N);
        solveGreedy(myD, ownSol, lambda);

        // Run Mosek's
        //cout << "Running Mosek Solution " << endl;
        Solution mosekSol(myD.N);
        mosekSol = solveGreedy2WithMosek(myD, lambda);

        // Compare Solutions
        // check if solutions are valid
        bool isOwnValid = true;
        if ( checkSolutionIsValid(myD, ownSol) == false )
        {
          cout << "Own Invalid" << endl;
          ownInvalid++;
          isOwnValid = false;
        }

        bool isMosekValid = true;
        if ( checkSolutionIsValid(myD, mosekSol) == false )
        {
          cout << "Mosek Invalid: " << lambda << endl;
          mosekInvalid++;
          isMosekValid = false;

          //cout << endl << "Count: " << count << " || Eucl: " << euclD << " || Max: " << maxD << " || Avg: " << avgD << endl;
          //cout << "OWN - MOSEK" << endl;
          //printSolDiff( ownSol, mosekSol, myD.N);
        }

        /*if ( isOwnValid == false )
        {
          //cout << endl << "Own Invalid -- Test x Run : " << test << " x " << run << endl;
        }

        if ( isMosekValid == false ) 
        {
          //cout << endl << "Mosek Invalid -- Test x Run : " << test << " x " << run << endl;
        }*/
        double euclD = euclideanDist(ownSol, mosekSol, myD.N);
        double maxD = maxDist(ownSol, mosekSol, myD.N);
        double avgD = avgDist(ownSol, mosekSol, myD.N);

        // Auxilary variable to determine the difference which metrics had greater than 1% difference
        // +1 to euclD ; +2 to maxD ; +4 to avgD
        int aux = 0; 

        double ownTotal = evalGreedySolWithLagr(myD, ownSol, lambda);
        double mosekTotal = evalGreedySolWithLagr(myD, mosekSol, lambda);
      
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
            cout << count << endl;
            cout << "TotalDiff: " << fabs(ownTotal-mosekTotal) << " || PercDif: " << percentDiff << " || Own - Mosek " << ownTotal << " - " << mosekTotal <<  endl;
            cout << "EuclDist: " << euclD << " || MaxDist: " << maxD << " || AvgDist: " << avgD << endl << endl;
          }
        }

        // Comment following if running new tests
        if ( count == -1 )
        {
          cout << endl << "OWN: " << ownTotal << endl;
          printSolWithLagr( ownSol, lambda, myD);
          cout << endl;

          cout << endl << "MOSEK: " << mosekTotal  << endl;
          printSolWithLagr(mosekSol, lambda, myD);
          
        }
        count++;
      }
      //cout << "Finished testing on all lambdas" << endl;

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
