#include "TestQuickSort.hpp"
#include <iostream>
#include "../testSolver.hpp"
#include <algorithm>
#include "../commonFuncs.hpp"

//string filepath = "exportedData.txt";


bool checkOrder( vector<double> &v )
{
  for( int i = 0 ; i < v.size()-1 ; ++i)
    if ( v[i] > v[i+1] )
      return false;

  return true;
}

void findFailure(Data &myD, Solution &sol)
{
  double tol = 0.001; // tolerance for errors

  if ( sol.isEmpty() == true )
  {
    cout << "Infeasible" << endl; 
  }
    
  double sumR1 = 0;       // sumR1 - resource1
  double sumR2 = 0;       // sumR2 - resource2
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    sumR1 += sol[i];
    sumR2 += sol[i]*myD.p[i];
  }
    
  // Checks B1 and B2 constraints
  double resB1 = fabs(myD.B1 - sumR1) ;
  double resB2 = fabs(myD.B2 - sumR2) ;
    
  //cout << "Checking B1 and B2 " << endl;
  if ( resB1 > tol )  cout << "Failed B1" << endl;
  if ( resB2 > tol )  cout << "Failed B2 (solB2, realB2): (" << sumR2 << ", " << myD.B2 << ")" << endl;
    
  // Checks nonnegativity
  //cout << "Checking nonnegativity" << endl;
  for ( int i = 0 ; i < myD.N ; ++i )
    if ( sol[i] < 0 )
    {
      cout << "Failed nonnegativity, index: " << i << endl;
      break;
    }
    
  // Checks nested constraints
  //cout << "Checking nested " << endl;
  sumR1 = 0;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    sumR1 += sol[i];
    if ( myD.nestCst[i] && sumR1 > myD.upperCst[i] + tol)
    {
      cout << "Failed nested, index: " << i << endl;
      break;
    }
  }
}

void printSolWithCst( Solution &x, Data &myD)
{
  int j = 0;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    cout << "Index: " << i << "\t- amount: " << x[i] << "\t\t- (cost,p): (" << myD.c[i] << ", " << myD.p[i] << ")\t- UBound: " << myD.upperCst[j] << endl;

    if ( myD.nestCst[i] == true ) 
      j++;
  }
}

bool checkResults(Data &myD, Solution &mosekSol, Solution &customSol, bool writeToOutput)
{
  // Check result
  bool isCustomSuccesful = myD.isFeasible(customSol);//checkOrder(v);
  bool isMosekSuccesful = myD.isFeasible(mosekSol);
      
  // Custom failed and mosek succeeded
  if ( !isCustomSuccesful && isMosekSuccesful )
  {
    findFailure(myD, customSol);

    cout << "Own Sol: " << endl;
    printSolWithCst(customSol, myD);
    cout << "Mosek Sol: " << endl;
    printSolWithCst(mosekSol, myD);

    //if ( writeToOutput )
    //{
    //  myD.exportData(filepath);
    //}

    return false;
  }

  return true;
}

bool checkProblemFeasibility(Data &myD)
{
  vector<double> lambda;

  for ( int i = 0 ; i < myD.N ; ++i )
  {
    for ( int j = i+1 ; j < myD.N ; ++j )
    {
      if ( myD.p[i] == myD.p[j] )
        continue;

      double intersection = (myD.c[i] - myD.c[j])/(myD.p[j] - myD.p[i]);
      lambda.push_back(intersection);
    }
  }

  sort(lambda.begin(), lambda.end());

  Solution left(myD.N), right(myD.N), middle(myD.N);
  
  greedyStrong(myD, lambda[0], left, right);
  convexSum(myD, left, right, middle);
  double subgrad0 = myD.subgradient(middle);

  greedyStrong(myD, lambda[lambda.size()-1], left, right);
  convexSum(myD, left, right, middle);
  double subgradN = myD.subgradient(middle);

  return subgrad0 * subgradN < 0;
}


bool isDataOk(Data &myD)
{
  vector<pair<double, double> > costP(myD.N);
  
  for ( int i = 0 ; i < myD.N ; ++i )
    costP.push_back(make_pair(myD.c[i], myD.p[i]));

  sort(costP.begin(), costP.end());

  for ( int i = 1 ; i < myD.N ; ++i )
  {
    if ( costP[i].first == costP[i-1].first && costP[i].second == costP[i-1].second )
      return false;
  }

  return true;
}

void runQuickSortTest(const int numTests, const int numRuns, const int inputInitialSize, const int inputGrowth)
{
  // Test Setup
  int sizeOfVector = inputInitialSize;
  int numFailures = 0;

  unsigned int seed = time(NULL);
  int runCount = 0;
  vector<int> failureRuns;

  for ( int testRun = 0 ; testRun < numTests ; ++testRun )
  {
    sizeOfVector += inputGrowth;

    cout << "Input Size: " << sizeOfVector << endl;
    for ( int run = 0 ; run < numRuns ; ++run )
    {
      cout << "Run: " << run+1 << endl;
      cout << "runCount: " << runCount << endl;
      runCount++;

      // Run Setup
      vector<double> v;
      Data myD(sizeOfVector, seed);

      if ( checkProblemFeasibility(myD) == false )
      {
        cout << "Infeasible" << endl;
        continue;
      }

      if ( isDataOk(myD) == false )
      {
        cout << "Data not ok: coinciding lines" << endl;
        //continue;
      }

      myD.printData();

      for ( int i = 0 ; i < sizeOfVector; ++i )
        v.push_back(i);

      Solution customSol(myD.N);
      Solution mosekSol(myD.N);

      cout << "Running custom algorithm..." << endl;
      // Run Algorithm
      resetLambdas();
      customSort(myD, v, 0, sizeOfVector-1, customSol) ;
      cout << "Finished custom algorithm" << endl;

      cout << "Running mosek algorithm..." << endl;
      // Run Algorithm
      mosekSol = solveMosek(myD);
      cout << "Finished mosek algorithm" << endl;

      bool isOk = checkResults(myD, mosekSol, customSol, true);
     
      if ( isOk == false )
      {
        numFailures++;
        std::cout << "Fail - (testRun, run): ( " << testRun << ", " << run << ")" << std::endl;
        failureRuns.push_back(runCount);
      }
          
      v.clear();
    }
  }

  if ( numFailures == 0 )
  {
    std::cout << "Test Succesful - No Errors" << std::endl;
  }
  else
  {
    std::cout << "Test UNSUCCESFUL - seed: " << seed << std::endl;
    std::cout << "Fail Runs: " << std::endl;
    for ( int i = 0 ; i < failureRuns.size() ; ++i )
      std::cout << "run: " << failureRuns[i] << std::endl;
  }
}

//void testOnSpecificData(ios::pos_type &position)
//{
//  Data myD(filepath, position);

//  Solution mosekSol(myD.N), customSol(myD.N);
//
//  // Run Setup
//  vector<double> v;
//  for ( int i = 0 ; i < myD.N; ++i )
//    v.push_back(i+1);
//
//  // Solve custom
//  customSort(myD, v, 0, myD.N-1, customSol);
//
//  // Solve Mosek
//  mosekSol = solveMosek(myD);
//
//  cout << endl << endl;
//  checkResults(myD, mosekSol, customSol, false);
//
//  cout << "\n\nData info: "<< endl;
//  cout << "B1, B2: " << myD.B1 << ", " << myD.B2 << endl;
//
//  cout << "Mosek: " ;
//  mosekSol.printSol();
//
//  cout << "Custom: " ;
//  customSol.printSol();
//
//} 
