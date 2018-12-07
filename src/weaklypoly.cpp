#include <iostream>
#include <limits>
#include "weaklypoly.hpp"
#include "commonFuncs.hpp"

double tolz = 0.001;



int countInversions(Data &myD, double lambda)
{
  // First is c_i + lambda *p_i, second is p_i
  vector<pair<double, double> > v;

  for ( int i = 0 ; i < myD.N ; ++i )
    v.push_back(make_pair(myD.c[i] + lambda*myD.p[i], myD.p[i]));

  sort(v.begin(), v.end());

  vector<double> pi(myD.N);
  for ( int i = 0 ; i < myD.N ; ++i )
    pi[i] = v[i].second;

  return inversion(pi, 0, myD.N-1);
}


void weaklyPoly(Data &myD, double lambdaMin, double lambdaMax, Solution &x)
{
  int countLeft = countInversions(myD, lambdaMin);
  int countRight = countInversions(myD, lambdaMax);

  while ( countLeft - countRight > 1 )
  {
    double lambdaMid = lambdaMin + (lambdaMax - lambdaMin)/2;

    greedy(myD, x, lambdaMid);

    double subgradient = myD.subgradient(x);

    if ( subgradient > tolz )
    {
      lambdaMin = lambdaMid;
      countLeft = countInversions(myD, lambdaMin);
    }
    else if ( subgradient < -tolz )
    {
      lambdaMax = lambdaMid;
      countRight = countInversions(myD, lambdaMax);
    }
    else
      return;
  }

  Solution xLeft(myD.N), xRight(myD.N);
  greedy(myD, xLeft, lambdaMin);
  greedy(myD, xRight, lambdaMax);
  convexSum(myD, xLeft, xRight, x);

  return;
}

/// ----------------------------------------------
/// ---------- Other exported functions ----------
/// ----------------------------------------------

pair<double, double> findBigLambdaInterval(Data &myD)
{
  vector<pair<double, int> > pAndIndex(myD.N);

  for ( int i = 0 ; i < myD.N ; ++i )
    pAndIndex[i] = make_pair(myD.p[i], i);

  sort(pAndIndex.begin(), pAndIndex.end());

  vector<double> intersections;

  for ( int i = 0 ; i < myD.N-1 ; ++i )
  {
    double value = (myD.c[i] - myD.c[i+1])/(myD.p[i+1] - myD.p[i]);
    intersections.push_back(value);
  }

  sort(intersections.begin(), intersections.end());

  pair<double, double> limits = make_pair(intersections[0]-1, intersections[myD.N-1]+1);

  return limits;
}

// Checks if problem is feasible
bool checkProblemFeasibility(Data &myD, pair<double, double> limits)
{
  double leftLambda = limits.first - 1;
  double rightLambda = limits.second + 1;

  Solution left(myD.N), right(myD.N);

  greedy(myD, left, leftLambda);
  greedy(myD, right, rightLambda);

  // checks if the signs are opposite. Also, if one of them is zero, then it also works
  if ( myD.subgradient(left) * myD.subgradient(right) <= 0 )
    return true;

  return false;
}

// Checks if the dual lines are intersecting at unique values
bool checkData(Data &myD)
{
  // find all the intersections
  // see that there are no repetitions
  
  vector<double> intersections;

  // Calculates all intersections
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    for ( int j = i + 1 ; j < myD.N ; ++j )
    {
      if ( myD.p[i] == myD.p[j] )
      {
        if ( myD.c[i] == myD.c[j] ) 
          return false; // Lines coincide
        else
          continue; // Lines are parallel
      }

      double value = (myD.c[i]-myD.c[j])/(myD.p[j]-myD.p[i]);
      intersections.push_back(value);
    }
  }

  // Sorts the values
  sort(intersections.begin(), intersections.end());
  
  // Checks for repeated ones
  for ( int i = 1 ; i < intersections.size() ; ++i )
  {
    if ( intersections[i] == intersections[i-1] )
      return false;
  }

  // Returns true when none found
  return true;
}

