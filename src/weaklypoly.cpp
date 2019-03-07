#include <iostream>
#include <limits>
#include "weaklypoly.hpp"
#include "commonFuncs.hpp"

double tolz = 0.00001;



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
  unsigned long countLeft = countInversions(myD, lambdaMin);
  unsigned long countRight = countInversions(myD, lambdaMax);

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
