#ifndef weaklypoly_hpp
#define weaklypoly_hpp

#include <vector>
#include "DataStructures/Data.hpp"
#include "DataStructures/Solution.hpp"

using namespace std;

void weaklyPoly(Data &myD, double lambdaMin, double lambdaMax, Solution &x);

bool checkProblemFeasibility(Data &myD, pair<double, double> limits);

bool checkData(Data &myD);

pair<double, double> findBigLambdaInterval(Data &myD);
#endif
