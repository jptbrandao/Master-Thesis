#ifndef commonFuncs_hpp
#define commonFuncs_hpp

#include <vector>
#include "DataStructures/Data.hpp"
#include "DataStructures/Solution.hpp"

#define MaxFloat numeric_limits<double>::max()

using namespace std;

void convexSum( Data &myD, Solution &left, Solution &right, Solution &cand);

void greedy( Data &myD, Solution &sol, double lambda);

void greedyStrong( Data &myD, double lambda, Solution &left, Solution &right);

unsigned long long inversion(vector<double> &v, int low, int high);

bool checkProblemFeasibility(Data &myD, pair<double, double> limits);

bool checkData(Data &myD);

pair<double, double> findBigLambdaInterval(Data &myD);


#endif
