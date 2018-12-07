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

void bruteForceStrong(Data &myD, Solution &sol);

int inversion(vector<double> &v, int low, int high);

double getLeftLambda();
double getRightLambda();
void setLeftLambda(double lambda);
void setRightLambda(double lambda);
void resetLambdas();
#endif
