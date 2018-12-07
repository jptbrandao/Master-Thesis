#ifndef quicksort_hpp
#define quicksort_hpp

#include <vector>
#include "DataStructures/Data.hpp"
#include "DataStructures/Solution.hpp"

using namespace std;

int customSort(Data &myD, vector<double> &v, int low, int high, Solution &candidate);
void printVector(vector<double> &v, int start, int end);

void greedyStrong(Data &myD, double lambda, Solution &left, Solution &right);


void simpleGreedy(Data &myD, Solution &sol);

void convexSum( Data &myD, Solution &left, Solution &right, Solution &cand);

#endif

