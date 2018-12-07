#ifndef testSolver_hpp
#define testSolver_hpp

#include <iostream>
#include "DataStructures/Data.hpp"
#include "DataStructures/Solution.hpp"

Solution test(Data &myD);
Solution solveMosek(Data &myData);

void solveMosekRepeat(Data &myData, double &totalDuration, int numRepetitions);

Solution solveGreedyWithMosek(Data &myData);
Solution solveGreedy2WithMosek(Data &myData, double lambda);
#endif
