#ifndef newQS_hpp
#define newQS_hpp

#include "DataStructures/Data.hpp"
#include "DataStructures/Solution.hpp"
#include <queue>
#include "SearchSpace.hpp"
#include <utility> 
#include <vector>

class RAPLBLSolver
{
  public:
    // Variables
    
    // Constructors
    RAPLBLSolver(Data &myD, double tolerance);

    // Member functions
    void solve_StronglyPoly(Solution &candidate);
    
  private:
    SearchSpace dualSpace;
    std::vector<double> v;
    std::queue<std::pair<int, int>> qsQ;
    std::queue<std::pair<int, int>> partitionQ;
    std::vector<double> batch;
    Data probData;
    double tol;

    // QuickSort parametric search
    int popAllqsQ(Solution &candidate);
    int popAllpartitionQ(Solution &candidate);
    void addToBatch(std::pair<int, int> qsData, const int pivot);
    double solveBatch();
    int solveIntersection(double lambda);
    int partition(const std::pair<int, int> partitionData);
    int isPivotGrter(int pivot, int notPivor, double lambda);
    double compare(double lambda, Solution &candidate);
};


#endif
