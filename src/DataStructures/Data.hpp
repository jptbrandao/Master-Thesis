#ifndef Data_hpp
#define Data_hpp

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "Solution.hpp"

class Data
{
  public:
    int N;            ///< Number of variables
 
    vector<double> c; ///< Linear weights
    double B1;        ///< Resource
    double B2;        ///< Linear Constraint
    vector<double> p; ///< Linear Constraint 

    double lambda;    ///< Lagrangean Dual

    int M;                    ///< Number of Nested Constraint
    vector<bool> nestCst;     ///< Boolean vector to determine which nested constraints are active
    //vector<int> indxToNest;   ///< Array that informs which nested constraint a variable index belongs to
    vector<double> upperCst;  ///< Upper bounds of nested constraint

    //vector<int> nestedLargestIdx;

    // Constructor with random initialization
    Data( int numVar, unsigned int seed);

    // Basic test case
    Data();

    Data(int numVar, int numCst);

    // Constructor that reads data from a file with several problems instances, and can determine which instance to read by the position in the file.
    Data(string filepath, ios::pos_type &position);

    Data(string filepath);

    // Randomly reassign values
    void restart();

    // Prints to console
    void printData();

    // Evaluates the objective function of the optimization instance
    double F(Solution &x);

    // Exports problem's instance to a text file
    void exportData(std::string filepath);

    // Exports problem's instance to a text file with additional notes
    void exportData(string filepath, string notes);


    // Determines if a solution is feasible
    bool isFeasible(Solution &x);

    double subgradient(Solution &x);

    // Destructor
    ~Data();
};

#endif
