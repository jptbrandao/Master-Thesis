//
//  Solution.hpp
//  tryone
//
//  Created by Joao Pedro Brandao on 1/10/17.
//  Copyright © 2017 Joao Pedro Brandao. All rights reserved.
//

#ifndef Solution_hpp
#define Solution_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;

class Solution
{
private:
    int numVar;
    double timeTaken;
    vector<double> sol;
    double dualVar;
    
public:
    
    Solution(int N, double* opt, double time);
    Solution();
    Solution(int N);
    
    Solution(Solution const &x);
    Solution(Solution &x);
    
    
    void printSol();
    
    
    void setIndexZero(int index);
    void setDual(double dual);
    double dual();
    
    double& operator[](int x);
    
    void setTime(double time);
    double Time();
    
    bool isEmpty();
};

#endif /* Solution_hpp */
