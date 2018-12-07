//
//  Solution.cpp
//  tryone
//
//  Created by Joao Pedro Brandao on 1/10/17.
//  Copyright © 2017 Joao Pedro Brandao. All rights reserved.
//

#include "Solution.hpp"

Solution::Solution()
{
  numVar = 0;
  timeTaken = 0;
}

Solution::Solution(int N, double* opt, double time)
{
    numVar = N;
    timeTaken = time;
    
    for ( int i = 0 ; i < N ; ++i )
        sol.push_back(opt[i]);
}

Solution::Solution(int N)
{
    sol.resize(N);
    numVar = N;
}

Solution::Solution(Solution const &x)
{
    numVar = x.numVar;
    timeTaken = x.timeTaken ;
    
    for ( int i = 0 ; i < numVar ; ++i )
        sol.push_back(x.sol[i]);
}

Solution::Solution(Solution &x)
{
  numVar = x.numVar;
  timeTaken = x.timeTaken ;

  for ( int i = 0 ; i < numVar ; ++i )
    sol.push_back(x.sol[i]);
}

void Solution::printSol()
{
    cout << "Solution: " << numVar << endl;
    
    for ( int i = 0 ; i < numVar ; ++i )
        cout << i << " " << sol[i] << endl;
}


double& Solution::operator[](int x)
{
    if ( x < 0 || x > sol.size() )
    {
        cout << "Invalid Index" << endl;
        exit(1);
    }
    
    return sol[x];
}

void Solution::setIndexZero(int index)
{
    sol[index] = 0;
}

void Solution::setDual(double dual)
{
    dualVar = dual;
}

double Solution::dual()
{
    return dualVar;
}

void Solution::setTime(double time)
{
    timeTaken = time;
}

double Solution::Time()
{
    return timeTaken;
}

bool Solution::isEmpty()
{
    return sol.empty();
}
