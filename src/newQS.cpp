#include <iostream>
#include <limits>
#include "newQS.hpp"
#include "commonFuncs.hpp"

// -----------------------------------------------------
// ---------------- New QuickSort ----------------------
// -----------------------------------------------------


/// -----------------------------------------------------------
/// ------------------ 1. CONSTRUCTORS ------------------------
/// -----------------------------------------------------------

RAPLBLSolver::RAPLBLSolver(Data &incomingData, double tolerance)
{
  //TODO: Check incoming data is ok, possibly with feasibility check
  this->dualSpace = SearchSpace();
  this->probData  = incomingData;
  this->tol = tolerance;

  this->v.resize(this->probData.N);
  for ( int i = 0 ; i < this->probData.N ; ++i )
    this->v[i] = i;
}

/// ---------------- 1. END - CONSTRUCTOR ---------------------
/// -----------------------------------------------------------


/// -----------------------------------------------------------
/// ------------------ 4. PUBLIC METHODS  ---------------------
/// -----------------------------------------------------------

void RAPLBLSolver::solve_StronglyPoly(Solution &candidate)
{
  qsQ.push(std::make_pair(0, this->probData.N-1));

  popAllqsQ(candidate);

  Solution leftSol(this->probData.N), rightSol(this->probData.N);
  double leftSubgrad = this->compare(this->dualSpace.getLeftLambda(), leftSol);
  double rightSubgrad = this->compare(this->dualSpace.getRightLambda(), rightSol);

  leftSubgrad = fabs(leftSubgrad);
  rightSubgrad = fabs(rightSubgrad);
  
  if ( leftSubgrad < rightSubgrad)
  {
    candidate = leftSol;
    return;
  }

  if ( rightSubgrad < leftSubgrad )
  {
    candidate = rightSol;
    return;
  }

  double leftObj = this->probData.F(leftSol);
  double rightObj = this->probData.F(rightSol);

  if ( leftObj < rightObj )
    candidate = leftSol;
  else
    candidate = rightSol;

  return;
}

/// ---------------- 4. END - PUBLIC METHODS ------------------
/// -----------------------------------------------------------

double RAPLBLSolver::compare( double lambda, Solution& candidate)
{
  // Get double solution (left and right solution at lambda)
  Solution left(this->probData.N), right(this->probData.N);
  greedyStrong( this->probData, lambda, left, right);

  // Get primal candidate
  convexSum(this->probData, left, right, candidate);
  
  return  this->probData.subgradient(candidate);
}

double getIntersection(Data &myD, int i, int j)
{
  if ( myD.p[i] == myD.p[j] )
    return numeric_limits<double>::infinity();

  return (myD.c[i] - myD.c[j])/(myD.p[j] - myD.p[i]);
}


int RAPLBLSolver::solveIntersection(double lambda)
{
  double leftLambda = this->dualSpace.getLeftLambda();
  double rightLambda = this->dualSpace.getRightLambda();

  if ( lambda <= leftLambda ) return 1;
  else if (lambda >= rightLambda ) return -1;

  Solution candidate(this->probData.N);
  double subgradient = compare(lambda, candidate);

  if ( subgradient > 0 )
  {
    this->dualSpace.setLeftLambda(lambda);
    return 1;
  }
  else if ( subgradient < 0 )
  {
    this->dualSpace.setRightLambda(lambda);
    return -1;
  }
  else
  {
    this->dualSpace.setLeftLambda(lambda);
    this->dualSpace.setRightLambda(lambda);
    return 0;
  }
}

void RAPLBLSolver::addToBatch(std::pair<int, int> qs, const int pivot)
{
  int low = qs.first;
  int high = qs.second;

  // get intersection
  for ( int i = low; i < high ; ++i )
  {
    double intersValue = getIntersection(this->probData, v[i], pivot);
    if ( intersValue > numeric_limits<double>::max() )
        continue;

    this->batch.push_back(intersValue);
  }
}

double RAPLBLSolver::solveBatch()
{
  int bot, mid, top;
  bot = 0;
  top = this->batch.size();
  mid = 0;

  while ( bot < top )
  {
    mid = bot + (top - bot)/2;

    std::nth_element(this->batch.begin(), this->batch.begin() + mid, this->batch.end());

    double subgradient = this->solveIntersection(this->batch[mid]); 

    if ( subgradient > this->tol )
      bot = mid+1;
    else if ( subgradient < -this->tol )
      top = mid;
    else
      break;
  }

  return this->batch[mid];

}



// We return +1 if pivot is greater at the dual optimal
// We return -1 if pivot is smaller at the dual optimal
int RAPLBLSolver::isPivotGrter(int pivot, int notPivot, double lambda)
{
  if ( this->probData.p[pivot] == this->probData.p[notPivot] )
    return this->probData.c[pivot] > this->probData.c[notPivot] ? 1 : -1 ;

  // Following the pattern: +1 if pivot has greater gradient (of the line c_i + \lambda p_i)
  //                        -1 if pivot has smaller gradient
  int pGradientGrter = this->probData.p[pivot] > this->probData.p[notPivot] ? 1 : -1;
  
  // By convention, negative subgradient implies that the dual optimal is to the left of intersection (+1)
  //                positive subgradient implies that the dual optimal is to the right of intersection (-1)
  // Thus, by multiplying them both together, we get all four possibilities.
  return pGradientGrter * this->solveIntersection( lambda);
}

int RAPLBLSolver::partition(const std::pair<int, int> partitionData)
{
  const int low = partitionData.first;
  const int high = partitionData.second;
  const int mid = low + (high - low)/2;
  const int pivot = this->v[mid];

  std::swap(this->v[mid], this->v[low]);

  int i = low; //low + 1;
  int j = high + 1; //high;

  while ( i < j )
  {
    double lambda;
    double aux;
    do
    {
      i++;
      aux = this->v[i];
      lambda = getIntersection(this->probData, aux, pivot);
    } while ( i < j && i < high && isPivotGrter(pivot, aux, lambda) < 0 ); 

    do
    {
      --j;
      aux = this->v[j];
      lambda = getIntersection(this->probData, aux, pivot);
    } while ( i <= j && j > low && isPivotGrter(pivot, aux, lambda) > 0);

    if ( i < j ) 
    {
      std::swap(this->v[i], this->v[j]);
    }
  }

  std::swap(this->v[j], this->v[low]);
  return j;
}

int RAPLBLSolver::popAllqsQ(Solution &candidate)
{
  while ( this->qsQ.empty() == false )
  {
    // pop queue
    std::pair<int, int> qsData = qsQ.front();
    qsQ.pop();
    int low = qsData.first;
    int high = qsData.second;

    if ( low >= high )
      continue;

    // choose pivot
    int mid = low + (high - low)/2;
    int pivot = this->v[mid];
    
    // create comparisons
    addToBatch(qsData, pivot);

    // add to partition queue
    partitionQ.push(qsData);
  }

  // Solve comparisons
  double lambdaPrime = solveBatch();

  // If optimal, return it
  double subgradient = compare(lambdaPrime, candidate);
  if ( subgradient == 0 )
    return -1;

  if ( this->partitionQ.empty() == true ) 
    return 0;

  return popAllpartitionQ(candidate);
}

int RAPLBLSolver::popAllpartitionQ(Solution &candidate)
{
  while ( this->partitionQ.empty() == false )
  {
    // pop queue
    std::pair<int, int> partitionData = this->partitionQ.front();
    this->partitionQ.pop();

    int low = partitionData.first;
    int high = partitionData.second;

    // partition data
    int pivotIndex  = this->partition(partitionData);

    if ( pivotIndex > high )
      continue;

    // add to qsQ
    std::pair<int, int> leftSort = std::make_pair(low, pivotIndex-1);
    std::pair<int, int> rightSort = std::make_pair(pivotIndex+1, high);

    this->qsQ.push(leftSort);
    this->qsQ.push(rightSort);
  }

  return popAllqsQ(candidate);
}

void startQuickSort(int low, int high, Solution &candidate)
{
  //qsData.push(make_std::pair<low, high>);

  //popQSQueue(this->probData, candidate, qsData, partData);
}


