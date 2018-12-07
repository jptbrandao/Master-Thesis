#include <iostream>
#include <limits>
#include "quicksort.hpp"
#include "commonFuncs.hpp"


double tol = 0.001;

void printVector( vector<double> &v, int start, int end )
{
  cout << "( " << v[start];
  for ( int i = start+1 ; i < end; ++i )
   cout << ", " << v[i];

  cout << " )" << endl;
  return;
}


double compare(Data &myD, double lambda, Solution& candidate)
{
  // Get double solution (left and right solution at lambda)
  Solution left(myD.N), right(myD.N);
  greedyStrong( myD, lambda, left, right);

  // Get primal candidate
  convexSum(myD, left, right, candidate);
  
  return  myD.subgradient(candidate);
}

double getIntersection(Data &myD, int i, int j)
{
  if ( myD.p[i] == myD.p[j] )
    return numeric_limits<double>::infinity();

  return (myD.c[i] - myD.c[j])/(myD.p[j] - myD.p[i]);
}



int solveIntersection(Data &myD, double lambda)
{
  double leftLambda = getLeftLambda();
  double rightLambda = getRightLambda();

  if ( lambda <= leftLambda ) return 1;
  else if (lambda >= rightLambda ) return -1;

  Solution candidate(myD.N);
  double subgradient = compare(myD, lambda, candidate);

  if ( subgradient > 0 )
  {
    setLeftLambda(lambda);
    return 1;
  }
  else if ( subgradient < 0 )
  {
    setRightLambda(lambda);
    return -1;
  }
  else
  {
    setLeftLambda(lambda);
    setRightLambda(lambda);
    return 0;
  }
}

double batchBinarySearch(Data &myD, vector<double> &v, const int low, const int high, const int pivot)
{
  vector<double> intersections;

  for ( int i = low; i < high ; ++i )
  {
    double intersValue = getIntersection(myD, v[i], pivot);
    if ( intersValue > numeric_limits<double>::max() )
        continue;

    intersections.push_back(intersValue);
  }

  if ( intersections.size() == 0 )
  {
    return getRightLambda() + 100; // Just returns a random value that is not dual optimal
  }

  sort(intersections.begin(), intersections.end());
  
  int bot, mid, top;
  bot = 0;
  top = intersections.size();
  mid = 0;

  while ( bot < top )
  {
    mid = bot + (top - bot)/2;

    double subgradient = solveIntersection(myD, intersections[mid]);    

    if ( subgradient > tol )
      bot = mid+1;
    else if ( subgradient < -tol )
      top = mid;
    else
      break;
  }

  return intersections[mid];
}


// We return +1 if pivot is greater at the dual optimal
// We return -1 if pivot is smaller at the dual optimal
int isPivotGrter(Data &myD, int pivot, int notPivot, double lambda)
{
  if ( myD.p[pivot] == myD.p[notPivot] )
    return myD.c[pivot] > myD.c[notPivot] ? 1 : -1 ;

  // Following the pattern: +1 if pivot has greater gradient (of the line c_i + \lambda p_i)
  //                        -1 if pivot has smaller gradient
  int pGradientGrter = myD.p[pivot] > myD.p[notPivot] ? 1 : -1;
  
  // By convention, negative subgradient implies that the dual optimal is to the left of intersection (+1)
  //                positive subgradient implies that the dual optimal is to the right of intersection (-1)
  // Thus, by multiplying them both together, we get all four possibilities.
  return pGradientGrter * solveIntersection(myD, lambda);
}

int partition(Data &myD, vector<double> &v, const int low, const int high, Solution &candidate)
{
  const int mid = low + (high - low)/2;
  const int pivot = v[mid];

  std::swap(v[mid], v[low]);

  int i = low; //low + 1;
  int j = high + 1; //high;

  
  //double lambdaPrime = batchBinarySearch(myD, v, low + 1, high, pivot);
  double lambdaPrime = batchBinarySearch(myD, v, i, j, pivot);

  double subgradient = compare(myD, lambdaPrime, candidate);
  if ( subgradient == 0 )
    return -1;

  while ( i < j )
  {
    double lambda;
    double aux;
    do
    {
      i++;
      aux = v[i];
      lambda = getIntersection(myD, aux, pivot);
    } while ( i < j && i < high && isPivotGrter(myD, pivot, aux, lambda) < 0 ); 

    do
    {
      --j;
      aux = v[j];
      lambda = getIntersection(myD, aux, pivot);
    } while ( i <= j && j > low && isPivotGrter(myD, pivot, aux, lambda) > 0);

    if ( i < j ) 
    {
      std::swap(v[i], v[j]);
    }
  }

  std::swap(v[j], v[low]);
  return j;
}

void chooseSolution(Data &myD, Solution &candidate)
{
  Solution leftSol(myD.N), rightSol(myD.N);
  double leftSubGrad = compare(myD, getLeftLambda(), leftSol);
  double rightSubGrad = compare(myD, getRightLambda(), rightSol);

  leftSubGrad = fabs(leftSubGrad);
  rightSubGrad = fabs(rightSubGrad);

  if ( leftSubGrad < rightSubGrad )
  {
    candidate = leftSol;
    return;
  }

  if ( leftSubGrad > rightSubGrad )
  {
    candidate = rightSol;
    return;
  }

  double leftObj = myD.F(leftSol);
  double rightObj = myD.F(rightSol);

  if ( leftObj < rightObj )
  {
    candidate = leftSol;
    return;
  }

  candidate = rightSol;
  return;
}

// I'm using the convention that if the return value is negative, the optimal has been found
int customSort(Data &myD, vector<double> &v, int low, int high, Solution &candidate)
{
  if ( low >= high )
    return 0;
  
  bool isLastRecursion = (high - low == v.size() -1 );
  int pivot = partition(myD, v, low, high, candidate); 

  if ( pivot < 0 )
    return pivot;

  // It stop the unnecessary recursion
  if ( pivot > high )
    return pivot;

  if ( customSort(myD, v, low, pivot-1, candidate) < 0 && (isLastRecursion == false) )
    return -1;
  
  if ( customSort(myD, v, pivot+1, high, candidate) < 0 && (isLastRecursion == false) )
    return -1;

  if ( isLastRecursion == true )
  {
    chooseSolution(myD, candidate);
  }

  return pivot;
}
