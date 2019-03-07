#include <limits>
#include "commonFuncs.hpp"


double toler = 0.00000001; // added 2 zeros 19/1/19

unsigned long long mergeAndCount(vector<double> &v, int low, int mid, int high)
{
  unsigned long long count = 0;

  vector<double> L, R;

  for ( int i = low ; i <=mid; ++i ) 
    L.push_back(v[i]);
  for ( int i = mid+1 ; i <= high; ++i )
    R.push_back(v[i]);

  L.push_back(numeric_limits<double>::max());
  R.push_back(numeric_limits<double>::max());

  int i = 0 ;
  int j = 0 ;

  for ( int k = low ; k <= high ; ++k ) 
  {
    if ( L[i] <= R[j] )
    {
      v[k] = L[i];
      ++i;
    }
    else
    {
      if ( L[i] != L.back() && R[j] != R.back() )
        count += mid - low + 1 - i ;

      v[k] = R[j];
      ++j;
    }
  }

  return count;
}

unsigned long long inversion(vector<double> &v, int low, int high)
{
  unsigned long long count = 0;
  if ( low < high )
  {
    int mid = (low+high)/2;
    count = inversion(v, low, mid);
    count += inversion(v, mid+1, high);
    count += mergeAndCount(v, low, mid, high);
  }

  return count;
}


// Performs a convex sum on solutions left and right, and returns it on solution cand
void convexSum( Data &myD, Solution &left, Solution &right, Solution &cand)
{
  double sumLeft = 0;
  double sumRight = 0;

  for ( int i = 0 ; i < myD.N ; ++i )
  {
    sumLeft += myD.p[i]*left[i];
    sumRight += myD.p[i]*right[i];
  }
  
  if ( sumLeft == sumRight )
  {
    cand = left ;
    return;
  }

  double deltaLeft = fabs(sumLeft - myD.B2); 
  double deltaRight = fabs(myD.B2 - sumRight); 

  double mu = deltaRight/(deltaRight + deltaLeft);

  for ( int i = 0 ; i < myD.N ; ++i )
    cand[i] = mu * left[i] + (1-mu) * right[i];

  return;
}

// TODO: add dual variable
// Performs a greedy solution on 
void greedy( Data &myD, Solution &sol, double lambda)
{
  int minIdx = myD.N-1;
  int cstIdx = myD.M-1;

  double resource = myD.upperCst[cstIdx] - myD.upperCst[cstIdx-1];
  cstIdx--;

  double c_Min = myD.c[minIdx] + lambda * myD.p[minIdx];

  sol[myD.N-1] = 0;
  for ( int i = myD.N-2 ; i >= 0 ; --i )
  {
    sol[i] = 0;

    if ( myD.nestCst[i] == true ) 
    {
      sol[minIdx] += resource;

      if ( cstIdx > 0 )
      {
        // Update resource
        resource = myD.upperCst[cstIdx] - myD.upperCst[cstIdx-1];
        cstIdx--;
      }
    }
    // Calculate current y value
    double c_i = myD.c[i] + lambda * myD.p[i]; 
    
    // Change minimum if c is lower than current minimum
    if (  c_i - c_Min  < 0 )
    {
      minIdx = i;
      c_Min = c_i;
    }
  }

  sol[minIdx] += myD.upperCst[0];

  return;
}

void greedyStrong( Data &myD, double lambda, Solution &left, Solution &right)
{
  int minIdxLeft = myD.N-1;       // minimum left Index
  int minIdxRight = myD.N-1;      // minimum right Index
  int cstIdx = myD.M-1;           // Constraint Index
  
  double resource = myD.upperCst[cstIdx] - myD.upperCst[cstIdx-1];  
  cstIdx--;
  //cout << "B_1 - b_" << cstIdx+1 << " = " << resource << endl;

  double y_Min = myD.c[myD.N-1] + myD.p[myD.N-1] * lambda;  // y = c + \lambda p

  left[minIdxLeft] = 0;
  right[minIdxRight] = 0;

  for ( int i = myD.N-2 ; i >= 0 ; --i )
  {
    // Initialize solutions with zero
    left[i] = 0;
    right[i] = 0;
   
    // Update values if entering a nested constraint
    if ( myD.nestCst[i] == true )
    {
      left[minIdxLeft] += resource;
      right[minIdxRight] += resource;

      if ( cstIdx > 0 )
      {
        // Update resource
        resource = myD.upperCst[cstIdx] - myD.upperCst[cstIdx-1];
        cstIdx--;
      }
    }

    // Calculate current y value
    double y_i = myD.c[i] + myD.p[i] * lambda;

    // Change minimum if y is lower than current minimum
    if ( fabs( y_i - y_Min ) < toler )
    {
      // Prioritize most positive gradient for left solution
      // Prioritize most negative gradient for right solution
      if ( myD.p[i] > myD.p[minIdxLeft] )
      {
        minIdxLeft = i;
      }
      else if ( myD.p[i] < myD.p[minIdxRight] )
      {
        minIdxRight = i;
      }

    }
    else if ( y_i < y_Min )
    {
      minIdxLeft = i;
      minIdxRight = i;

      y_Min = myD.c[i] + myD.p[i] * lambda;
    }
  }

  left[minIdxLeft] += myD.upperCst[0];
  right[minIdxRight] += myD.upperCst[0];

  return;
}


// Checks if problem is feasible
bool checkProblemFeasibility(Data &myD, pair<double, double> limits)
{
  double leftLambda = limits.first - 1;
  double rightLambda = limits.second + 1;

  Solution left(myD.N), right(myD.N);

  greedy(myD, left, leftLambda);
  greedy(myD, right, rightLambda);

  // checks if the signs are opposite. Also, if one of them is zero, then it also works
  if ( myD.subgradient(left) * myD.subgradient(right) <= 0 )
    return true;

  return false;
}

// Checks if the dual lines are intersecting at unique values
bool checkData(Data &myD)
{
  // find all the intersections
  // see that there are no repetitions
  
  vector<double> intersections;

  // Calculates all intersections
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    for ( int j = i + 1 ; j < myD.N ; ++j )
    {
      if ( myD.p[i] == myD.p[j] )
      {
        if ( myD.c[i] == myD.c[j] ) 
          return false; // Lines coincide
        else
          continue; // Lines are parallel
      }

      double value = (myD.c[i]-myD.c[j])/(myD.p[j]-myD.p[i]);
      intersections.push_back(value);
    }
  }

  // Sorts the values
  sort(intersections.begin(), intersections.end());
  
  // Checks for repeated ones
  for ( int i = 1 ; i < intersections.size() ; ++i )
  {
    if ( intersections[i] == intersections[i-1] )
      return false;
  }

  // Returns true when none found
  return true;
}


pair<double, double> findBigLambdaInterval(Data &myD)
{
  vector<pair<double, int> > pAndIndex(myD.N);

  for ( int i = 0 ; i < myD.N ; ++i )
    pAndIndex[i] = make_pair(myD.p[i], i);

  sort(pAndIndex.begin(), pAndIndex.end());

  vector<double> intersections;

  for ( int i = 0 ; i < myD.N-1 ; ++i )
  {
    double value = (myD.c[i] - myD.c[i+1])/(myD.p[i+1] - myD.p[i]);
    intersections.push_back(value);
  }

  sort(intersections.begin(), intersections.end());

  pair<double, double> limits = make_pair(intersections[0]-1, intersections[myD.N-1]+1);

  return limits;
}
