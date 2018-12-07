#include <limits>
#include "commonFuncs.hpp"


static double leftLambda = -numeric_limits<double>::max();
static double rightLambda = numeric_limits<double>::max();
double toler = 0.000001;


double getLeftLambda()
{
  return leftLambda;
}

double getRightLambda()
{
  return rightLambda;
}

void setLeftLambda(double lambda)
{
  leftLambda = lambda;
}

void setRightLambda(double lambda)
{
  rightLambda = lambda;
}

void resetLambdas()
{
  leftLambda = -numeric_limits<double>::max();
  rightLambda = numeric_limits<double>::max();
}

int mergeAndCount(vector<double> &v, int low, int mid, int high)
{
  int count = 0;

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

int inversion(vector<double> &v, int low, int high)
{
  int count = 0;
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



vector<double> generateIntersections(Data &myD)
{
  vector<double> intersections;
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    for ( int j = i+1 ; j < myD.N ; ++j )
    {
      if ( myD.p[i] == myD.p[j] ) continue;
      if ( myD.c[i] == myD.c[j] ) continue;

      double lambda = (myD.c[i] - myD.c[j])/(myD.p[j] - myD.p[i]);

      intersections.push_back(lambda);
    }
  }

  return intersections;
}

void getPrimalCandidate(Data &myD, double lambda, Solution &sol)
{
  Solution left(myD.N), right(myD.N);
  greedyStrong(myD, lambda, left, right);

  convexSum(myD, left, right, sol);
  return;
}

double medianOf5(vector<double> &v, int left, int right )
{
  // TODO: check iterator indices
  sort(v.begin() + left, v.begin() + right + 1);

  int mid = left + (right-left)/2;
  return v[mid];
}

double medianOfMedians(vector<double> &v, int left, int right);

double medianOfMedians(vector<double> &v, int left, int right)
{
  if ( right - left < 5 )
    return medianOf5(v, left, right);

  for ( int i = left ; i <= right ; i += 5 )
  {
    int subRight = i + 4;
    if ( subRight > right )
      subRight = right;

    double median5 = medianOf5(v, left, subRight);
    int idx = left + (i-left)/5;
    std::swap(v[median5], v[idx]);
  }
  
}

void bruteForceStrong( Data &myD, Solution &sol )
{
  // generate breakpoints
  vector<double> intersections = generateIntersections(myD);

  sort(intersections.begin(), intersections.end());
  intersections.erase( unique (intersections.begin(), intersections.end()), intersections.end() );
  
  // binary search
  int bot, mid, top;
  bot = 0;
  top = intersections.size();
  mid = 0;

  while ( bot < top )
  {
    mid = bot + (top-bot)/2;

    getPrimalCandidate(myD, intersections[mid], sol);
    double subgradient = myD.subgradient(sol);

    if ( subgradient > toler )
      bot = mid + 1;
    else if ( subgradient < -toler )
      top = mid;
    else
      break;
  }
}
