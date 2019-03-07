#include "SearchSpace.hpp"
#include <limits>


/// -----------------------------------------------------------
/// ------------------ 1. CONSTRUCTORS ------------------------
/// -----------------------------------------------------------

SearchSpace::SearchSpace()
{
  this->resetLambdas();
}
 

/// ---------------- 1. END - CONSTRUCTOR ---------------------
/// -----------------------------------------------------------


/// -----------------------------------------------------------
/// ------------------ 4. PUBLIC METHODS  ---------------------
/// -----------------------------------------------------------


void SearchSpace::resetLambdas()
{
  this->leftLambda = -std::numeric_limits<double>::max();
  this->rightLambda = std::numeric_limits<double>::max();
}

double SearchSpace::getLeftLambda()
{
  return this->leftLambda;
}

double SearchSpace::getRightLambda()
{
  return this->rightLambda;
}

void SearchSpace::setLeftLambda(double lambda)
{
  this->leftLambda = lambda;
}

void SearchSpace::setRightLambda(double lambda)
{
  this->rightLambda = lambda;
}
