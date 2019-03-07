#ifndef SearchSpace_hpp
#define SearchSpace_hpp


class SearchSpace
{
  public:
    // Constructors
    SearchSpace();

    // Member functions
    
    // Get/Sets
    double getLeftLambda();
    double getRightLambda();
    void setLeftLambda(double lambda);
    void setRightLambda(double lambda);
    void resetLambdas();

  private:
    double leftLambda, rightLambda; 
};


#endif
