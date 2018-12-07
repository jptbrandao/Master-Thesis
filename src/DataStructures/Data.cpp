#include "Data.hpp"
#include <string>

void normalizeUpperCst(Data &myD);

/// -----------------------------------------------------------
/// ------------------ 1. CONSTRUCTORS ------------------------
/// -----------------------------------------------------------

Data::Data(int numVar, int numCst)
{
  N = numVar;
  c.resize(N);
  p.resize(N);
  nestCst.resize(N);

  M = numCst;
  upperCst.resize(M);
  lambda = 0;
}

Data::Data(string filepath)
{
  ifstream file(filepath, ios::in);
  
  string line;

  getline(file, line);  // Program Name
  getline(file, line);  // Problem Type

  getline(file, line, ':'); // Number of Variables
  file >> N;
  M = N;
  //cout << "numVar: " << N << endl;
  getline(file, line);  // Empty line

  getline(file, line, ':'); // Linear Term
  file >> B2;
  //cout << "B2: " << B2 << endl;
  getline(file, line); // Empty line

  getline(file, line);  // Bounds Section
  
  nestCst.resize(N);
  upperCst.resize(N);
  c.resize(N);
  p.resize(N);
  for ( int i = 0 ; i < N ; ++i )
  {
    int index;
    file >> index >> upperCst[i] >> p[i];
    nestCst[i] = true;
  }

  B1 = upperCst[N-1];

  getline(file, line);  // Empty line
  getline(file, line);  // Function Profile

  for ( int i = 0 ; i < N ; ++i )
  {
    int index;
    file >> index >> c[i];
  } 
  
    normalizeUpperCst(*this);
}

/**
 *   @brief  Constructor that reads data from a file with several problems instances
 *
 *   @param  filepath Filepath to the file with the instances of problems
 *   @param  position Position of ifstream in the file.
 */
Data::Data(string filepath, ios::pos_type &position)
{
    ifstream arq;
    arq.open(filepath);
    
    // Sets the ifstream's position to the proper problem
    arq.seekg(position);
    
    if (arq.is_open())
    {
        // N, M
        arq >> N >> M;
        
        // Linear term
        c.resize(N);
        for ( int i = 0 ; i < N ; ++i )
            arq >> c[i] ;
        
        // Resource
        arq >> B1 >> B2;
        
        // Linear constraints
        p.resize(N);
        for ( int i = 0 ; i < N ; ++i )
            arq >> p[i] ;
        
        // Upper nested constraints
        nestCst.resize(N);
        for ( int i = 0 ; i < N ; ++i )
        {
            int val;
            arq >> val ;
            
            if ( val == 1 )
                nestCst[i] = true;
            else
                nestCst[i] = false;
        }
        
        upperCst.resize(M);
        for ( int i = 0 ; i < M ; ++i )
            arq >> upperCst[i] ;
    }
    
    //indxToNest.resize(N);
    //int nestCount = 0;
    //for ( int i = 0 ; i < N ; ++i )
    //{
    //    indxToNest[i] = nestCount;
    //    
    //    if ( nestCst[i] == true )
    //        nestCount++;
   // }
    
    // updates position in file
    position = arq.tellg();
    arq.close();

    normalizeUpperCst(*this);
}

bool hasSeedBeenSet = false;
// Randomly assigned instance
// If seed is positive, set seed, otherwise, use seed
Data::Data(int numVar, unsigned int seed)
{
    N = numVar ;
    M = numVar;
    c.resize(N);
    p.resize(N);
    upperCst.resize(M);
    nestCst.resize(M);
     
    if ( !hasSeedBeenSet )
    {
      if ( seed == 0 )
      {
        unsigned int newSeed = time(NULL);
        cout << "Seed: " << seed << endl;
        srand(newSeed);
      }
      else
        srand(seed);

      hasSeedBeenSet = true;
    }
    
    restart();
    
    //indxToNest.resize(N);
    //int nestCount = 0;
    
    //for ( int i = 0 ; i < N ; ++i )
    //{
     //   indxToNest[i] = nestCount;
      //  if ( nestCst[i] == true )
       //     nestCount++;
    //}
}

Data::Data()
{
    N = 3 ;
    c.resize(N); //= new double[this->N] ;
    p.resize(N); //= new double[this->N] ;
    
    M = 3;
    upperCst.resize(M);
    lambda = 0;
    nestCst.resize(M);
    
    c[0] = -5 ;
    c[1] = 7;
    c[2] = 30;
    
    p[0] = 1;
    p[1] = -1;
    p[2] = 2;
   
    upperCst[0] = 2;
    upperCst[1] = 3;
    upperCst[2] = 5;
    
    B1 = upperCst[2];
    B2 = 5;
    
    nestCst[0] = true;
    nestCst[1] = true;
    nestCst[2] = true;
    //B = 50;
}

Data::~Data()
{
}


/// ---------------- 1. END - CONSTRUCTOR ---------------------
/// -----------------------------------------------------------


/// -----------------------------------------------------------
/// ------------------ 2. I/O METHODS -------------------------
/// -----------------------------------------------------------

/**
 *   @brief  Exports data to text file
 *
 *   @param  filepath Filepath to the text file
 */
void Data::exportData(string filepath)
{
    fstream arq;
    arq.open(filepath, std::ios_base::app);
    
    // Number of variables and number of constraints
    arq << N << " " << M << endl;
    
    // Linear term
    for ( int i = 0 ; i < N ; ++i )
        arq << c[i] << " " ;
    arq << endl;
    
    // Resource B1 and B2
    arq << B1 << " " << B2 << endl ;
    
    // Linear constraints
    for ( int i = 0 ; i < N ; ++i )
        arq << p[i] << " " ;
    arq << endl;
    
    // Upper nested constraints
    for ( int i = 0 ; i < N ; ++i )
    {
        if ( nestCst[i] == true )
            arq << "1" << " " ;
        else
            arq << "0" << " " ;
    }
    arq << endl;
    
    for ( int i = 0 ; i < M ; ++i )
        arq << upperCst[i] << " ";
    arq << endl;
    
    arq << endl;
    arq.close();
}

void Data::printData()
{
    cout << endl;
    cout << "\t ------ Instance Data -----" << endl;
    cout << "\t --- Objective Function --- " << endl;
    cout << "Number of Variables:\t " << N << endl;
    cout << "Number of Constraints:\t " << M << endl;
    
    cout << "Linear Terms: " << endl;
    cout << c[0];
    for ( int i = 1 ; i < N ; ++i )
        cout << ", " << c[i] ;
    cout << endl;
    
    cout << endl << "\t --- Constraints ---" << endl;
    cout << "Resource:\t "<< B1 << endl;
    cout << "HyperPlane:\t " << B2 << endl;
    
    cout << "Hyperplane Terms: " << endl;
    cout << p[0];
    for ( int i = 1 ; i < N ; ++i )
        cout << ", " << p[i] ;
    cout << endl;
    
    cout << "Nested Constraints: " << endl;
    cout << upperCst[0];
    for ( int i = 1 ; i < M ; ++i )
        cout << ", " << upperCst[i] ;
    cout << endl;
    cout << endl;
}

void Data::exportData(string filepath, string notes)
{
    fstream arq ;
    arq.open(filepath, std::ios_base::app);

    arq << " ------------ \t -------------- \t --------------- " << endl;

    arq << notes << endl;

    // Number of variables and number of constraints
    arq << N << " " << M << endl;

    // Linear term
    for ( int i = 0 ; i < N ; ++i )
        arq << c[i] << " " ;
    arq << endl;

    // Resource
    arq << B1 << " " << B2 << endl ;

    // Hyperplane term
    for ( int i = 0 ; i < N ; ++i )
        arq << p[i] << " " ;
    arq << endl;

    // Upper nested constraints
    for ( int i = 0 ; i < M ; ++i )
        arq << upperCst[i] << " ";
    arq << endl;

    arq << endl << endl;
    arq.close();
}


/// ---------------- 2. END - I/O METHODS ---------------------
/// -----------------------------------------------------------



/// -----------------------------------------------------------
/// ------------------ 3. INSTANCE METHODS --------------------
/// -----------------------------------------------------------

void Data::restart()
{
    this->lambda = 0;
    double upperLim = 0;
    
    for ( int i = 0 ; i < M ; i++ )
    {
        upperCst[i] = (double)(rand()%(N*2)) + upperLim;
        nestCst[i] = true;
        
        upperLim = upperCst[i];
    }
    
    B1 = upperCst[N-1] ;
    B2 = rand()%100 + 1 + B1;
    
    for ( int i = 0 ; i < N ; i++ )
    {
        c[i] = (double)(rand()%10);
        c[i] = (rand()%2 == 0 ? c[i] : -c[i]);
        
        p[i] = (double)(rand()%10)+1;
        p[i] = (rand()%2 == 0 ? p[i] : -p[i]);
        
//        nestCst[i] = rand()%2 == 0;
    }
}


/// ---------------- 3. END - INSTANCE METHODS ----------------
/// -----------------------------------------------------------

/// -----------------------------------------------------------
/// ------------------ 4. PUBLIC METHODS  ---------------------
/// -----------------------------------------------------------

/**
 *   @brief  Checks if solution satisfies the problem's constraints
 *
 *   @param  x Problem's solution
 *   @return True if x satisfies constraints, false otherwise
 */
bool Data::isFeasible(Solution &x)
{
    double tol = 0.001;     // Tolerance defined, as floating point numbers aren't always exact
    
    // Checks if solution is empty
    if ( x.isEmpty() == true )
    {
        cout << "Infeasible" << endl; 
    }
    
    double sumR1 = 0;       // sumR1 - resource1
    double sumR2 = 0;       // sumR2 - resource2
    for ( int i = 0 ; i < N ; ++i )
    {
        sumR1 += x[i];
        sumR2 += x[i]*p[i];
    }
    
    // Checks B1 and B2 constraints
    double resB1 = fabs(B1 - sumR1) ;
    double resB2 = fabs(B2 - sumR2) ;
    
    //cout << "Checking B1 and B2 " << endl;
    if ( resB1 > tol )  return false; //0.00000000000011368683772161603
    if ( resB2 > tol )  return false; // 0
    
    // Checks nonnegativity
    //cout << "Checking nonnegativity" << endl;
    for ( int i = 0 ; i < N ; ++i )
        if ( x[i] < 0 )
            return false;
    
    // Checks nested constraints
    //cout << "Checking nested " << endl;
    sumR1 = 0;
    for ( int i = 0 ; i < N ; ++i )
    {
        sumR1 += x[i];
        if ( nestCst[i] && sumR1 > upperCst[i] + tol)
            return false;
    }
    
    // If all constraints are satisfied, returns true
    return true;
}


double Data::subgradient(Solution &x)
{
  double subgradient = 0;

  for ( int i = 0 ; i < N ; ++i )
  {
    subgradient += p[i] * x[i];
  }

  subgradient -= B2;
  return subgradient;
}


/**
 *   @brief  Evaluates objective function
 *
 *   @param  x Problem's solution
 *   @return Objective function's value
 */
double Data::F(Solution &x)
{
    double funcSum = 0;
    
    for (int i = 0 ; i < N ; ++i )
        funcSum += c[i]*x[i] ;
    
    return funcSum;
}

/// ---------------- 4. END - PUBLIC METHODS ------------------
/// -----------------------------------------------------------

void normalizeUpperCst(Data &myD)
{
  for ( int i = myD.N-2 ; i >= 0 ; --i )
  {
    if ( myD.upperCst[i] > myD.upperCst[i+1] )
      myD.upperCst[i] = myD.upperCst[i+1] ;
  }
}
