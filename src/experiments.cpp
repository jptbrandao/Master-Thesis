#include "experiments.hpp"
#include <algorithm>
#include "commonFuncs.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include "newQS.hpp"
#include "weaklypoly.hpp"
#include <random>
#include <unordered_set>
#include <iomanip>
#include <chrono>
#include <ctime>

#define maxValue numeric_limits<double>::max()

using namespace std;

// read data
// transform into Data structure (including filling the boolean vector on upper nested cst)
// generate random B2 values
// test for feasibility
// save data

void runExperiment(string problemsFilePath, string resultsFilePath);

std::pair<vector<double>, vector<double> > readVidalData(string filepath);
void generateData(Data &myD);
void formatIntoDataStructure(Data &myD, vector<double> cost, vector<double> upperCst);
void saveData(Data &myD, string filepath, string filename);
void manageExperiments(string readDir, string writeName);

void createAndSaveAllInstances();
void createIncreasingInstances(Data &myD1, Data &myD2, std::mt19937 &gen);
void testIncreasingData(Data &myD1, Data &myD2);

void runOneExperiment(string readDir, string writeName, int probSize, int instance);
void manageExperiments(string readDir, string writeName)

void quicktest()
{
  cout << "quick test begun" << endl;
  //runExperiment("../Data/Brandao Instances/RAP-NC-L_0000.txt", "../Data/Results/0000.txt");
  //manageExperiments();
  
  //manageExperiments("Brandao Instances", "Results_Original");
  //manageExperiments("Increasing Instances", "Results_Increasing");
  //manageExperiments("Decreasing Instances", "Results_Decreasing");
  runOneExperiment("Increasing Instances", "extraExp", 2, 5);
  runOneExperiment("Increasing Instances", "extraExp", 2, 6);
  runOneExperiment("Increasing Instances", "extraExp", 2, 7);
  runOneExperiment("Increasing Instances", "extraExp", 2, 8);

  //createAndSaveAllInstances();

  /*// Testing generating data
  std::random_device rd;
  std::mt19937 gen(rd());

  int numVar = 10;
  Data myD1(numVar, numVar);
  Data myD2(numVar, numVar);

  createIncreasingInstances(myD1, myD2, gen);
  testIncreasingData(myD1, myD2);*/

  cout << "quick test has ended" << endl;
}

void writeToOutput(string filepath, string text);
void manageExperiments(string readDir, string writeName)
{
  string readfilepath = "../../Data/" + readDir + "/RAP-NC-L_";
  string writefilepath = "Results/" + writeName + ".txt";

  string header = "Date: 24/1/19\nFormat:\nProblem ID\nWeakly: \t Time - Obj - PercentDif\nStrongly: \t Time - Obj - PercentDif\nMosek: \t\t Time - Obj\n\n";
  writeToOutput(writefilepath, header);
  
  for ( int probSize = 0 ; probSize < 16 ; probSize++ )
  {
    for ( int instance = 0 ; instance < 10 ; instance++ )
    {
      string extraProbSize = ( probSize < 10 ) ? "0" : "";
      string probID = extraProbSize + to_string(probSize) + "0" + to_string(instance) + ".txt";
      
      cout << endl << "----------------------------------" << endl;
      cout << "ProbID: " << probID << endl;
      writeToOutput(writefilepath, probID);
      runExperiment(readfilepath + probID, writefilepath);
    }
  }
}

void runOneExperiment(string readDir, string writeName, int probSize, int instance)
{
  string readfilepath = "../../Data/" + readDir + "/RAP-NC-L_";
  //string readfilepath = "../Data/Increasing Instances/RAP-NC-L_";
  string writefilepath = "Results/" + writeName + ".txt";

  string header = "Date: 25/1/19\nFormat:\nProblem ID\nWeakly: \t Time - Obj - PercentDif\nStrongly: \t Time - Obj - PercentDif\nMosek: \t\t Time - Obj\n\n";
  writeToOutput(writefilepath, header);
  
  string extraProbSize = ( probSize < 10 ) ? "0" : "";
  string probID = extraProbSize + to_string(probSize) + "0" + to_string(instance) + ".txt";
      
  cout << endl << "----------------------------------" << endl;
  cout << "ProbID: " << probID << endl;
  writeToOutput(writefilepath, probID);
  runExperiment(readfilepath + probID, writefilepath);
}

void readAndGenerateProblemInstances()
{
  string readfilepath = "../Data/Vidal Instances/";
  string writefilepath = "../Data/Brandao Instances/";

  for ( int probSize = 0 ; probSize < 16 ; probSize++ )
  {
    for ( int instance = 0 ; instance < 10 ; instance++ )
    {
      cout << endl << "ProbSize: " << probSize << " Instance: " << instance << endl;
      string readfilename = "RAP-NC_";
      string extraProbSize = ( probSize < 10 ) ? "0" : "";
      string fileID = extraProbSize + to_string(probSize) + "0" + to_string(instance) + ".txt";
      readfilename += fileID;

      pair<vector<double>, vector<double> > vidalData = readVidalData(readfilepath + readfilename);
      int n = vidalData.first.size();
      Data myD(n, n);
      formatIntoDataStructure(myD, vidalData.first, vidalData.second);
      generateData(myD);

      string writefilename = "RAP-NC-L_";
      writefilename += fileID;
      saveData(myD, writefilepath + writefilename, writefilename); 
    }
  }
}

void saveData(Data &myD, string filepath, string filename)
{
  ofstream file(filepath, ios::out);

  file << "PROBLEM_NAME: " << filename << endl;
  file << "PROBLEM_TYPE: UB" << endl;
  file << "NB_RESOURCES: " << myD.N << endl;
  file << "LINEAR_TERM: " << myD.B2 << endl;
  file << "BOUNDS SECTION: INDEX \t UB \t LINEAR_TERM" << endl;
  
  std::cout << std::fixed;

  for ( int i = 0 ; i < myD.N ; ++i )
  {
    string extra = (myD.p[i] > 0 ) ? " " : "";
    file << i << " \t " << std::setprecision(8) << myD.upperCst[i] << " \t " + extra << std::setprecision(8) << myD.p[i] << endl;
  }

  file << "FUNCTION PROFILE: " << endl;

  for ( int i = 0 ; i < myD.N ; ++i )
    file << i << " \t " << std::setprecision(12) << myD.c[i] << endl;

  file << "# END OF FILE #";
}

void formatIntoDataStructure(Data &myD, vector<double> cost, vector<double> upperCst)
{
  for ( int i = 0 ; i < myD.N ; ++i )
  {
    myD.c[i] = cost[i];
    myD.upperCst[i] = upperCst[i];
    myD.nestCst[i] = true;
  }

  myD.B1 = upperCst[myD.N-1];
}

pair<double, double> generateLinearTerms(Data &myD, std::mt19937 &gen)
{
  int numVar = myD.N;
  std::uniform_real_distribution<> pDis(-1.0, 1.0);

  std::unordered_set<double> pSet;
  pair<double, double> limits(maxValue, -maxValue);
  for ( int i = 0 ; i < numVar ; ++i )
  {
    double randomNumber = pDis(gen);
    std::unordered_set<double>::const_iterator it = pSet.find(randomNumber);
    while ( it == pSet.end() && randomNumber == 0 )
      randomNumber = pDis(gen);

    if ( limits.first > randomNumber )
      limits.first = randomNumber;

    if ( limits.second < randomNumber )
      limits.second = randomNumber;

    pSet.insert(randomNumber);
    myD.p[i] = randomNumber;
  }

  return limits;
}

bool generateB2(Data &myD, pair<double, double> limits, std::mt19937 &gen)
{
  bool feasibilityFlag = false;
  double minRange = limits.first*myD.B1;
  double maxRange = limits.second*myD.B1;

  ////cout << "initial limits: (" << minRange << ", " << maxRange << ")" << endl;
  while ( feasibilityFlag == false && fabs(maxRange - minRange) > 1 )
  {
    std::uniform_real_distribution<> B2Dis(minRange, maxRange);
    myD.B2 = B2Dis(gen);
    pair<double, double> limits = findBigLambdaInterval(myD);
    feasibilityFlag = checkProblemFeasibility(myD, limits);

    double subgrad = 0;
    if ( feasibilityFlag == false )
    {
      Solution x(myD.N);
      greedy(myD, x, limits.first);
      subgrad = myD.subgradient(x);

      if ( subgrad < 0 )
        maxRange = myD.B2;
      
      if ( subgrad > 0 )
        minRange = myD.B2;
    }

    //cout << "Limits: (" << minRange << ", " << maxRange << ")" << "  -  B2: " << myD.B2 << "  -  subgrad: " << subgrad << endl;
  }

  if ( fabs(maxRange - minRange) < 1 )
    return false;

  return true;
}

// First is B2, second is the linear terms p
void generateData(Data &myD)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  bool feasibilityFlag = false;
  while ( feasibilityFlag == false )
  {
    pair<double, double> limits =  generateLinearTerms(myD, gen);
    feasibilityFlag = generateB2(myD, limits, gen);
  }
}

// grib2
void generateB2Linear(Data &myD, std::mt19937 &gen)
{
  bool feasibilityFlag = false;
  while ( feasibilityFlag == false )
  {
    pair<double, double> limits =  generateLinearTerms(myD, gen);
    feasibilityFlag = generateB2(myD, limits, gen);
  }
}

// grib4
void createAndSaveAllInstances()
{
  string writefilepath1 = "../Data/Increasing Instances/";
  string writefilepath2 = "../Data/Decreasing Instances/";

  std::random_device rd;
  std::mt19937 gen(rd());

  int nSize[] = {10, 20, 50};
  int power = 0;
  int multPower = 1;
  for ( int probSize = 0 ; probSize < 16 ; probSize++ )
  {
    if ( probSize/3 > power )
    {
      multPower = multPower*10;
      power++;
    }
    int n = nSize[probSize%3] * multPower; 
    cout << "Problem Size: " << n << endl;

    for ( int instance = 0 ; instance < 10 ; instance++ )
    {
      cout << endl << "ProbSize: " << probSize << " Instance: " << instance << endl;
      string extraProbSize = ( probSize < 10 ) ? "0" : "";
      string fileID = extraProbSize + to_string(probSize) + "0" + to_string(instance) + ".txt";

      Data myD1(n, n);
      Data myD2(n, n);
      createIncreasingInstances(myD1, myD2, gen);

      string writefilename = "RAP-NC-L_";
      writefilename += fileID;
      saveData(myD1, writefilepath1 + writefilename, writefilename); 
      saveData(myD2, writefilepath2 + writefilename, writefilename);
    }
  }
}

// grib1
void createIncreasingInstances(Data &myD1, Data &myD2, std::mt19937 &gen)
{
  std::uniform_real_distribution<> alphaDis(0, 0.5);
  std::uniform_real_distribution<> cDis(0, 1);

  vector<double> alpha1(myD1.N);
  vector<double> alpha2(myD2.N);
  for ( int i = 0 ; i < myD1.N ; ++i )
  {
    alpha1[i] = alphaDis(gen);
    alpha2[i] = alphaDis(gen);
  }

  sort(alpha1.begin(), alpha1.end());
  sort(alpha2.begin(), alpha2.end());

  double cst_prev1 = 0;
  double cst_prev2 = 0;
  for ( int i = 0 ; i < myD1.N ; ++i )
  {
    myD1.c[i] = cDis(gen);
    myD1.upperCst[i] = cst_prev1 + alpha1[i];
    myD1.nestCst[i] = true;

    myD2.c[i] = cDis(gen);
    int j = myD1.N - i - 1;
    myD2.upperCst[i] = cst_prev2 + alpha2[j];
    myD2.nestCst[i] = true;

    cst_prev1 = myD1.upperCst[i];
    cst_prev2 = myD2.upperCst[i];
  }

  myD1.B1 = myD1.upperCst[myD1.N-1];
  myD2.B1 = myD2.upperCst[myD2.N-1];

  generateB2Linear(myD1, gen);
  generateB2Linear(myD2, gen);
}

// grib3
void testIncreasingData(Data &myD1, Data &myD2)
{
  pair<double, double> limits = findBigLambdaInterval(myD1);
  bool isFeasible = checkProblemFeasibility(myD1, limits);
  
  if ( isFeasible == false )
  {
    cout << "Problem1 not feasible " << endl;
    exit(1);
  }

  limits = findBigLambdaInterval(myD2);
  isFeasible = checkProblemFeasibility(myD2, limits);

  if ( isFeasible == false )
  {
    cout << "Problem2 not feasible " << endl;
    exit(1);
  }

  cout << "Problems feasible " << endl;
  cout << "1 - Distance between constraints: no. variables: " << myD1.N << endl;
  for ( int i = 0 ; i < myD1.N-1 ; ++i )
  {
    double dist = myD1.upperCst[i+1] - myD1.upperCst[i];
    cout << "b[" << i+1 << "] - b[" << i << "] = " << dist << endl;
  }

  cout << "2 - Distance between constraints: no. variables: " << myD2.N << endl;
  for ( int i = 0 ; i < myD2.N-1 ; ++i )
  {
    double dist = myD2.upperCst[i+1] - myD2.upperCst[i];
    cout << "b[" << i+1 << "] - b[" << i << "] = " << dist << endl;
  }

}

// First is cost, second is Upper Nested Constraint
pair<vector<double>, vector<double> > readVidalData(string filepath)
{
  ifstream file(filepath, ios::in);

  string line;

  getline(file, line);  // Program Name
  getline(file, line);  // Problem Type

  int numVar;
  getline(file, line, ':'); // Number of Variables
  file >> numVar;
  getline(file, line);  // Empty line

  getline(file, line);  // Function Type
  getline(file, line);  // Bounds Section

  vector<double> dataUpperCst(numVar);
  for ( int i = 0 ; i < numVar ; ++i )
  {
    int index;
    double  lowCst, uppCst, lowBound, uppBound;
    file >> index >> lowCst >> uppCst >> lowBound >> uppBound;
    dataUpperCst[index] = uppCst;
  }

  getline(file, line);  // Empty line
  getline(file, line);  // Function Profile

  vector<double> dataCost(numVar);
  for ( int i = 0 ; i < numVar ; ++i )
  {
    int index;
    double cost;
    file >> index >> cost;
    dataCost[index] = cost;
  } 

  return make_pair(dataCost, dataUpperCst);
}



// 0.1% tolerance
double tolerance = 0.001;

bool feasibilityCheck(Data &myD, vector<Solution> sols)
{
  for ( int i = 0 ; i < sols.size() ; ++i )
  {
    if ( myD.isFeasible(sols[i]) == false )
     return false; 
  }

  return true;
}

double objectiveValueCheck(Data &myD, vector<Solution> sols)
{
  // Checks if the objective values are all the same
  int size = sols.size();
  double value = myD.F(sols[0]);
  double avg = value/size;
  for ( int i = 1 ; i < size ; ++i )
  {
    value = myD.F(sols[i]);
    if ( fabs(value - myD.F(sols[i-1])) > tolerance )
      return maxValue;

    avg += value/size;
  }

  return avg;
}

string getCurrentTime()
{
  time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  return ctime(&now);
}

double runWeaklyAlgorithm(Data &myD, double &totalDuration, int numRepetitions) 
{
  // Setup
  /**numRepetitions = 1;
    sols.emplace_back(myD.N); // added 16/1/19
  //for (int i = 0 ; i < numRepetitions ; ++i )
  //{
  //  sols.emplace_back(myD.N);
  //}

  beginTime = clock();
  // Runs algorithm
  cout << "Weakly repetition no.: 0 - Time: " << getCurrentTime();
  weaklyPoly(myD, limits.first-1, limits.second+1, sols[0]);
  endTime = clock();
  totalDuration = ((double)(endTime - beginTime)/CLOCKS_PER_SEC) ;

  if ( totalDuration > 1 )
  {
    if ( feasibilityCheck(myD, sols) == false )
      return maxValue;

    return objectiveValueCheck(myD, sols);
  }*/
  // Setup
  pair<double, double> limits = findBigLambdaInterval(myD);
  vector<Solution> sols;
  
  for ( int i = 0 ; i < numRepetitions; ++i )
      sols.emplace_back(myD.N);

  clock_t beginTime, endTime;       // Snapshot of current times
  beginTime = clock();
  for ( int i = 0 ; i < numRepetitions; ++i )
  {
    // Runs algorithm
    cout << "Weakly repetition no.: " << i << " - Time: " << getCurrentTime();
    weaklyPoly(myD, limits.first-1, limits.second+1, sols[i]);
  }
  endTime = clock();
  totalDuration += ((double)(endTime - beginTime)/CLOCKS_PER_SEC) ;

  totalDuration = totalDuration/numRepetitions;

  // Returns maxValue is there is an infeasible solution
  if ( feasibilityCheck(myD, sols) == false )
    return maxValue;

  return objectiveValueCheck(myD, sols);
}


double runStronglyAlgorithm(Data &myD, double &totalDuration, int numRepetitions)
{
  if ( numRepetitions == 0 )
    return 0;
 
  // Setup
  vector<Solution> sols;
  vector<RAPLBLSolver> solver;
  double tol = 0.001;

  for ( int i = 0 ; i < numRepetitions; ++i )
    sols.emplace_back(myD.N);

  clock_t beginTime, endTime;       // Snapshot of current times
  beginTime = clock();
  for ( int i = 0 ; i < numRepetitions; ++i )
  {
    // Runs algorithm
    cout << "Strongly repetition no.: " << i << " - Time: " << getCurrentTime();
    RAPLBLSolver solver(myD, tol) ;
    solver.solve_StronglyPoly(sols[i]);
  }
  endTime = clock();
  totalDuration += ((double) endTime - beginTime)/CLOCKS_PER_SEC ;

  totalDuration = totalDuration/numRepetitions;

  // Returns maxValue is there is an infeasible solution
  if ( feasibilityCheck(myD, sols) == false )
    return maxValue;

  return objectiveValueCheck(myD, sols);
}


string formatOutput( vector<double> times, vector<double> objVal, double weakPer, double strongPer)
{
  string text      = ""       ;
  string delimiter = "\n" ;

  string weakError = "";
  if ( fabs(weakPer) > tolerance )
    weakError = "OFF" ;

  string strongError = "";
  if ( fabs(strongPer) > tolerance )
    strongError = "OFF" ;

  stringstream weakTime, strongTime;
  weakTime << fixed << setprecision(12) << times[0];
  strongTime << fixed << setprecision(12) << times[1];

  text += "Weakly: \t "   + weakTime.str()   + " - " + to_string(objVal[0]) + " - " + to_string(weakPer) + " \t " + weakError + delimiter ;
  text += "Strongly: \t " + strongTime.str() + " - " + to_string(objVal[1]) + " - " + to_string(strongPer) + " \t " + strongError + delimiter ;

  text += delimiter ;

  return text ;
}

void writeToOutput(string filepath, string text)
{
  fstream arq;
  arq.open(filepath, std::ios_base::app);
  arq << text << endl;
  arq.close();
}

double percentDiff(double target, double base) 
{
  if ( base == 0 )
    base = std::numeric_limits<double>::min();
  
  return (target - base)/base;
}

void runExperiment(string problemsFilePath, string resultsFilePath)
{
  string objError = "";
  
  // Run Setup
  // It is assumed data is ok, i.e. problem instance is feasible and follow the assumptions stated in the dissertation
  Data myD(problemsFilePath); 

  pair<double, double> limits = findBigLambdaInterval(myD);
  if ( checkProblemFeasibility(myD, limits) == false )
  {
    cout << "Problem not feasible: " << problemsFilePath << endl;
    string output = "##### PROBLEM NOT FEASIBLE #####\n";
    writeToOutput(resultsFilePath, output);
    return;
  }

  vector<double> avgTimes  (4, 0) ; // 1. Weakly, 2. Strongly, 3. Mosek
  vector<double> objVal (4)    ; // Same as times array
  int repetitions[2] = {10, 10};

  if ( myD.N >= 50000 )
  {
    repetitions[0] = 1;
    repetitions[1] = 1;
  }
  //if ( myD.N >= 500000 )
  //  repetitions[1] = 0;

  // Runs through the same problem several times
  cout << "Running Weakly"    << endl ;
  objVal[0] = runWeaklyAlgorithm   (myD, avgTimes[0], repetitions[0]) ;
  cout << "Finished Weakly"   << endl;

  cout << "Running Strongly"  << endl ;
  objVal[1] = runStronglyAlgorithm (myD, avgTimes[1], repetitions[1]) ;
  cout << "Finished Strongly" << endl ;

    
  double weakPer = percentDiff(objVal[0], objVal[1]);
  double strongPer = percentDiff(objVal[1], objVal[0]);

  if ( repetitions[1] == 0)
  {
    weakPer = 0;
    strongPer = 0;
  }
  else if ( repetitions[2] == 0 )
  {
    weakPer = percentDiff(objVal[0], objVal[1]);
    strongPer = percentDiff(objVal[1], objVal[0]);
  }
  
  if ( fabs(weakPer)  > tolerance )
  {
    objError += "Weakly Off by " + to_string(weakPer) + "\n";
    cout << objError << endl;
  }

  if ( fabs(strongPer) > tolerance )
  {
    objError += "Strongly Off by " + to_string(strongPer) + "\n";
    cout << objError << endl;
  }


  // Writes to output
  string output = formatOutput(avgTimes, objVal, weakPer, strongPer);
  writeToOutput(resultsFilePath, output);
  cout << "Written times and objVal to file" << endl;
  cout << "----------------------------------" << endl << endl;
  cout << endl << endl << objError << endl;
}
