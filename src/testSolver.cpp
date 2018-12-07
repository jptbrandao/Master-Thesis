#include "testSolver.hpp"
#include "scopt-ext.h"
#include <algorithm>
#include "mosek.h"

/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
    printf("%s",str);
} /* printstr */

Solution test(Data &myD)
{
  std::cout << "testing"<< std::endl;

  Solution sol(3);
  return sol;
}  

Solution solveGreedy2WithMosek(Data &myData, double lambda)
{
  Solution opt(myData.N);

  char buffer[MSK_MAX_STR_LEN];
  MSKenv_t env;
  MSKrescodee r;
  MSKtask_t task;
  void* sch;

  /* Specify nonlinear terms in the objective. */
  double * oprfo = new double [myData.N] ;
  double * oprgo = new double [myData.N] ;
  double * oprho = new double [myData.N] ;
  int * opro = new int [myData.N] ;
  int * oprjo = new int [myData.N] ;
  for (int i=0 ; i < myData.N ; i++)
  {
    // first term is the type of funtion 3 for polynoms
    // second term is the constant in front of polynom
    // third term is the exponent
    // fourth i don't remember, stays zero
    opro[i]  = 3 ;
    oprfo[i] = 0 ;
    oprgo[i] = 2 ;
    oprho[i] = 0 ;
    oprjo[i] = i ;
  }

  // Replace this by definition of the hyperplanes
  int *aptrb, *aptre, *asub;
  double *aval;
    
  // First Non-zero element in column
  // Start from the top left element of the matrix, count downwards
  // Let COUNT be the number of non-zero elements until some point in the matrix
  // The first non-zero element will have COUNT == 0 -> aptrb[0] = COUNT
  // Go on to the next column, find the first non-zero element -> aptrb[1] = COUNT
  // i.e. aptrb[i] = COUNT
  aptrb = new int [myData.N] ;

  // counts how many non zero numbers in each column
  aptre = new int [myData.N] ;
    
  int count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    bool isFirstNonZero = true;
    // the + 2 refers to the Resource B1 and Hyperpane B2 restrictions
    for ( int row = 0 ; row < myData.M /*+ 1*/ ; ++row )
    {
      //if ( row == myData.M )
      //{
      //  count++;
      //  continue;
      //}

      if ( myData.nestCst[row] && row >= col )
      {
        if ( isFirstNonZero == true )
        {
          aptrb[col] = count;
          isFirstNonZero = false;
        }
        count++;
      }
    }
    aptre[col] = count;
  }
  
  //TODO figure out what asub and aval mean
  asub = new int[count];
  aval = new double[count];
  count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    for ( int row = 0 ; row < myData.M /*+ 1*/ ; ++row )
    {
      if ( myData.nestCst[row] && row >= col )
      {
        aval[count] = 1;
        asub[count] = row;
        count++;
      }
    }
  }

  double * blc = new double [myData.M] ;
  double * buc = new double [myData.M] ;
  MSKboundkeye * bkc = new MSKboundkeye [myData.M] ;

  for ( int i = 0 ; i < myData.M-1; i++ )
  {
    bkc[i] = MSK_BK_UP;
    blc[i] = -MSK_INFINITY;
    buc[i] = myData.upperCst[i];
  }

  bkc[myData.M-1] = MSK_BK_FX;
  blc[myData.M-1] = myData.B1;
  buc[myData.M-1] = myData.B1;

    
  // linear term
  double *c = new double [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    c[i] = myData.c[i] + lambda * myData.p[i];
  }

  double * blx = new double [myData.N] ;
  double * bux = new double [myData.N] ;
  MSKboundkeye * bkx = new MSKboundkeye [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    blx[i] = 0;
    bux[i] = MSK_INFINITY;
    bkx[i] = MSK_BK_LO;
  }

  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);
    
  //    r = MSK_putexitfunc(env, exitFunc, NULL);
    
  if ( r==MSK_RES_OK )
  {
    /* Make the optimization task. */
    r = MSK_makeemptytask(env,&task);
    if ( r==MSK_RES_OK )
      MSK_linkfunctotaskstream(task,MSK_STREAM_MSG,NULL,printstr);
        
    if ( r==MSK_RES_OK )
    {
      //            checkData(myData);
      r = MSK_inputdata(task,
                        myData.M, myData.N,
                        myData.M, myData.N,
                        c, 0,
                        aptrb,aptre,
                        asub,aval,
                        bkc,blc,buc,
                        bkx,blx,bux);
    }
        
    if ( r== MSK_RES_OK )
    {
      /* Set-up of nonlinear expressions. */
      r = MSK_scbegin(task,
                      myData.N,opro,oprjo,oprfo,oprgo,oprho,
                      0,NULL,NULL,NULL,NULL,NULL,NULL,
                      &sch);
        
      if ( r==MSK_RES_OK )
      {
        // No log to avoid any overhead related to outputs
        MSK_putnaintparam(task, "MSK_IPAR_LOG", 0);
            
        // Only one single thread
        MSK_putnaintparam(task, "MSK_IPAR_NUM_THREADS", 1);
                
        // Setting a maximum time value for MOSEK
        // cannot set float time for MOSEK
        MSK_putnadouparam( task, "MSK_DPAR_OPTIMIZER_MAX_TIME", (double)1000000/*((double) myData.cpu_time_limit + (double)myData.time_StartComput - (double)clock())/(double) CLOCKS_PER_SEC */);
                
        MSK_writedata(task,"taskdump.opf");
        clock_t mosekBeg = clock();
        r = MSK_optimize(task);
        clock_t mosekEnd = clock();
             
        opt.setTime( double ( mosekEnd - mosekBeg ) / CLOCKS_PER_SEC );
                
        // if the run has been aborted (in fact its because of time)
        // then we don't print the solution
        if (r!=MSK_RES_OK)
        {
          std::cout << "didn't optimize: " << r << std::endl;
          return opt; //throw exceptionTime() ;
        }
                
        MSK_solutionsummary(task,MSK_STREAM_MSG);
           
        MSKstakeye useless ;
        MSKrealt resulti ;
        MSKrealt temp ;
                
        for (int i=0 ; i < myData.N ; i++)
        {
          MSK_getsolutioni(task,MSK_ACC_VAR,i,MSK_SOL_ITR,&useless,&resulti,&temp,&temp,&temp);
                  
          opt[i] = (double)resulti ;
        }
      }
            
      /* The nonlinear expressions are no longer needed. */
      MSK_scend(task,&sch);
    }
    MSK_deletetask(&task);
  }

  MSK_deleteenv(&env);
  delete [] oprfo ;
  delete [] oprgo ;
  delete [] oprho ;
  delete [] opro ;
  delete [] oprjo ;
  delete [] aptrb ;
  delete [] aptre ;
  delete [] asub ;
  delete [] aval ;
  delete [] blc ;
  delete [] buc ;
  delete [] bkc ;
  delete [] blx ;
  delete [] bux ;
  delete [] bkx ;
    
  if ( r!=MSK_RES_OK )
  {
    MSK_getcodedesc(r,buffer,NULL); 
    printf("Description: %s\n",buffer); 
  }
  
  return opt;
}


Solution solveGreedyWithMosek(Data &myData)
{
  Solution opt(myData.N);

  char buffer[MSK_MAX_STR_LEN];
  MSKenv_t env;
  MSKrescodee r;
  MSKtask_t task;
  void* sch;

  /* Specify nonlinear terms in the objective. */
  double * oprfo = new double [myData.N] ;
  double * oprgo = new double [myData.N] ;
  double * oprho = new double [myData.N] ;
  int * opro = new int [myData.N] ;
  int * oprjo = new int [myData.N] ;
  for (int i=0 ; i < myData.N ; i++)
  {
    // first term is the type of funtion 3 for polynoms
    // second term is the constant in front of polynom
    // third term is the exponent
    // fourth i don't remember, stays zero
    opro[i]  = 3 ;
    oprfo[i] = 0 ;
    oprgo[i] = 2 ;
    oprho[i] = 0 ;
    oprjo[i] = i ;
  }

  // Replace this by definition of the hyperplanes
  int *aptrb, *aptre, *asub;
  double *aval;
    
  // First Non-zero element in column
  // Start from the top left element of the matrix, count downwards
  // Let COUNT be the number of non-zero elements until some point in the matrix
  // The first non-zero element will have COUNT == 0 -> aptrb[0] = COUNT
  // Go on to the next column, find the first non-zero element -> aptrb[1] = COUNT
  // i.e. aptrb[i] = COUNT
  aptrb = new int [myData.N] ;

  // counts how many non zero numbers in each column
  aptre = new int [myData.N] ;
    
  int count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    bool isFirstNonZero = true;
    // the + 2 refers to the Resource B1 and Hyperpane B2 restrictions
    for ( int row = 0 ; row < myData.M /*+ 1*/ ; ++row )
    {
      //if ( row == myData.M )
      //{
      //  count++;
      //  continue;
      //}

      if ( myData.nestCst[row] && row >= col )
      {
        if ( isFirstNonZero == true )
        {
          aptrb[col] = count;
          isFirstNonZero = false;
        }
        count++;
      }
    }
    aptre[col] = count;
  }
  
  //TODO figure out what asub and aval mean
  asub = new int[count];
  aval = new double[count];
  count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    for ( int row = 0 ; row < myData.M /*+ 1*/ ; ++row )
    {
      if ( myData.nestCst[row] && row >= col )
      {
        aval[count] = 1;
        asub[count] = row;
        count++;
      }
    }
  }

  double * blc = new double [myData.M] ;
  double * buc = new double [myData.M] ;
  MSKboundkeye * bkc = new MSKboundkeye [myData.M] ;

  for ( int i = 0 ; i < myData.M-1; i++ )
  {
    bkc[i] = MSK_BK_UP;
    blc[i] = -MSK_INFINITY;
    buc[i] = myData.upperCst[i];
  }

  bkc[myData.M-1] = MSK_BK_FX;
  blc[myData.M-1] = myData.B1;
  buc[myData.M-1] = myData.B1;

    
  // linear term
  double *c = new double [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    c[i] = myData.c[i];
  }

  double * blx = new double [myData.N] ;
  double * bux = new double [myData.N] ;
  MSKboundkeye * bkx = new MSKboundkeye [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    blx[i] = 0;
    bux[i] = MSK_INFINITY;
    bkx[i] = MSK_BK_LO;
  }

  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);
    
  //    r = MSK_putexitfunc(env, exitFunc, NULL);
    
  if ( r==MSK_RES_OK )
  {
    /* Make the optimization task. */
    r = MSK_makeemptytask(env,&task);
    if ( r==MSK_RES_OK )
      MSK_linkfunctotaskstream(task,MSK_STREAM_MSG,NULL,printstr);
        
    if ( r==MSK_RES_OK )
    {
      //            checkData(myData);
      r = MSK_inputdata(task,
                        myData.M, myData.N,
                        myData.M, myData.N,
                        c, 0,
                        aptrb,aptre,
                        asub,aval,
                        bkc,blc,buc,
                        bkx,blx,bux);
    }
        
    if ( r== MSK_RES_OK )
    {
      /* Set-up of nonlinear expressions. */
      r = MSK_scbegin(task,
                      myData.N,opro,oprjo,oprfo,oprgo,oprho,
                      0,NULL,NULL,NULL,NULL,NULL,NULL,
                      &sch);
        
      if ( r==MSK_RES_OK )
      {
        // No log to avoid any overhead related to outputs
        MSK_putnaintparam(task, "MSK_IPAR_LOG", 0);
            
        // Only one single thread
        MSK_putnaintparam(task, "MSK_IPAR_NUM_THREADS", 1);
                
        // Setting a maximum time value for MOSEK
        // cannot set float time for MOSEK
        MSK_putnadouparam( task, "MSK_DPAR_OPTIMIZER_MAX_TIME", (double)1000000/*((double) myData.cpu_time_limit + (double)myData.time_StartComput - (double)clock())/(double) CLOCKS_PER_SEC */);
                
        MSK_writedata(task,"taskdump.opf");
        clock_t mosekBeg = clock();
        r = MSK_optimize(task);
        clock_t mosekEnd = clock();
             
        opt.setTime( double ( mosekEnd - mosekBeg ) / CLOCKS_PER_SEC );
                
        // if the run has been aborted (in fact its because of time)
        // then we don't print the solution
        if (r!=MSK_RES_OK)
        {
          std::cout << "didn't optimize: " << r << std::endl;
          return opt; //throw exceptionTime() ;
        }
                
        MSK_solutionsummary(task,MSK_STREAM_MSG);
           
        MSKstakeye useless ;
        MSKrealt resulti ;
        MSKrealt temp ;
                
        for (int i=0 ; i < myData.N ; i++)
        {
          MSK_getsolutioni(task,MSK_ACC_VAR,i,MSK_SOL_ITR,&useless,&resulti,&temp,&temp,&temp);
                  
          opt[i] = (double)resulti ;
        }
      }
            
      /* The nonlinear expressions are no longer needed. */
      MSK_scend(task,&sch);
    }
    MSK_deletetask(&task);
  }

  MSK_deleteenv(&env);
  delete [] oprfo ;
  delete [] oprgo ;
  delete [] oprho ;
  delete [] opro ;
  delete [] oprjo ;
  delete [] aptrb ;
  delete [] aptre ;
  delete [] asub ;
  delete [] aval ;
  delete [] blc ;
  delete [] buc ;
  delete [] bkc ;
  delete [] blx ;
  delete [] bux ;
  delete [] bkx ;
    
  if ( r!=MSK_RES_OK )
  {
    MSK_getcodedesc(r,buffer,NULL); 
    printf("Description: %s\n",buffer); 
  }
  
  return opt;
}


/**
 * Solving With MOSEK
*/
Solution solveMosek(Data &myData)
{
  Solution opt(myData.N);

  char buffer[MSK_MAX_STR_LEN];
  MSKenv_t env;
  MSKrescodee r;
  MSKtask_t task;
  void* sch;

  /* Specify nonlinear terms in the objective. */
  double * oprfo = new double [myData.N] ;
  double * oprgo = new double [myData.N] ;
  double * oprho = new double [myData.N] ;
  int * opro = new int [myData.N] ;
  int * oprjo = new int [myData.N] ;
  for (int i=0 ; i < myData.N ; i++)
  {
    // first term is the type of funtion 3 for polynoms
    // second term is the constant in front of polynom
    // third term is the exponent
    // fourth i don't remember, stays zero
    opro[i]  = 3 ;
    oprfo[i] = 0 ;
    oprgo[i] = 2 ;
    oprho[i] = 0 ;
    oprjo[i] = i ;
  }

  // Replace this by definition of the hyperplanes
  int *aptrb, *aptre, *asub;
  double *aval;
    
  // First Non-zero element in column
  // Start from the top left element of the matrix, count downwards
  // Let COUNT be the number of non-zero elements until some point in the matrix
  // The first non-zero element will have COUNT == 0 -> aptrb[0] = COUNT
  // Go on to the next column, find the first non-zero element -> aptrb[1] = COUNT
  // i.e. aptrb[i] = COUNT
  aptrb = new int [myData.N] ;

  // counts how many non zero numbers in each column
  aptre = new int [myData.N] ;
    
  int count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    bool isFirstNonZero = true;
    // the + 2 refers to the Resource B1 and Hyperpane B2 restrictions
    for ( int row = 0 ; row < myData.M + 1 ; ++row )
    {
      if ( row == myData.M )
      {
        count++;
        continue;
      }

      if ( myData.nestCst[row] && row >= col )
      {
        if ( isFirstNonZero == true )
        {
          aptrb[col] = count;
          isFirstNonZero = false;
        }
        count++;
      }
    }
    aptre[col] = count;
  }
  

  asub = new int[count];
  aval = new double[count];
  count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    for ( int row = 0 ; row < myData.M + 1 ; ++row )
    {
      if ( row == myData.M )
      {
        aval[count] = myData.p[col];
        asub[count] = row;
        count++;
        continue;
      }
      if ( myData.nestCst[row] && row >= col )
      {
        aval[count] = 1;
        asub[count] = row;
        count++;
      }
    }
  }

  double * blc = new double [myData.M+1] ;
  double * buc = new double [myData.M+1] ;
  MSKboundkeye * bkc = new MSKboundkeye [myData.M+2] ;

  for ( int i = 0 ; i < myData.M-1; i++ )
  {
    bkc[i] = MSK_BK_UP;
    blc[i] = -MSK_INFINITY;
    buc[i] = myData.upperCst[i];
  }

  bkc[myData.M-1] = MSK_BK_FX;
  blc[myData.M-1] = myData.B1;
  buc[myData.M-1] = myData.B1;

  bkc[myData.M] = MSK_BK_FX;
  blc[myData.M] = myData.B2;
  buc[myData.M] = myData.B2;
    
  // linear term
  double *c = new double [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    c[i] = myData.c[i];
  }

  double * blx = new double [myData.N] ;
  double * bux = new double [myData.N] ;
  MSKboundkeye * bkx = new MSKboundkeye [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    blx[i] = 0;
    bux[i] = MSK_INFINITY;
    bkx[i] = MSK_BK_LO;
  }



  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);
    
  //    r = MSK_putexitfunc(env, exitFunc, NULL);
    
  if ( r==MSK_RES_OK )
  {
    /* Make the optimization task. */
    r = MSK_makeemptytask(env,&task);
    if ( r==MSK_RES_OK )
      MSK_linkfunctotaskstream(task,MSK_STREAM_MSG,NULL,printstr);
        
    if ( r==MSK_RES_OK )
    {
      //            checkData(myData);
      r = MSK_inputdata(task,
                        myData.M+1,myData.N,
                        myData.M+1,myData.N,
                        c, 0,
                        aptrb,aptre,
                        asub,aval,
                        bkc,blc,buc,
                        bkx,blx,bux);
    }
        
    if ( r== MSK_RES_OK )
    {
      /* Set-up of nonlinear expressions. */
      r = MSK_scbegin(task,
                      myData.N,opro,oprjo,oprfo,oprgo,oprho,
                      0,NULL,NULL,NULL,NULL,NULL,NULL,
                      &sch);
        
      if ( r==MSK_RES_OK )
      {
        // No log to avoid any overhead related to outputs
        MSK_putnaintparam(task, "MSK_IPAR_LOG", 0);
            
        // Only one single thread
        MSK_putnaintparam(task, "MSK_IPAR_NUM_THREADS", 1);
                
        // Setting a maximum time value for MOSEK
        // cannot set float time for MOSEK
        MSK_putnadouparam( task, "MSK_DPAR_OPTIMIZER_MAX_TIME", (double)100000/*((double) myData.cpu_time_limit + (double)myData.time_StartComput - (double)clock())/(double) CLOCKS_PER_SEC */);
                
        MSK_writedata(task,"taskdump.opf");
        clock_t mosekBeg = clock();
        r = MSK_optimize(task);
        clock_t mosekEnd = clock();
             
        opt.setTime( double ( mosekEnd - mosekBeg ) / CLOCKS_PER_SEC );
                
        // if the run has been aborted (in fact its because of time)
        // then we don't print the solution
        if (r!=MSK_RES_OK)
        {
          std::cout << "didn't optimize: " << r << std::endl;
          return opt; //throw exceptionTime() ;
        }
                
        MSK_solutionsummary(task,MSK_STREAM_MSG);
           
        MSKstakeye useless ;
        MSKrealt resulti ;
        MSKrealt temp ;
                
        for (int i=0 ; i < myData.N ; i++)
        {
          MSK_getsolutioni(task,MSK_ACC_VAR,i,MSK_SOL_ITR,&useless,&resulti,&temp,&temp,&temp);
                  
          opt[i] = (double)resulti ;
        }
      }
            
      /* The nonlinear expressions are no longer needed. */
      MSK_scend(task,&sch);
    }
    MSK_deletetask(&task);
  }

  MSK_deleteenv(&env);
  delete [] oprfo ;
  delete [] oprgo ;
  delete [] oprho ;
  delete [] opro ;
  delete [] oprjo ;
  delete [] aptrb ;
  delete [] aptre ;
  delete [] asub ;
  delete [] aval ;
  delete [] blc ;
  delete [] buc ;
  delete [] bkc ;
  delete [] blx ;
  delete [] bux ;
  delete [] bkx ;
    
  if ( r!=MSK_RES_OK )
  {
    MSK_getcodedesc(r,buffer,NULL); 
    printf("Description: %s\n",buffer); 
  }
  
  return opt;
}

void solveMosekRepeat(Data &myData, double &totalDuration, int numRepetitions) 
{
  char buffer[MSK_MAX_STR_LEN];
  MSKenv_t env;
  MSKrescodee r;
  MSKtask_t task;
  void* sch;

  /* Specify nonlinear terms in the objective. */
  double * oprfo = new double [myData.N] ;
  double * oprgo = new double [myData.N] ;
  double * oprho = new double [myData.N] ;
  int * opro = new int [myData.N] ;
  int * oprjo = new int [myData.N] ;
  for (int i=0 ; i < myData.N ; i++)
  {
    // first term is the type of funtion 3 for polynoms
    // second term is the constant in front of polynom
    // third term is the exponent
    // fourth i don't remember, stays zero
    opro[i]  = 3 ;
    oprfo[i] = 0 ;
    oprgo[i] = 2 ;
    oprho[i] = 0 ;
    oprjo[i] = i ;
  }

  // Replace this by definition of the hyperplanes
  int *aptrb, *aptre, *asub;
  double *aval;
    
  // First Non-zero element in column
  // Start from the top left element of the matrix, count downwards
  // Let COUNT be the number of non-zero elements until some point in the matrix
  // The first non-zero element will have COUNT == 0 -> aptrb[0] = COUNT
  // Go on to the next column, find the first non-zero element -> aptrb[1] = COUNT
  // i.e. aptrb[i] = COUNT
  aptrb = new int [myData.N] ;

  // counts how many non zero numbers in each column
  aptre = new int [myData.N] ;
    
  int count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    bool isFirstNonZero = true;
    // the + 2 refers to the Resource B1 and Hyperpane B2 restrictions
    for ( int row = 0 ; row < myData.M + 1 ; ++row )
    {
      if ( row == myData.M )
      {
        count++;
        continue;
      }

      if ( myData.nestCst[row] && row >= col )
      {
        if ( isFirstNonZero == true )
        {
          aptrb[col] = count;
          isFirstNonZero = false;
        }
        count++;
      }
    }
    aptre[col] = count;
  }
  

  asub = new int[count];
  aval = new double[count];
  count = 0;
  for ( int col = 0 ; col < myData.N ; ++col )
  {
    for ( int row = 0 ; row < myData.M + 1 ; ++row )
    {
      if ( row == myData.M )
      {
        aval[count] = myData.p[col];
        asub[count] = row;
        count++;
        continue;
      }
      if ( myData.nestCst[row] && row >= col )
      {
        aval[count] = 1;
        asub[count] = row;
        count++;
      }
    }
  }

  double * blc = new double [myData.M+1] ;
  double * buc = new double [myData.M+1] ;
  MSKboundkeye * bkc = new MSKboundkeye [myData.M+2] ;

  for ( int i = 0 ; i < myData.M-1; i++ )
  {
    bkc[i] = MSK_BK_UP;
    blc[i] = -MSK_INFINITY;
    buc[i] = myData.upperCst[i];
  }

  bkc[myData.M-1] = MSK_BK_FX;
  blc[myData.M-1] = myData.B1;
  buc[myData.M-1] = myData.B1;

  bkc[myData.M] = MSK_BK_FX;
  blc[myData.M] = myData.B2;
  buc[myData.M] = myData.B2;
    
  // linear term
  double *c = new double [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    c[i] = myData.c[i];
  }

  double * blx = new double [myData.N] ;
  double * bux = new double [myData.N] ;
  MSKboundkeye * bkx = new MSKboundkeye [myData.N] ;
  for ( int i = 0 ; i < myData.N ; ++i )
  {
    blx[i] = 0;
    bux[i] = MSK_INFINITY;
    bkx[i] = MSK_BK_LO;
  }



  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);
    
  //    r = MSK_putexitfunc(env, exitFunc, NULL);
    
  if ( r==MSK_RES_OK )
  {
    /* Make the optimization task. */
    r = MSK_makeemptytask(env,&task);
    if ( r==MSK_RES_OK )
      MSK_linkfunctotaskstream(task,MSK_STREAM_MSG,NULL,printstr);
        
    if ( r==MSK_RES_OK )
    {
      //            checkData(myData);
      r = MSK_inputdata(task,
                        myData.M+1,myData.N,
                        myData.M+1,myData.N,
                        c, 0,
                        aptrb,aptre,
                        asub,aval,
                        bkc,blc,buc,
                        bkx,blx,bux);
    }
        
    if ( r== MSK_RES_OK )
    {
      /* Set-up of nonlinear expressions. */
      r = MSK_scbegin(task,
                      myData.N,opro,oprjo,oprfo,oprgo,oprho,
                      0,NULL,NULL,NULL,NULL,NULL,NULL,
                      &sch);
        
      if ( r==MSK_RES_OK )
      {
        // No log to avoid any overhead related to outputs
        MSK_putnaintparam(task, "MSK_IPAR_LOG", 0);
            
        // Only one single thread
        MSK_putnaintparam(task, "MSK_IPAR_NUM_THREADS", 1);
                
        // Setting a maximum time value for MOSEK
        // cannot set float time for MOSEK
        MSK_putnadouparam( task, "MSK_DPAR_OPTIMIZER_MAX_TIME", (double)100000/*((double) myData.cpu_time_limit + (double)myData.time_StartComput - (double)clock())/(double) CLOCKS_PER_SEC */);
                
        MSK_writedata(task,"taskdump.opf");
        clock_t mosekBeg = clock();
        for ( int i = 0 ; i < numRepetitions; ++i )
          r = MSK_optimize(task);
        clock_t mosekEnd = clock();
             
        totalDuration = (double) ( mosekEnd - mosekBeg ) / CLOCKS_PER_SEC ;
                
        // if the run has been aborted (in fact its because of time)
        // then we don't print the solution
        if (r!=MSK_RES_OK)
        {
          std::cout << "didn't optimize: " << r << std::endl;
          return; //throw exceptionTime() ;
        }
                
        MSK_solutionsummary(task,MSK_STREAM_MSG);
           
        MSKstakeye useless ;
        MSKrealt resulti ;
        MSKrealt temp ;
        
      }
            
      /* The nonlinear expressions are no longer needed. */
      MSK_scend(task,&sch);
    }
    MSK_deletetask(&task);
  }

  MSK_deleteenv(&env);
  delete [] oprfo ;
  delete [] oprgo ;
  delete [] oprho ;
  delete [] opro ;
  delete [] oprjo ;
  delete [] aptrb ;
  delete [] aptre ;
  delete [] asub ;
  delete [] aval ;
  delete [] blc ;
  delete [] buc ;
  delete [] bkc ;
  delete [] blx ;
  delete [] bux ;
  delete [] bkx ;
    
  if ( r!=MSK_RES_OK )
  {
    MSK_getcodedesc(r,buffer,NULL); 
    printf("Description: %s\n",buffer); 
  }
}
