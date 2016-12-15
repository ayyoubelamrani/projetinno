/*!
*  @file DemoPararealHeatPDE.cpp
*  @brief Demonstration code for parareal time integration of heat-like PDE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

// C++ headers
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdlib>

// MRG TIMEE headers
#include <TIMEE/Parareal.hpp>
#include <TIMEE/TIMEEIO.hpp>

// MRG third-party headers

// Demo global constants
#define MIN_ARGC 9

// Demo datatypes
typedef double T;
typedef int U;
const char* FMT_T = "%lf";
const char* FMT_U = "%d";


#define HELP "\n\
BRIEF\n\
  Integrates a FEM-discrete time-dependent PDE of the form\n\
          du/dt - Lu = f in Omega,\n\
                   u = g on Gamma1,\n\
               du/dn = h on Gamma2,\n\
  using the Parareal time discretization scheme.\n\
ARGS\n\
  <mat_M_fname>        # Mass matrix filename (MTX file).\n\
  <mat_K_fname>        # Stiffness matrix filename (MTX file).\n\
  <vec_F_fname>        # Source vector filename (MTX file).\n\
  <vec_H_fname>        # Neumann BC vector filename (MTX file).\n\
  <dirich_dof_fname>   # Dirichlet boundary DOF filename (MRG A1D file).\n\
  <dirich_val_fname>   # Dirichlet boundary DOF value filename\n\
                         (MRG A1D file).\n\
  <vec_U_fpref>        # Solution vector filename prefix for all time steps.\n\
  <vec_U_fsuff>        # Solution vector filename suffix for all time steps.\n\
  [-i <vec_U0_fname>]  # Initial value vector filename (MTX file).\n\
  [-f {MTX|ENSI|A1D}]    # Solution vector file format (default: MTX).\n\
  [-t <fine_time_step>]     # Fine time step (default: 0.001).\n\
  [-n <numb_step_per_win>]  # Number of time steps per window (default: 10).\n\
  [-s {BEULER|TRULE}]    # Fine time discretization scheme (default: BEULER).\n\
  [-S {BEULER|TRULE}]  # Coarse time discretization scheme (default: BEULER).\n\
  [-e <res_thresh>]    # Residual threshold (default: 1e-6).\n\
  [-h]                 # Displays this description.\n\
"
int main (int argc, char *argv[]) {

  // ---------------------------------------------------------------------------
  // -- variables
  // ---------------------------------------------------------------------------

  // -- error
  int err;

  // -- measurement
  clock_t time_begin;
  clock_t time_end;
  double time_exec;

  // -- files
  char mat_M_fname[TIMEE_BUFSIZE_L];
  const char* mat_M_fformat = "mtx";
  char mat_K_fname[TIMEE_BUFSIZE_L];
  const char* mat_K_fformat = "mtx";
  
  char vec_F_fname[TIMEE_BUFSIZE_L];
  const char* vec_F_fformat = "mtx";
  char vec_H_fname[TIMEE_BUFSIZE_L];
  const char* vec_H_fformat = "mtx";
  
  char dirich_dof_fname[TIMEE_BUFSIZE_L];
  const char* dirich_dof_fformat = "a1d";
  const char* dirich_dof_ftype = "ascii";
  char dirich_val_fname[TIMEE_BUFSIZE_L];
  const char* dirich_val_fformat = "a1d";
  const char* dirich_val_ftype = "ascii";
  
  char vec_U0_fname[TIMEE_BUFSIZE_L];
  const char* vec_U0_fformat = "mtx";
  char vec_U_fpref[TIMEE_BUFSIZE_L];
  char vec_U_fsuff[TIMEE_BUFSIZE_L];
  char vec_U_fformat[TIMEE_BUFSIZE_XXS];
  const char* vec_U_ftype = "ascii";

  // -- PDE
//  Matrix<T,U> mat_M;
//  Matrix<T,U> mat_K;
//  Vector<T,U> vec_F;
//  Vector<T,U> vec_H;
  U numb_dirich;
  U* dirich_dof;
  T* dirich_val;
//  Vector<T,U> vec_U0;
//  Vector<T,U> vec_U;
  
  // -- Parareal
  Parareal<T,U> solver;
  U numb_iter;
  T res_norm;
  T res_thresh;
  T fine_time_step;
  U numb_step_per_win;
  char fine_scheme[TIMEE_BUFSIZE_XXS];
  char coarse_scheme[TIMEE_BUFSIZE_XXS];

  // -- mpi
  int rank;
  

  // ---------------------------------------------------------------------------
  // -- parameters
  // ---------------------------------------------------------------------------

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  {
    int j;

    // -- check minimum number of arguments
    if (argc < MIN_ARGC) {
      if (rank == 0) {
        puts(HELP);
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    // -- set matrix M filename
    sscanf(argv[1], "%s", mat_M_fname);

    // -- set matrix K filename
    sscanf(argv[2], "%s", mat_K_fname);

    // -- set vector F filename
    sscanf(argv[3], "%s", vec_F_fname);

    // -- set vector H filename
    sscanf(argv[4], "%s", vec_H_fname);

    // -- set Dirichlet DOF filename
    sscanf(argv[5], "%s", dirich_dof_fname);

    // -- set Dirichlet values filename
    sscanf(argv[6], "%s", dirich_val_fname);

    // -- set vector U filename prefix
    sscanf(argv[7], "%s", vec_U_fpref);

    // -- set vector U filename suffix
    sscanf(argv[8], "%s", vec_U_fsuff);

    // -- set vector U0 filename
    strcpy(vec_U0_fname, "");
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-i")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], "%s", vec_U0_fname);
    }

    // -- set vector U file format
    strcpy(vec_U_fformat, "MTX");
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-f")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], "%s", vec_U_fformat);
    }

    // -- set fine time step
    fine_time_step = 0.001;
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-t")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], FMT_T, &fine_time_step);
    }

    // -- set number of steps per window
    numb_step_per_win = 10;
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-n")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], FMT_U, &numb_step_per_win);
    }

    // -- set fine discretization scheme
    strcpy(fine_scheme, "BEULER");
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-s")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], "%s", fine_scheme);
    }

    // -- set coarse discretization scheme
    strcpy(coarse_scheme, "BEULER");
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-S")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], "%s", coarse_scheme);
    }

    // -- set residual threshold
    res_thresh = 1e-6;
    j = MIN_ARGC;
    while ((j < (argc-1)) && strcmp(argv[j], "-e")) {
      j++;
    }
    if (j < (argc-1)) {
      sscanf(argv[j+1], FMT_T, &res_thresh);
    }

    // -- print help
    j = MIN_ARGC;
    while ((j < argc) && strcmp(argv[j], "-h")) {
      j++;
    }
    if (j < argc) {
      if (rank == 0) {
        puts(HELP);
      }
      MPI_Finalize();
      return EXIT_SUCCESS;
    }
  }

  // ---------------------------------------------------------------------------
  // -- input
  // ---------------------------------------------------------------------------

  // -- read PDE from file
  
  // mass matrix
  printf("--- Reading mass matrix from %s\n", mat_M_fname);
//  mat_M.ReadFromFile(mat_M_fname, mat_M_fformat);
  
  // stiffness matrix
  printf("--- Reading stiffness matrix from %s\n", mat_K_fname);
//  mat_K.ReadFromFile(mat_K_fname, mat_K_fformat);
  
  // source vector
  printf("--- Reading source vector from %s\n", vec_F_fname);
//  vec_F.ReadFromFile(vec_F_fname, vec_F_fformat);
  
  // Neumann BC vector
  printf("--- Reading Neumann BC vector from %s\n", vec_H_fname);
//  vec_H.ReadFromFile(vec_H_fname, vec_H_fformat);
  
  // Dirichlet DOF
  printf("--- Reading Dirichlet DOF from %s\n", dirich_dof_fname);
  err = TIMEEIO<U,U>::ReadFromFile(&numb_dirich, &dirich_dof, dirich_dof_fname,
                                   dirich_dof_fformat, dirich_dof_ftype);
  if (err != TIMEE_SUCCESS) {
    printf("[%d] %s\n", rank, TIMEE_ERR_MSG);
    fflush(stdout);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  // Dirichlet values
  printf("--- Reading Dirichlet values from %s\n", dirich_val_fname);
  err = TIMEEIO<T,U>::ReadFromFile(&numb_dirich, &dirich_val, dirich_val_fname,
                                   dirich_val_fformat, dirich_val_ftype);
  if (err != TIMEE_SUCCESS) {
    printf("[%d] %s\n", rank, TIMEE_ERR_MSG);
    fflush(stdout);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  // initial value vector
  if (strcmp(vec_U0_fname, "") != 0) {
    printf("--- Reading initial value vector from %s\n", vec_U0_fname);
//    vec_U0.ReadFromFile(vec_U0_fname, vec_U0_fformat);
  }

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------
  
  printf("--- Pre-processing\n");

  // -- initialize
  
  // initial value
  if (strcmp(vec_U0_fname, "") == 0) {
//    vec_U0.Allocate(vec_F.GetSize());
//    vec_U0 = 0;
  }
  
  // Parareal
//  err = solver.Init(&vec_U, &numb_iter, &res_norm,
//                    fine_time_step, numb_step_per_win, res_thresh,
//                    MPI_COMM_WORLD, async_flag);
  if (err != TIMEE_SUCCESS) {
    puts(TIMEE_ERR_MSG);
    fflush(stdout);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  // PDE
//  err = solver.InitHeatPDE(&mat_M, &mat_K, &vec_F, &vec_H,
//                           numb_dirich, dirich_dof, dirich_val, &vec_U0,
//                           coarse_scheme, fine_scheme);
  if (err != TIMEE_SUCCESS) {
    puts(TIMEE_ERR_MSG);
    fflush(stdout);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  // -- configure file output
  err = solver.ConfigHeatPDEFileOutput(vec_U_fpref, vec_U_fsuff, vec_U_fformat,
                                       vec_U_ftype);
  if (err != TIMEE_SUCCESS) {
    puts(TIMEE_ERR_MSG);
    fflush(stdout);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------
  
  printf("--- Processing\n");

  // -- start timer
  time_begin = clock();

  // -- integrate
  err = solver.IntegrateHeatPDE();
  if (err != TIMEE_SUCCESS) {
    puts(TIMEE_ERR_MSG);
    fflush(stdout);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // -- end timer
  time_end = clock();

  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------
  
  printf("--- Post-processing\n");

  // -- finalize
  solver.FinalizeHeatPDE();
  solver.Finalize();
  delete[] dirich_dof;
  delete[] dirich_val;

  // ---------------------------------------------------------------------------
  // -- output
  // ---------------------------------------------------------------------------
  
  // -- solution
  printf("[%d] *** Number of iterations: %d  ;  Residual: %.6e\n",
         rank, numb_iter, res_norm);

  // -- processing time
  time_exec = ((double)(time_end - time_begin)) / CLOCKS_PER_SEC;
  printf("[%d] *** Processing time: %.3f sec\n", rank, time_exec);

  return EXIT_SUCCESS;
}