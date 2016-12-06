/*!
*  @file HeatPDE.cpp
*  @brief Source of class HeatPDE from package TIMEE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

// C++ headers
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cstdio>

// MRG TIMEE headers
#include <TIMEE/HeatPDE.hpp>
#include <TIMEE/TIMEEIO.hpp>

// MRG third-party headers
//#include <Alinea/Blap.hpp>

// Including Eigen library

#include </Users/TM/Travail/Programming/C/libraries/Eigen/Dense>

/* __________________________________________________________________________ */

  // ---------------------------------------------------------------------------
  // -- Pre/Post-processing methods
  // ---------------------------------------------------------------------------

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
HeatPDE<T,U>::HeatPDE ( void ) {

  // -- PDE
//  m_mat_M = NULL;
//  m_mat_K = NULL;
//  m_vec_F = NULL;
//  m_vec_H = NULL;
  m_numb_dirich = 0;
  m_dirich_dof = NULL;
  m_dirich_val = NULL;
//  m_vec_U0 = NULL;
//  m_vec_U = NULL;

  // -- discretization
  m_numb_step = 0;
  strcpy(m_discret_scheme, "BEULER");

  // -- file I/O
  m_flag_fout = 0;

}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int HeatPDE<T,U>::Init (
//        Vector<T,U>* vec_U,
//        Matrix<T,U>* mat_M,
//        Matrix<T,U>* mat_K,
//        Vector<T,U>* vec_F,
//        Vector<T,U>* vec_H,
        U numb_dirich,
        U* dirich_dof,
        T* dirich_val,
//        Vector<T,U>* vec_U0,
        T time_step,
        U numb_step,
        const char* discret_scheme ) {

  // -- PDE
//  m_vec_U = vec_U;
//  m_mat_M = mat_M;
//  m_mat_K = mat_K;
//  m_vec_F = vec_F;
//  m_vec_H = vec_H;
  m_numb_dirich = numb_dirich;
  m_dirich_dof = dirich_dof;
  m_dirich_val = dirich_val;
//  m_vec_U0 = vec_U0;

  // -- discretization
  m_time_step = time_step;
  m_numb_step = numb_step;
  if (discret_scheme != NULL) {
    strcpy(m_discret_scheme, discret_scheme);
  }

  // -- linear system

  // uppercase keywords
  for (unsigned int i = 0; i < strlen(m_discret_scheme); i++) {
    m_discret_scheme[i] = toupper(m_discret_scheme[i]);
  }

  // apply discretization scheme
  if (strcmp(m_discret_scheme, "BEULER") == 0) {
    this->InitBEuler();
  }
  else if (strcmp(m_discret_scheme, "TRULE") == 0) {
    this->InitTRule();
  }
  else {
    TIMEE_ERR_CODE = TIMEE_ERR_NOT_IMPLEMENTED;
    sprintf(TIMEE_ERR_MSG,
            "*** [TIMEE] Error %d: %s scheme not implemented.\n",
            TIMEE_ERR_CODE, m_discret_scheme);
    return TIMEE_ERR_CODE;
  }

  // apply Dirichlet boundary conditions
  this->InitDirichlet();

  // factorize
  this->InitLDLt();

  // -- solver
//  m_solver.SetSolver("LDLT");
//  m_solver.Setup(m_mat_A);

  return TIMEE_SUCCESS;;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int HeatPDE<T,U>::Finalize ( void ) {

  // -- PDE
//  m_mat_M = NULL;
//  m_mat_K = NULL;
//  m_vec_F = NULL;
//  m_vec_H = NULL;
  m_numb_dirich = 0;
  m_dirich_dof = NULL;
  m_dirich_val = NULL;
//  m_vec_U0 = NULL;
  m_numb_step = 0;
//  m_vec_U = NULL;

  // -- file I/O
  m_flag_fout = 0;

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int HeatPDE<T,U>::ConfigFileOutput (
        const char* vec_U_fpref,
        const char* vec_U_fsuff,
        const char* vec_U_fformat,
        const char* vec_U_ftype ) {

  if ((vec_U_fpref == NULL) || (vec_U_fsuff == NULL)
      || (vec_U_fformat == NULL) || (vec_U_ftype == NULL)) {
    m_flag_fout = 0;
  }
  else {
    // -- set attributes
    strcpy(m_vec_U_fpref, vec_U_fpref);
    strcpy(m_vec_U_fsuff, vec_U_fsuff);
    strcpy(m_vec_U_fformat, vec_U_fformat);
    strcpy(m_vec_U_ftype, vec_U_ftype);
    m_flag_fout = 1;

    // -- uppercase keywords
    for (unsigned int i = 0; i < strlen(m_vec_U_fformat); i++) {
      m_vec_U_fformat[i] = toupper(m_vec_U_fformat[i]);
    }
    for (unsigned int i = 0; i < strlen(m_vec_U_ftype); i++) {
      m_vec_U_ftype[i] = toupper(m_vec_U_ftype[i]);
    }
  }

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal Backward Euler: (M + dt K)Un+1 = M Un + dt F + H
template <typename T, typename U>
int HeatPDE<T,U>::InitBEuler ( void ) {

  // -- set A = M + dt K
//  m_mat_A = (*m_mat_M) + (*m_mat_K) * m_time_step;

  // -- set B = M
//  m_mat_B = (*m_mat_M);

  // -- set C = dt F + H
//  m_vec_C = (*m_vec_H);
//  Blap1::Saxpy<T,U>(m_time_step, (*m_vec_F), m_vec_C);

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal Trapezoidal rule: (M + dt/2 K)Un+1 = (M - dt/2 K)Un + dt F + H
template <typename T, typename U>
int HeatPDE<T,U>::InitTRule ( void ) {

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal
//!    | Aii 0 |        | Bii Bif |      | Ci - Aif * d |
//!    |  0  1 | Un+1 = |  0   0  | Un + |       d      |
template <typename T, typename U>
int HeatPDE<T,U>::InitDirichlet ( void ) {

  // -- get number of DOF
//  U numb_dof = m_vec_C.GetSize();
  U numb_dof;

  // -- build constraint values and mask on all DOF

  // initializei
  int* dof_mask = new int[numb_dof];
  T* dof_val = new T[numb_dof];
  for (U i = 0; i < numb_dof; i++) {
    dof_mask[i] = 0;
  }

  // build
  for (U i = 0; i < m_numb_dirich; i++) {
    dof_mask[m_dirich_dof[i]] = 1;
    dof_val[m_dirich_dof[i]] = m_dirich_val[i];
  }

  // -- update linear system
  //    | Aii Aif |        | Bii Bif |      | Ci |
  //    | Afi Aff | Un+1 = | Bfi Bff | Un + | Cf |
  for (U i = 0; i < numb_dof; i++) {
    if (dof_mask[i] == 1) {
      // set Aif = 0
      U jj;
//      for (U j = m_mat_A.GetRowsIndices()[i]; j < m_mat_A.GetRowsIndices()[i+1];
//           j++) {
//        jj = m_mat_A.GetColumnsNumb()[j];
        if (dof_mask[jj] == 1) {
//          m_mat_A.GetCoef()[j] = T(0);
        }
//      }
    }
  }

  // -- finalize
  delete[] dof_mask;
  delete[] dof_val;

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal Factorizes the linear system matrix given by
//!           A Un+1 = B Un + C
template <typename T, typename U>
int HeatPDE<T,U>::InitLDLt ( void ) {

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

  // ---------------------------------------------------------------------------
  // -- Processing methods
  // ---------------------------------------------------------------------------

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int HeatPDE<T,U>::Integrate ( void ) {

  // -- init

  // error code
  int err = TIMEE_SUCCESS;

  // solution
//  (*m_vec_U) = (*m_vec_U0);

  // -- integrate without file output
  if (!m_flag_fout) {
  }

  // -- integrate with file output
  else {
  }

  return err;
}

/* __________________________________________________________________________ */

template class HeatPDE<double,int>;
