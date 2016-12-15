/*!
*  @file Parareal.cpp
*  @brief Source of class Parareal from package TIMEE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

// C++ headers
#include <cstdlib>

// MRG TIMEE headers
#include <TIMEE/Parareal.hpp>

// MRG third-party headers
//#include <Alinea/Blap.hpp>

/* __________________________________________________________________________ */

  // ---------------------------------------------------------------------------
  // -- Pre/Post-processing methods
  // ---------------------------------------------------------------------------

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
Parareal<T,U>::Parareal ( void ) {

  // -- solution
//  m_vec_U = NULL;
  m_numb_iter = NULL;
  m_res_norm = NULL;
  
  // -- communication
  m_mpi_comm = MPI_COMM_WORLD;
  m_sneighb_rank = NULL;
  m_rneighb_rank = NULL;
  m_sbuf_size = NULL;
  m_rbuf_size = NULL;
  m_send_buf = NULL;
  m_recv_buf = NULL;
  m_res_buf = NULL;
  
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int Parareal<T,U>::Init (
//        Vector<T,U>* vec_U,
        U* numb_iter,
        T* res_norm,
        T time_step,
        U numb_step,
        T res_thresh,
        MPI_Comm mpi_comm ) {

  // -- solution
//  m_vec_U = vec_U;
  m_numb_iter = numb_iter;
  m_res_norm = res_norm;
  m_res_thresh = res_thresh;
  
  // -- discretization
  m_time_step = time_step;
  m_numb_step = numb_step;
  
  // -- communication
  
  // network
  U size;
  m_numb_sneighb = 1;
  m_numb_rneighb = 1;
  m_sneighb_rank = new U[1];
  m_rneighb_rank = new U[1];
  m_mpi_comm = mpi_comm;
  MPI_Comm_rank(mpi_comm, &m_rank);
  MPI_Comm_size(mpi_comm, &size);
  m_sneighb_rank[0] = m_rank+1;
  m_rneighb_rank[0] = m_rank-1;
  if (m_rank == (size-1)) {
    m_numb_sneighb = 0;
  }
  if (m_rank == 0) {
    m_numb_rneighb = 0;
  }
  
  // buffers
  m_sbuf_size = new U[1];
  m_rbuf_size = new U[1];
//  m_sbuf_size[0] = (*m_vec_U).GetSize();
//  m_rbuf_size[0] = (*m_vec_U).GetSize();
  m_send_buf = new T*[1];
  m_recv_buf = new T*[1];
  m_res_buf = new T[1];
//  m_send_buf[0] = (*m_vec_U).GetCoef();

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int Parareal<T,U>::InitHeatPDE (
//        Matrix<T,U>* mat_M,
//        Matrix<T,U>* mat_K,
//        Vector<T,U>* vec_F,
//        Vector<T,U>* vec_H,
        U numb_dirich,
        U* dirich_dof,
        T* dirich_val,
//        Vector<T,U>* vec_U0,
        const char* coarse_scheme,
        const char* fine_scheme ) {

  // -- initial value
//  m_vec_U0 = (*vec_U0);
  
  // -- coarse integrator
//  m_coarse_vec_U.Allocate((*m_vec_U).GetSize());
//  m_coarse_heatpde.Init(&m_coarse_vec_U, mat_M, mat_K, vec_F, vec_H,
//                        numb_dirich, dirich_dof, dirich_val, &m_vec_U0,
//                        m_numb_step * m_time_step, 1,
//                        coarse_scheme);

  // -- fine integrator
//  m_fine_vec_U.Allocate((*m_vec_U).GetSize());
//  m_fine_heatpde.Init(&m_fine_vec_U, mat_M, mat_K, vec_F, vec_H,
//                      numb_dirich, dirich_dof, dirich_val, &m_vec_U0,
//                      m_time_step, m_numb_step,
//                      fine_scheme);

  // -- communicator
//  m_recv_buf[0] = m_vec_U0.GetCoef();
/*  m_jack_comm.Init(m_recv_buf, m_res_norm,
                   m_sbuf_size, m_rbuf_size, m_send_buf, 1, m_res_buf,
                   m_numb_sneighb, m_numb_rneighb,
                   m_sneighb_rank, m_rneighb_rank,
                   0,
                   m_mpi_comm);*/

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int Parareal<T,U>::FinalizeHeatPDE ( void ) {

  // -- communicator
//  m_jack_comm.Finalize() ;

  // -- integrators
  m_fine_heatpde.Finalize();
  m_coarse_heatpde.Finalize();
//  m_fine_vec_U.Deallocate();
//  m_coarse_vec_U.Deallocate();

  // -- initial value
//  m_vec_U0.Deallocate();

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int Parareal<T,U>::Finalize ( void ) {
  
  // -- communication
  
  // network
  delete[] m_sneighb_rank;
  m_sneighb_rank = NULL;
  delete[] m_rneighb_rank;
  m_rneighb_rank = NULL;
  
  // buffers
  delete[] m_sbuf_size;
  m_sbuf_size = NULL;
  delete[] m_rbuf_size;
  m_rbuf_size = NULL;
  delete[] m_send_buf;
  m_send_buf = NULL;
  delete[] m_recv_buf;
  m_recv_buf = NULL;
  delete[] m_res_buf;
  m_res_buf = NULL;

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int Parareal<T,U>::ConfigHeatPDEFileOutput (
        const char* vec_U_fpref,
        const char* vec_U_fsuff,
        const char* vec_U_fformat,
        const char* vec_U_ftype ) {

  m_fine_heatpde.ConfigFileOutput(vec_U_fpref, vec_U_fsuff,
                                  vec_U_fformat, vec_U_ftype);

  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

  // ---------------------------------------------------------------------------
  // -- Processing methods
  // ---------------------------------------------------------------------------

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int Parareal<T,U>::IntegrateHeatPDE ( void ) {
  
  // -- init
  
  // previous iteration coarse prediction
//  Vector<T,U> coarse_vec_U_prev;
  
  // previous iteration solution: Un+1<k>
//  Vector<T,U> vec_U_prev;

  // local residual vector: Un+1<k+1> - Un+1<k>
//  Vector<T,U> vec_local_res;
  
  // -- iterate
  (*m_numb_iter) = 0;
  (*m_res_norm) = m_res_thresh;
  while (((*m_res_norm) >= m_res_thresh) && ((*m_numb_iter) < m_rank)) {
    
    // -- recv Un<k+1>
//    m_jack_comm.Recv();
    
    // -- send Un+1<k+1>
//    m_jack_comm.Send();
    
    // -- eval |Un+1<k+1> - Un+1<k>|
//    (*m_res_buf) = vec_local_res.NormL2();
//    m_jack_comm.UpdateResidual();
    
    // -- k := k+1
    (*m_numb_iter)++;
  }
  
  return TIMEE_SUCCESS;
}

/* __________________________________________________________________________ */

template class Parareal<double,int>;