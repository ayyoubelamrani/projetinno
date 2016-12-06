/*!
*  @file Parareal.hpp
*  @brief Header of class Parareal from package TIMEE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

#ifndef GUARD_TIMEE_PARAREAL_HPP_
#define GUARD_TIMEE_PARAREAL_HPP_

// C++ headers

// MRG TIMEE headers
#include <TIMEE/HeatPDE.hpp>

// MRG third-party headers
#include <JACKComm.hpp>

// Other third-party headers

//! @class Parareal
//! @brief Integrates time-dependent equations using the Parareal scheme.
//! @brief This class handles one time window and communicates with neighbors.
template <typename T, typename U>
class Parareal {

  // ---------------------------------------------------------------------------
  // -- Input attributes
  // ---------------------------------------------------------------------------

  private:

    //! Fine time step.
    T m_time_step;

    //! Number of time steps in the time window.
    U m_numb_step;

    //! Residual threshold (default: 1e-6).
    T m_res_thresh;

  // ---------------------------------------------------------------------------
  // -- Output attributes
  // ---------------------------------------------------------------------------

  private:

    //! Solution vector at the end of the time window.
//    Vector<T,U>* m_vec_U;
    
    //! Number of iterations.
    U* m_numb_iter;
    
    //! Residual norm.
    T* m_res_norm;

  // ---------------------------------------------------------------------------
  // -- Local attributes
  // ---------------------------------------------------------------------------

  private:

    //! Coarse integrator.
    HeatPDE<T,U> m_coarse_heatpde;
    
    //! Coarse solution vector.
//    Vector<T,U> m_coarse_vec_U;

    //! Fine integrator.
    HeatPDE<T,U> m_fine_heatpde;
    
    //! Fine solution vector.
//    Vector<T,U> m_fine_vec_U;

    //! Initial value vector for the window.
//    Vector<T,U> m_vec_U0;
    
    //! Communication tool.
    //! @see JACKComm::Init.
    JACKComm<T,U> m_jack_comm;
    MPI_Comm m_mpi_comm;
    U m_rank;
    U m_numb_sneighb;
    U m_numb_rneighb;
    U* m_sneighb_rank;
    U* m_rneighb_rank;
    U* m_sbuf_size;
    U* m_rbuf_size;
    T** m_send_buf;
    T** m_recv_buf;
    T* m_res_buf;

  // ---------------------------------------------------------------------------
  // -- Pre/Post-processing methods
  // ---------------------------------------------------------------------------

  public:

    //! @brief Initializes attributes.
    Parareal ( void ) ;

    //! @brief Initializes attributes.
    int Init (
//        Vector<T,U>* vec_U,
        U* numb_iter,
        T* res_norm,
        T time_step,
        U numb_step,
        T res_thresh = 1e-6,
        MPI_Comm mpi_comm = MPI_COMM_WORLD ) ;

    //! @brief Initializes coarse and fine HeatPDE integrators.
    //! @see HeatPDE::Init
    int InitHeatPDE (
//        Matrix<T,U>* mat_M,
//        Matrix<T,U>* mat_K,
//        Vector<T,U>* vec_F,
//        Vector<T,U>* vec_H,
        U numb_dirich,
        U* dirich_dof,
        T* dirich_val,
//        Vector<T,U>* vec_U0,
        const char* coarse_scheme = "BEULER",
        const char* fine_scheme = "BEULER" ) ;

    //! @brief Finalizes coarse and fine HeatPDE integrators.
    //! @see HeatPDE::Finalize
    int FinalizeHeatPDE ( void ) ;

    //! @brief Restores attributes to default values.
    int Finalize ( void ) ;

    //! @brief Configures file output for fine HeatPDE integrator.
    //! @see HeatPDE::ConfigFileOutput
    int ConfigHeatPDEFileOutput (
        const char* vec_U_fpref,
        const char* vec_U_fsuff,
        const char* vec_U_fformat = "MTX",
        const char* vec_U_ftype = "ASCII" ) ;

  // ---------------------------------------------------------------------------
  // -- Processing methods
  // ---------------------------------------------------------------------------

  public:

    //! @brief Computes HeatPDE approximated solution at each time step.
    int IntegrateHeatPDE ( void ) ;

}; // class Parareal {

#endif // #ifndef GUARD_TIMEE_PARAREAL_HPP_