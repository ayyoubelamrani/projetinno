/*!
*  @file HeatPDE.hpp
*  @brief Header of class HeatPDE from package TIMEE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

#ifndef GUARD_TIMEE_HEATPDE_HPP_
#define GUARD_TIMEE_HEATPDE_HPP_

// C++ headers

// MRG TIMEE headers
#include <TIMEE/TIMEEError.hpp>
#include <TIMEE/TIMEEConst.hpp>

// MRG third-party headers
//#include <Alinea/Matrix.hpp>
//#include <Alinea/Vector.hpp>
//#include <Alinea/Solver.hpp>

// Other third-party headers

//! @class HeatPDE
//! @brief Integrates FEM-discrete time-dependent PDEs of the form
//!        du/dt - Lu = f in Omega
//!                 u = g on Gamma1
//!             du/dn = h on Gamma2
//! @brief Schemes: Backward Euler, Trapezoidal rule.
template <typename T, typename U>
class HeatPDE {

  // ---------------------------------------------------------------------------
  // -- Input attributes
  // ---------------------------------------------------------------------------

  private:

    //! Mass matrix.
//    Matrix<T,U>* m_mat_M;

    //! Stiffness matrix.
//    Matrix<T,U>* m_mat_K;

    //! Source vector.
//    Vector<T,U>* m_vec_F;

    //! Neumann boundary conditions vector.
//    Vector<T,U>* m_vec_H; 

    //! Number of Dirichlet boundary DOF.
    U m_numb_dirich;

    //! List of Dirichlet boundary DOF.
    U* m_dirich_dof;

    //! Constraint value of each Dirichlet boundary DOF.
    T* m_dirich_val;

    //! Initial value vector.
//    Vector<T,U>* m_vec_U0;

    //! Time step.
    T m_time_step;

    //! Number of time steps.
    U m_numb_step;
    
    //! Time discretization scheme (BEULER [default], TRULE)
    char m_discret_scheme[TIMEE_BUFSIZE_XXS];

  // ---------------------------------------------------------------------------
  // -- Output attributes
  // ---------------------------------------------------------------------------

  private:

    //! Solution vector at the final time step.
//    Vector<T,U>* m_vec_U;

  // ---------------------------------------------------------------------------
  // -- Configuration attributes
  // ---------------------------------------------------------------------------

  private:

    //! Solution vector filename prefix for all time steps (1023 car. max.).
    char m_vec_U_fpref[TIMEE_BUFSIZE_L];

    //! Solution vector filename suffix for all time steps (1023 car. max.).
    char m_vec_U_fsuff[TIMEE_BUFSIZE_L];

    //! Solution vector file format (MTX [default], A1D, ENSI).
    char m_vec_U_fformat[TIMEE_BUFSIZE_XXS];

    //! Solution vector file type (ASCII [default]).
    char m_vec_U_ftype[TIMEE_BUFSIZE_XXS];

  // ---------------------------------------------------------------------------
  // -- Local attributes
  // ---------------------------------------------------------------------------

  private:

    //! Linear system for each time step resolution.
    //! A Un+1 = B Un + C
//    Matrix<T,U> m_mat_A;
//    Matrix<T,U> m_mat_B;
//    Vector<T,U> m_vec_C;
    
    //! Solver.
//    Solver<T,U> m_solver;
    
    //! File output activation flag.
    short m_flag_fout;

  // ---------------------------------------------------------------------------
  // -- Pre/Post-processing methods
  // ---------------------------------------------------------------------------

  public:

    //! @brief Initializes attributes.
    HeatPDE ( void ) ;

    //! @brief Initializes attributes.
    int Init (
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
        const char* discret_scheme = "BEULER") ;

    //! @brief Restores attributes to default values.
    int Finalize ( void ) ;

    //! @brief Configures file output for all time steps solutions.
    //! @remarks Passing NULL pointer disables file output.
    int ConfigFileOutput (
        const char* vec_U_fpref,
        const char* vec_U_fsuff,
        const char* vec_U_fformat = "MTX",
        const char* vec_U_ftype = "ASCII" ) ;
    
  private:

    //! @brief Initializes and enables Backward Euler integration scheme.
    int InitBEuler ( void ) ;

    //! @brief Initializes and enables Trapezoidal rule integration scheme.
    int InitTRule ( void ) ;

    //! @brief Applies Dirichlet boundary conditions.
    int InitDirichlet ( void ) ;

    //! @brief Factorizes the linear system matrix.
    int InitLDLt ( void ) ;

  // ---------------------------------------------------------------------------
  // -- Processing methods
  // ---------------------------------------------------------------------------

  public:

    //! @brief Computes approximated solution at each time step.
    int Integrate ( void ) ;

}; // class HeatPDE {

#endif // #ifndef GUARD_TIMEE_HEATPDE_HPP_
