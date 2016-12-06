/*!
*  @file TIMEEIO.hpp
*  @brief Header of class IO from package TIMEE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

#ifndef GUARD_TIMEE_IO_HPP_
#define GUARD_TIMEE_IO_HPP_

// C++ headers

// MRG TIMEE headers

// MRG third-party headers

//! @class TIMEEIO
//! @brief Data Input/Output utilities.
template <typename T, typename U>
class TIMEEIO {

  // ---------------------------------------------------------------------------
  // -- File IO class methods
  // ---------------------------------------------------------------------------

  public:

    //! @brief Reads a 1D array from a file.
    //! @param array_size [out] size of the array
    //! @param array_buf [out] reference to the array buffer
    //! @param file_name [in] full name of the file
    //! @param file_format [in] format of the contents (a1d)
    //! @param file_type [in] type of the storage (ascii)
    static int ReadFromFile (
        U* array_size,
        T** array_buf,
        const char* file_name,
        const char* file_format="a1d",
        const char* file_type="ascii" ) ;

  private:

    //! @brief Reads a 1D array from an MRG A1D ascii file.
    //! @see ReadFromFile
    static int ReadFromFileA1DAscii (
        U* array_size,
        T** array_buf,
        const char* file_name ) ;

  public:

    //! @brief Writes a 1D array into a file.
    //! @see ReadFromFile
    static int WriteToFile (
        U array_size,
        T* array_buf,
        const char* file_name,
        const char* file_format="a1d",
        const char* file_type="ascii" ) ;

  private:

    //! @brief Writes a 1D array into an MRG A1D ascii file.
    //! @see ReadFromFile
    static int WriteToFileA1DAscii (
        U array_size,
        T* array_buf,
        const char* file_name ) ;

    //! @brief Writes a 1D array into an Ensight scalar per-node ascii file.
    //! @see ReadFromFile
    static int WriteToFileEnsiAscii (
        U array_size,
        T* array_buf,
        const char* file_name ) ;

}; // class TIMEEIO {

#endif // #ifndef GUARD_TIMEE_IO_HPP_
