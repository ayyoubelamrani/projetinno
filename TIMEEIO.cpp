/*!
*  @file TIMEEIO.cpp
*  @brief Source of class IO from package TIMEE.
*  @author CentraleSupelec, MICS, Pr. Magoules HPC Research Group (MRG).
*  @author G. G.-Benissan
*  @date 2016-10-13, 2016-10-13
*/

// C++ headers
#include <cstring>
#include <cctype>
#include <cstdio>
#include <typeinfo>
#include <cstdlib>

// MRG TIMEE headers
#include <TIMEE/TIMEEIO.hpp>
#include <TIMEE/TIMEEError.hpp>
#include <TIMEE/TIMEEConst.hpp>

// MRG third-party headers

/* __________________________________________________________________________ */

  // ---------------------------------------------------------------------------
  // -- File IO class methods
  // ---------------------------------------------------------------------------

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int TIMEEIO<T,U>::ReadFromFile (
        U* array_size,
        T** array_buf,
        const char* file_name,
        const char* file_format,
        const char* file_type ) {

  // -- return code
  int err = TIMEE_SUCCESS;

  // -- keywords
  char* fformat = new char[strlen(file_format)+1];
  char* ftype = new char[strlen(file_type)+1];

  // -- init
  strcpy(fformat, file_format);
  strcpy(ftype, file_type);

  // -- uppercase keywords
  for (unsigned int i = 0; i < strlen(fformat); i++) {
    fformat[i] = toupper(fformat[i]);
  }
  for (unsigned int i = 0; i < strlen(ftype); i++) {
    ftype[i] = toupper(ftype[i]);
  }

  // -- read according to keywords
  if (strcmp(fformat, "A1D") == 0) {
    if (strcmp(ftype, "ASCII") == 0) {
      err = TIMEEIO<T,U>::ReadFromFileA1DAscii(array_size, array_buf,
                                                     file_name);
    }
    else {
      err = TIMEE_ERR_NOT_IMPLEMENTED;
      sprintf(TIMEE_ERR_MSG,
              "*** [TIMEE] Error %d: %s %s format not implemented.\n",
              err, ftype, fformat);
      TIMEE_ERR_CODE = err;
    }
  }
  else {
    err = TIMEE_ERR_NOT_IMPLEMENTED;
    sprintf(TIMEE_ERR_MSG,
            "*** [TIMEE] Error %d: %s %s format not implemented.\n",
            err, ftype, fformat);
    TIMEE_ERR_CODE = err;
  }

  // -- finalize
  delete[] fformat;
  delete[] ftype;

  return err;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int TIMEEIO<T,U>::ReadFromFileA1DAscii (
        U* array_size,
        T** array_buf,
        const char* file_name ) {

  // -- return code
  int err = TIMEE_SUCCESS;

  // -- file
  FILE* p_file;
  char header[TIMEE_BUFSIZE_S];

  //-- string format specifier
  char fmt_u[TIMEE_BUFSIZE_XXS];
  char fmt_t[TIMEE_BUFSIZE_XXS];

  // -- init
  if (typeid(U) == typeid(long)) {
    strcpy(fmt_u, "%ld");
  }
  else {
    strcpy(fmt_u, "%d");
  }
  if (typeid(T) == typeid(double)) {
    strcpy(fmt_t, "%lf");
  }
  else if (typeid(T) == typeid(float)) {
    strcpy(fmt_t, "%f");
  }
  else if (typeid(T) == typeid(long)) {
    strcpy(fmt_t, "%ld");
  }
  else {
    strcpy(fmt_t, "%d");
  }

  // -- open file
  p_file = fopen(file_name, "r");
  if (p_file == NULL) {
    err = TIMEE_ERR_IO_FILE_OPEN;
    sprintf(TIMEE_ERR_MSG,
            "*** [TIMEE] Error %d: can not open %s.\n",
            err, file_name);
    TIMEE_ERR_CODE = err;
  }
  else {
    // -- skip comments
    fscanf(p_file, "%s", header);
    while (!feof(p_file) && (header[0] == '%')) {
      fscanf(p_file, "%[^\n]", header);
      fscanf(p_file, "%s", header);
    }

    // -- read size
    sscanf(header, fmt_u, array_size);

    // -- read array
    (*array_buf) = new T[(*array_size)];
    for (U i = 0; i < (*array_size); i++) {
      fscanf(p_file, fmt_t, &(*array_buf)[i]);
    }

    // -- close file
    fclose(p_file);
  }

  return err;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int TIMEEIO<T,U>::WriteToFile (
        U array_size,
        T* array_buf,
        const char* file_name,
        const char* file_format,
        const char* file_type ) {

  // -- return code
  int err = TIMEE_SUCCESS;

  // -- keywords
  char* fformat = new char[strlen(file_format)+1];
  char* ftype = new char[strlen(file_type)+1];

  // -- init
  strcpy(fformat, file_format);
  strcpy(ftype, file_type);

  // -- uppercase keywords
  for (unsigned int i = 0; i < strlen(fformat); i++) {
    fformat[i] = toupper(fformat[i]);
  }
  for (unsigned int i = 0; i < strlen(ftype); i++) {
    ftype[i] = toupper(ftype[i]);
  }

  // -- read according to keywords
  if (strcmp(fformat, "A1D") == 0) {
    if (strcmp(ftype, "ASCII") == 0) {
      err = TIMEEIO<T,U>::WriteToFileA1DAscii(array_size, array_buf,
                                                  file_name);
    }
    else {
      err = TIMEE_ERR_NOT_IMPLEMENTED;
      sprintf(TIMEE_ERR_MSG,
              "*** [TIMEE] Error %d: %s %s format not implemented.\n",
              err, ftype, fformat);
      TIMEE_ERR_CODE = err;
    }
  }
  else if (strcmp(fformat, "ENSI") == 0) {
    if (strcmp(ftype, "ASCII") == 0) {
      err = TIMEEIO<T,U>::WriteToFileEnsiAscii(array_size, array_buf,
                                                  file_name);
    }
    else {
      err = TIMEE_ERR_NOT_IMPLEMENTED;
      sprintf(TIMEE_ERR_MSG,
              "*** [TIMEE] Error %d: %s %s format not implemented.\n",
              err, ftype, fformat);
      TIMEE_ERR_CODE = err;
    }
  }
  else {
    err = TIMEE_ERR_NOT_IMPLEMENTED;
    sprintf(TIMEE_ERR_MSG,
            "*** [TIMEE] Error %d: %s %s format not implemented.\n",
            err, ftype, fformat);
    TIMEE_ERR_CODE = err;
  }

  // -- finalize
  delete[] fformat;
  delete[] ftype;

  return err;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int TIMEEIO<T,U>::WriteToFileA1DAscii (
        U array_size,
        T* array_buf,
        const char* file_name ) {

  // -- return code
  int err = TIMEE_SUCCESS;

  // -- file
  FILE* p_file;
  char header[TIMEE_BUFSIZE_S];

  //-- string format specifier
  char fmt_u[TIMEE_BUFSIZE_XXS];
  char fmt_t[TIMEE_BUFSIZE_XXS] = "%.6f\n";

  // -- init
  strcpy(fmt_u, "%d\n");
  if (typeid(T) == typeid(double)) {
    strcpy(fmt_t, "%.12f\n");
  }
  else if (typeid(T) == typeid(float)) {
    strcpy(fmt_t, "%.6f\n");
  }
  else {
    strcpy(fmt_t, "%d\n");
  }

  // -- open file
  p_file = fopen(file_name, "w");
  if (p_file == NULL) {
    err = TIMEE_ERR_IO_FILE_OPEN;
    sprintf(TIMEE_ERR_MSG,
            "*** [TIMEE] Error %d: can not open %s.\n",
            err, file_name);
    TIMEE_ERR_CODE = err;
  }
  else {
    // -- write comments
    sprintf(header, "%% Pr. Magoules HPC Research Group (MRG) A1D ascii file.");
    fprintf(p_file, "%s\n", header);

    // -- write size
    fprintf(p_file, fmt_u, array_size);

    // -- write array
    for (U i = 0; i < array_size; i++) {
      fprintf(p_file, fmt_t, array_buf[i]);
    }

    // -- close file
    fclose(p_file);
  }

  return err;
}

/* __________________________________________________________________________ */

//! @internal See header file.
template <typename T, typename U>
int TIMEEIO<T,U>::WriteToFileEnsiAscii (
        U array_size,
        T* array_buf,
        const char* file_name ) {

  // -- return code
  int err = TIMEE_SUCCESS;

  // -- file
  FILE* p_file;
  char header[TIMEE_BUFSIZE_S];

  //-- string format specifier
  char fmt_u[TIMEE_BUFSIZE_XXS] = "%10d\n";
  char fmt_t[TIMEE_BUFSIZE_XXS] = "%12.5e\n";

  // -- init
  if ((typeid(T) == typeid(int)) || (typeid(T) == typeid(long))) {
    strcpy(fmt_t, "%10d\n");
  }

  // -- open file
  p_file = fopen(file_name, "w");
  if (p_file == NULL) {
    err = TIMEE_ERR_IO_FILE_OPEN;
    sprintf(TIMEE_ERR_MSG,
            "*** [TIMEE] Error %d: can not open %s.\n",
            err, file_name);
    TIMEE_ERR_CODE = err;
  }
  else {
    // -- write comments
    sprintf(header, "Alya Ensight Gold - Scalar per-node variables file");
    sprintf(header, "%s --- Generated by MRG TIMEE lib.", header);
    fprintf(p_file, "%s\n", header);

    // -- write part number (only one part)
    fprintf(p_file, "part\n");
    fprintf(p_file, fmt_u, 1);

    // -- write array
    fprintf(p_file, "coordinates\n");
    for (U i = 0; i < array_size; i++) {
      fprintf(p_file, fmt_t, array_buf[i]);
    }

    // -- close file
    fclose(p_file);
  }

  return err;
}

/* __________________________________________________________________________ */

template class TIMEEIO<double,int>;
template class TIMEEIO<int,int>;