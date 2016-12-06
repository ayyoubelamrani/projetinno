/*!
*  @file TIMEEError.cpp
*  @brief Source of specialized error for TIMEE.
*  @author Magoules Research Group ( MRG )
*  @author GBG, CAAK, MF
*  @date 2016-07-16, 2016-07-27
*  @version 4.1
*  @post 1996, 2005, 2007, 2012
*  @remarks 0000-00000000-XXXXX-0-XX
*/

// C++ packages

// MRG TIMEE packages
#include <TIMEE/TIMEEError.hpp>

// MRG third-party packages

// -- last error
int TIMEE_ERR_CODE = TIMEE_SUCCESS;
// -- global error message string (TIMEE)
char TIMEE_ERR_MSG[TIMEE_ERR_MSG_SIZE];
