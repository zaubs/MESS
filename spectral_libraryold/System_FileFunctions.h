//#########################################################################
//#########################################################################
//
//  System File Functions for Windows or Linux
//
//#########################################################################
//#########################################################################

#ifndef _H_GUARD_SYSTEM_FILEFUNCTIONS
#define _H_GUARD_SYSTEM_FILEFUNCTIONS

#ifdef _WIN32 /************* WINDOWS ******************************/

  #include <windows.h>    // string fncts, malloc/free, system

#else /********************** LINUX *******************************/

  #include <dirent.h>
  #include <pthread.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>

  extern unsigned int sleep( unsigned int __seconds );

#endif /***********************************************************/


int  GenerateImageryFileListing( char *folder_pathname, char *listing_pathname );
void Delay_msec( int milliseconds );
void SplitPathname( char* fullpathname, char* folder_pathname, char* filename );
void BuildPathname( char* datafolder_pathname, char* imagery_filename, char* imagery_pathname );
int  GetRandom();


#endif /* _H_GUARD_SYSTEM_FILEFUNCTIONS */
