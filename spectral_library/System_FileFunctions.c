//#########################################################################
//#########################################################################
//
//  System File Functions for Windows or Linux
//
//#########################################################################
//#########################################################################

#pragma warning(disable: 4996)  // disable warning on strcpy, fopen, ...
#include "System_FileFunctions.h"

#ifndef _WIN32 /************* LINUX ******************************/
    /*
     * System_FileFunctions requires some global state to track init that only runs
     * once and thread local storage for random seeds so we can't include this in a
     * single header library
     */
    static pthread_once_t getrandom_init_once = PTHREAD_ONCE_INIT;
    static pthread_key_t getrandom_seed_key;

#endif /************* End of WINDOWS or LINUX *********************************/


//################################################################################
//################################################################################
//
//    Get the directory listing of all the *.fits or *.vid files in the user
//    specified folder and pipe it to a text file "listing_pathname".
//
//################################################################################
//################################################################################

int  GenerateImageryFileListing( char *folder_pathname, char *listing_pathname )
{

#ifdef _WIN32 /************* WINDOWS *******************************************/

    char    windows_system_command[512];

    strcpy( windows_system_command, "dir /b \"" );

    strcat( windows_system_command, folder_pathname );

    strcat( windows_system_command, "*.fits\" > \"" );   //... listing of FITS files

    strcat( windows_system_command, listing_pathname );

    strcat( windows_system_command, "\"" );

    system( windows_system_command );

    return(0);


#else /********************** LINUX *******************************************/

    struct dirent *entry;
    DIR           *directoryfolder;
    FILE          *filelisting;


    if( ( filelisting = fopen( listing_pathname, "w" ) ) == NULL )  {
        fprintf( stdout, "ERROR ===> FITS file directory listing %s cannot be open for write\n\n", listing_pathname );
        fprintf( stderr, "ERROR ===> FITS file directory listing %s cannot be open for write\n\n", listing_pathname );
        return(1);
    }

    if( ( directoryfolder = opendir( folder_pathname ) ) == NULL )  {
        fprintf( stdout, "ERROR ===> FITS file folder %s cannot be opened\n\n", folder_pathname );
        fprintf( stderr, "ERROR ===> FITS file folder %s cannot be opened\n\n", folder_pathname );
        return(1);
    }


    while( entry = readdir( directoryfolder ) ) {      //... listing of VID files

        if( strstr( entry->d_name, ".vid" ) != NULL )  fprintf( filelisting, "%s\n", entry->d_name );

    }

    fclose( filelisting );

    closedir( directoryfolder );

    return(0);


#endif /************* End of WINDOWS or LINUX *********************************/

}


//################################################################################
//################################################################################
//
//                Sleep functions for Windows and Linux
//
//################################################################################
//################################################################################

void  Delay_msec( int milliseconds )
{

#ifdef _WIN32 /************* WINDOWS ******************************/

    Sleep( milliseconds );

#else /********************** LINUX *******************************/

    sleep( (unsigned int)( milliseconds / 1000 ) );

#endif /***********************************************************/

}

//################################################################################
//################################################################################
//
//             Split a full pathname in to a folder and filename
//
//################################################################################
//################################################################################

void   SplitPathname( char* fullpathname, char* folder_pathname, char* filename )
{
long   k, namelength, lastslash, slash;


    namelength = (long)strlen( fullpathname );

    #ifdef _WIN32
        slash = 92;
    #else  //UNIX wants forward slashes
        slash = 47;
    #endif

    lastslash = (long)( strrchr( fullpathname, (int)slash ) - fullpathname );

    strcpy( folder_pathname, fullpathname ); //... include last slash

    folder_pathname[lastslash+1] = '\0';
    folder_pathname[lastslash+2] = '\0';

    for( k=lastslash+1; k<namelength; k++ ) filename[k-lastslash-1] = fullpathname[k];

    filename[namelength-lastslash-1] = '\0';
    filename[namelength-lastslash-0] = '\0';

}


//################################################################################
//################################################################################
//
//             Build a full pathname from a folder and filename
//
//################################################################################
//################################################################################

void    BuildPathname( char* datafolder_pathname, char* imagery_filename, char* imagery_pathname )
{
    long    kNewLine;


    //... make sure newline is removed

    kNewLine = (long)( strrchr( imagery_filename, '\n' ) - imagery_filename );

    imagery_filename[kNewLine+0] = '\0';
    imagery_filename[kNewLine+1] = '\0';

    strcpy( imagery_pathname, datafolder_pathname );  //... build full pathname to image

    strcat( imagery_pathname, imagery_filename       );

}



//################################################################################
//################################################################################
//
//    Initialize thread local random seed state
//    This function is internal and only used on Linux platforms.
//
//    GetRandom_Init_State sets up thread-local storage for random state and must
//    be called exactly once in a program. It is intended to be wrapped in a call
//    to pthread_once.
//
//################################################################################
//################################################################################

#ifndef _WIN32 /************* LINUX *******************************************/
static void GetRandom_Init_State()
{
    pthread_key_create(&getrandom_seed_key, &free);
}
#endif

//################################################################################
//################################################################################
//
//    Get a uniformly distributed psueodorandom number = [0,RAND_MAX].
//
//    On Windows, this is simply a wrapper for rand() since it supplies its own
//    thread-local random state.
//
//    On Linux, this checks & initializes thread-local random state if needed
//    before calling a reentrant random number generator function to avoid
//    waiting on global locks. This substantially improves threaded performance
//
//################################################################################
//################################################################################

int GetRandom()
{

#ifdef _WIN32 /************* WINDOWS ******************************/

    return rand();

#else /********************** LINUX *******************************/

    unsigned int *seed;

    pthread_once(&getrandom_init_once, &GetRandom_Init_State);
    seed = pthread_getspecific(getrandom_seed_key);

    if (seed == NULL) {
        seed = malloc(sizeof(*seed));
        *seed = 0;
        pthread_setspecific(getrandom_seed_key, seed);
    }

    return rand_r(seed);

#endif /***********************************************************/

}
