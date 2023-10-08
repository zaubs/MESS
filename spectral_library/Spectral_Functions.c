
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Spectral functions  -  C++ file
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 
// Input functions for reading configs, element emission lines, partition functions
// IO and processing functions for extinction and responsivity
// Processing functions for computing an emission line spectrum and abundances 
//
//    Date     Ver   Author        Description
// ----------  ----  ------------  ---------------------------------------------------------
// 2021-05-12  1.00  Gural	       Derived from CAMS-Spectral code base
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double  SPversion = 1.00;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%   Include statements for C library functions   %%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#pragma warning(disable: 4996)  // disable warning on strcpy, fopen, ...

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Spectral_Functions.h"
#include "System_FileFunctions.h"


#ifdef _WIN32 /************* WINDOWS ******************************/

    #include <windows.h>    // string fncts, malloc/free, system

#else  /********************  LINUX  ******************************/

    #include <dirent.h>

#endif /***********************************************************/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%       Standard meteor spectral lines       %%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define  NLINES  61

double   meteorwave_nm[NLINES] = { 382.1509, 382.6966,                // Fe
								   383.0442, 383.3391, 383.9381,      // Mg
								   386.1005,                          // Fe
								   393.4777, 396.9592,                // Ca+
								   422.7918,                          // Ca
								   427.2961, 430.9114, 432.6978,      // Fe
								   438.4776, 440.5987, 441.6362,      // Fe
								   516.8760, 517.4125, 518.5048,      // Mg
								   527.1004, 532.9520,                // Fe
								   537.2983, 539.8628, 540.7277,      // Fe
								   543.1205, 544.8430,                // Fe
								   557.7340,                          // Forbidden O
								   558.8307, 561.7202,                // Fe
								   589.1583, 589.7558,                // Na
								   615.7675, 615.8482, 615.9891,      // O
								   630.0300, 636.3770,                // Forbidden O
								   634.8864, 637.3132,                // Si+
								   742.5686, 744.4348, 747.0369,      // N
								   766.7021, 770.1093,                // K
								   777.4083, 777.6305, 777.7528,      // O
								   818.7111, 819.0263, 820.2611,      // N
								   821.2972, 822.5389, 824.4655,      // N
								   844.8568, 844.8680, 844.9079,      // O
								   857.0089, 859.6361, 863.1605,      // N
								   865.8256, 868.2666, 868.5788,      // N
								   868.8535 };                        // N


//======== Display color palette for element specific fiducials

typedef  int  color[3];

color   DODGERBLUE = { 30, 144, 255 };   // Fe
color        GREEN = { 0, 255,   0 };   // Mg
color       PURPLE = { 160,  32, 240 };   // Ca+
color         CYAN = { 0, 255, 255 };   // Forbidden O
color       YELLOW = { 255, 255,   0 };   // Na
color    ORANGERED = { 255,  69,   0 };   // O
color          TAN = { 210, 180, 140 };   // Si+
color      DARKRED = { 139,   0,   0 };   // N
color       VIOLET = { 238, 130, 238 };   // K

/*
int   meteorwave_color[NLINES][3] = { DODGERBLUE,  DODGERBLUE,               // Fe
									  GREEN,       GREEN,       GREEN,       // Mg
									  DODGERBLUE,                            // Fe
									  PURPLE,      PURPLE,                   // Ca+
									  PURPLE,                                // Ca
									  DODGERBLUE,  DODGERBLUE,  DODGERBLUE,  // Fe
									  DODGERBLUE,  DODGERBLUE,  DODGERBLUE,  // Fe
									  GREEN,       GREEN,       GREEN,       // Mg
									  DODGERBLUE,  DODGERBLUE,               // Fe
									  DODGERBLUE,  DODGERBLUE,  DODGERBLUE,  // Fe
									  DODGERBLUE,  DODGERBLUE,               // Fe
									  CYAN,                                  // Forbidden O
									  DODGERBLUE,  DODGERBLUE,               // Fe
									  YELLOW,      YELLOW,                   // Na
									  ORANGERED,   ORANGERED,   ORANGERED,   // O
									  CYAN,        CYAN,                     // Forbidden O
									  TAN,         TAN,                      // Si+
									  DARKRED,     DARKRED,     DARKRED,     // N
									  VIOLET,      VIOLET,                   // K
									  ORANGERED,   ORANGERED,   ORANGERED,   // O
									  DARKRED,     DARKRED,     DARKRED,     // N
									  DARKRED,     DARKRED,     DARKRED,     // N
									  ORANGERED,   ORANGERED,   ORANGERED,   // O
									  DARKRED,     DARKRED,     DARKRED,     // N
									  DARKRED,     DARKRED,     DARKRED,     // N
									  DARKRED };                             // N
*/




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%                   Functions                    %%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Read the spectral configuration file's parameters
//-----------------------------------------------------------------------------------------

void  ReadSpectralConfigfile( char *pathname, struct specconfiguration *spconfig )
{
int    k;
FILE  *spconfigfile;
char   text[128];

                           
                           
      //========== open configuration file
            
      if( ( spconfigfile = fopen(pathname, "rt" )) == NULL )  {
           printf(" Cannot open spectral config file %s for reading in ReadSpectralConfigfile \n", pathname);
		   Delay_msec(15000);
		   exit(1);
      }
 
	  
	  //========== read the configuration parameters applicable to ALL cameras   	  

      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->version               );

	  if( spconfig->version != 1.00 )  {
           printf(" Spectral config file %s must be version 1.00 \n", pathname);
		   fclose( spconfigfile );
		   Delay_msec(15000);
		   exit(1);
      }
 

	  //========== read the configuration parameters applicable to ALL cameras   	  

      fscanf( spconfigfile, "%[^=]= %ld",  text, &spconfig->order4spcalib         );
      fscanf( spconfigfile, "%[^=]= %ld",  text, &spconfig->rowdelta              );

      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->min_cal_wavelength_nm );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->max_cal_wavelength_nm ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->del_wavelength_nm     ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->min_fit_wavelength_nm );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->max_fit_wavelength_nm ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->minspan_wavelength_nm );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->smooth_bandwidth_nm   );
	  fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->faintest_star_vmag    );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->airmass_limit         ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->fading_coef           ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->coin_time_tolerance   ); 

      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->min_lo_exc_temp       ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->max_lo_exc_temp       ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->step_lo_exc_temp      ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->nominal_lo_exc_temp   ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->nominal_hi_exc_temp   ); 
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->nominal_sigma0        ); 
	  
      fscanf(spconfigfile,  "%[^=]= %lf", text, &spconfig->grating_offnormal_deg );
	  fscanf(spconfigfile,  "%[^=]= %lf", text, &spconfig->default_roll_deg      );
	  fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->default_pitch_deg     );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->default_yaw_deg       );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->default_ne            );
      fscanf( spconfigfile, "%[^=]= %lf", text, &spconfig->default_hot2warm      );

      fscanf( spconfigfile, "%[^=]= %ld", text, &spconfig->ncams_grating         );
      fscanf( spconfigfile, "\n"                                                 );
      fscanf( spconfigfile, "%[^\n]\n",   text                                   );
      fscanf( spconfigfile, "%[^\n]\n",   text                                   );
      fscanf( spconfigfile, "%[^\n]\n",   text                                   );

	  if( spconfig->ncams_grating > MAXGRATINGS )  {
           printf(" Increase MAXGRATINGS dimension in SP_CalibrationFunctions.h \n" );
		   Delay_msec(15000);
		   exit(1);
	  }

	  //========== read the configuration parameters applicable to EACH camera

	  for( k=0; k<spconfig->ncams_grating; k++ )  {
	       fscanf( spconfigfile, "%ld %lf", &spconfig->camnum[k], &spconfig->linespermm[k] );
	  }


      //========== close the configuration file

      fclose( spconfigfile );


	  //========== check wavelength ranges

	  if( spconfig->min_fit_wavelength_nm < spconfig->min_cal_wavelength_nm )  
		  spconfig->min_fit_wavelength_nm = spconfig->min_cal_wavelength_nm;

	  if( spconfig->max_fit_wavelength_nm > spconfig->max_cal_wavelength_nm )  
		  spconfig->max_fit_wavelength_nm = spconfig->max_cal_wavelength_nm;


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Allocate memory for the gratinginfo structure for number of grating cameras in use.
//   Use FreeMemorySpectralCALfile to release the memory.
//-----------------------------------------------------------------------------------------

int   AllocMemorySpectralCalWavelengths( struct speccalibration *spcalib, int ncams_spcalib, int nwavelengths, double startwavelength_nm, double deltawavelength_nm)
{

	//========== Allocate memory for grating information 
	  
	spcalib->ncams_spcalib = ncams_spcalib;

	spcalib->gratinfo = (struct gratinginfo*) malloc( ncams_spcalib * sizeof( struct gratinginfo ) );

	if( spcalib->gratinfo == NULL )  {
		printf(" Memory not allocated for spcalib.gratinfo structures in AllocMemorySpectralCALfile\n");
		return(2);
	}


	//========== Allocate memory for wavelengths, responsivity, and extinction

	spcalib->wavelength_nm  = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->ref_star_spec  = (double*)malloc(nwavelengths * sizeof(double));

	spcalib->wcum_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->prev_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->cumm_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->modl_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->prev_extn_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->cumm_extn_spec = (double*)malloc(nwavelengths * sizeof(double));

	spcalib->ord1_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->ord2_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->modl_extn_spec = (double*)malloc(nwavelengths * sizeof(double));

	spcalib->esti_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->wsum_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->aver_resp_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->esti_extn_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->wsum_extn_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->aver_extn_spec = (double*)malloc(nwavelengths * sizeof(double));
	spcalib->exp_kernel     = (double*)malloc(nwavelengths * sizeof(double));

	if (spcalib->wavelength_nm  == NULL || spcalib->ref_star_spec  == NULL || spcalib->wcum_resp_spec == NULL ||
		spcalib->prev_resp_spec == NULL || spcalib->cumm_resp_spec == NULL || spcalib->modl_resp_spec == NULL ||
		spcalib->prev_extn_spec == NULL || spcalib->cumm_extn_spec == NULL || spcalib->ord1_resp_spec == NULL ||
		spcalib->ord2_resp_spec == NULL || spcalib->modl_extn_spec == NULL || spcalib->esti_resp_spec == NULL ||
		spcalib->wsum_resp_spec == NULL || spcalib->aver_resp_spec == NULL || spcalib->esti_extn_spec == NULL ||
		spcalib->wsum_extn_spec == NULL || spcalib->aver_extn_spec == NULL || spcalib->exp_kernel     == NULL ) {
		printf(" Memory not allocated for spcalib structures in AllocMemorySpectralCALfile\n");
		return(3);
	}


	//========== Set the user defined wavelengths

	spcalib->nwavelengths = nwavelengths;

	for (int kwave = 0; kwave < nwavelengths; kwave++) {
		spcalib->wavelength_nm[kwave] = startwavelength_nm + (double)kwave * deltawavelength_nm;
	}

	  
	return(0);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Read the responsivity and extinction file contents interpolated to the user's
//   defined "nwave" wavelengths vector "wave_spec" given in nanometers.
//   User must have previously called AllocMemorySpectralCALfile.
//-----------------------------------------------------------------------------------------



int   ReadSpectralCALfile( char                     *pathname, 
	                       struct specconfiguration *spconfig,
                           struct speccalibration   *spcalib )

{
FILE   *spcalibfile;
long    kcam_file, kthcam, ncams_file, kwave, kwave_file, nwave_file, kdelta;
long    camnum, year, month, day, hour, minute, second;
char    text[128];
double  roll_deg, pitch_deg, yaw_deg, FileVersion, term;
double  mslope, *wave, *wcumR, *prevR, *cummR, *modlR, *prevE, *cummE;

	


      //==========================================================================================
      //     ASSUMES memory for spcalib->gratinfo is allocated and each gratinfo structure element
      //         has an assignment of camera numbers. The only assigned camera specific info will
      //         be defined for the number of cameras "ncams_spcalib" and their camera numbers
      //         specified in spcalib->gratinfo[*].camnum.
      //==========================================================================================
                           
      //========== open calibration file

	 
            
      if( ( spcalibfile = fopen(pathname, "rt" )) == NULL )  {
           printf(" Cannot open spectral calibration file %s for reading in ReadSpectralCALfile\n", pathname);

	       for( kwave=0; kwave<spcalib->nwavelengths; kwave++ )  {
		        spcalib->wcum_resp_spec[kwave] = 0.0;
				spcalib->prev_resp_spec[kwave] = 0.0;
				spcalib->cumm_resp_spec[kwave] = 0.0;
				spcalib->prev_extn_spec[kwave] = 0.0;
				spcalib->cumm_extn_spec[kwave] = 0.0;
	       }

		   return(1);
	  }

	  

 

	  //========== Read some header info 

      fscanf( spcalibfile, "Version = %lf\n", &FileVersion );
      fscanf( spcalibfile, "NumCams = %ld\n",  &ncams_file );
	  fscanf( spcalibfile, "%[^\n]\n", text );
	  fscanf( spcalibfile, "%[^\n]\n", text );


	  //========== Loop over each camera to read latest date, time, orientation

	  for( kcam_file=0; kcam_file<ncams_file; kcam_file++ )  {
     	  
           fscanf( spcalibfile, "%ld %ld/%ld/%ld %ld:%ld:%ld %lf %lf %lf\n",
			       &camnum, &month, &day, &year, &hour, &minute, &second, &roll_deg, &pitch_deg, &yaw_deg );

		   for( kthcam=0; kthcam< spcalib->ncams_spcalib; kthcam++ )  {

			   if( spcalib->gratinfo[kthcam].camnum == camnum )  {  //... fill in values from previous spec cal

		           spcalib->gratinfo[kthcam].year               = year;
		           spcalib->gratinfo[kthcam].month              = month;
		           spcalib->gratinfo[kthcam].day                = day;
		           spcalib->gratinfo[kthcam].hour               = hour;
		           spcalib->gratinfo[kthcam].minute             = minute;
		           spcalib->gratinfo[kthcam].second             = second;
		           spcalib->gratinfo[kthcam].grating_roll       = (3.141592654 / 180.0) * roll_deg;
		           spcalib->gratinfo[kthcam].grating_pitch      = (3.141592654 / 180.0) * pitch_deg;
		           spcalib->gratinfo[kthcam].grating_yaw        = (3.141592654 / 180.0) * yaw_deg;

		           break;
			   }
		   }
	  }

	  


	  //========== Read the latest spectral responsivity info and extinction coefs

 	  fscanf( spcalibfile, "%[^\n]\n",     text );
	  fscanf( spcalibfile, "%[^=]=%d",     text, &spcalib->extn_fittype  );
	  fscanf( spcalibfile, "%[^=]=%lf",    text, &spcalib->extn_scale    );
	  fscanf( spcalibfile, "%[^=]=%lf",    text, &spcalib->extn_rayleigh );
	  fscanf( spcalibfile, "%[^=]=%lf",    text, &spcalib->extn_oxygen2  );
	  fscanf( spcalibfile, "%[^=]=%lf",    text, &spcalib->extn_ozone    );
	  fscanf( spcalibfile, "%[^=]=%lf",    text, &spcalib->extn_dust     );
	  fscanf( spcalibfile, "%[^=]=%lf",    text, &spcalib->extn_water    );
	  fscanf( spcalibfile, "%[^\n]\n",     text );
	  fscanf( spcalibfile, "%[^=]= %ld\n",  text, &nwave_file );
	  fscanf( spcalibfile, "%[^=]= %lf\n", text, &spcalib->resp_normalization );
	  fscanf( spcalibfile, "%[^\n]\n",     text );
	  fscanf( spcalibfile, "%[^\n]\n",     text );


	  //========== Allocate memory to read file data (this is done because the requested 
	  //              wavelengths "wave_spec" may not match the file wavelengths "wav"
	  //              and the data will need to get interpolated.

      wave  = (double*) malloc( nwave_file * sizeof(double) );
      wcumR = (double*) malloc( nwave_file * sizeof(double) );
      prevR = (double*) malloc( nwave_file * sizeof(double) );
      cummR = (double*) malloc( nwave_file * sizeof(double) );
      modlR = (double*) malloc( nwave_file * sizeof(double) );
      prevE = (double*) malloc( nwave_file * sizeof(double) );
      cummE = (double*) malloc( nwave_file * sizeof(double) );

	  if( wave == NULL   || wcumR == NULL  || prevR == NULL  || cummR == NULL  ||
		  modlR == NULL  || prevE == NULL  || cummE == NULL )  {
		  printf("====> ERROR: Memory not allocated in ReadSpectralCALfile \n");
		  exit(20);
	  } 

	 



	  //========== Read each wavelength and the spectral response

	  for( kwave_file=0; kwave_file<nwave_file; kwave_file++ )  {

		   fscanf( spcalibfile, " %lf %lf %lf %lf %lf %lf %lf", &wave[kwave_file], 
			                                                    &wcumR[kwave_file], 
			                                                    &prevR[kwave_file], 
			                                                    &cummR[kwave_file], 
			                                                    &modlR[kwave_file], 
			                                                    &prevE[kwave_file], 
			                                                    &cummE[kwave_file]  );

		   prevR[kwave_file] *= spcalib->resp_normalization;
		   cummR[kwave_file] *= spcalib->resp_normalization;
		   modlR[kwave_file] *= spcalib->resp_normalization;
	  }

	  fclose( spcalibfile );

	  //========== Interpolate the responsivity as needed to match desired wavelengths

	  kwave_file = 0;

	  for( kwave=0; kwave< spcalib->nwavelengths; kwave++ )  {

		    //-------- Are we below or above the star spectral wavelength band

		    if( spcalib->wavelength_nm[kwave] < wave[0]  ||
				spcalib->wavelength_nm[kwave] > wave[nwave_file-1] )  {

				spcalib->wcum_resp_spec[kwave] = 0.0;
				spcalib->prev_resp_spec[kwave] = 0.0;
				spcalib->cumm_resp_spec[kwave] = 0.0;
				spcalib->modl_resp_spec[kwave] = 0.0;
				spcalib->prev_extn_spec[kwave] = 0.0;
				spcalib->cumm_extn_spec[kwave] = 0.0;

			    continue;
		    }

		    
			//-------- Find the wavelength bin in the star spectrum vector 
		    //           (assumed monotonically increasing in wavelength)

            while( kwave_file < nwave_file-2 )  {

			       if(spcalib->wavelength_nm[kwave] <= wave[kwave_file+1] ) break;

			       kwave_file++;
		    }


		    //-------- Compute spectral slope in the wavelength bin and do linear interpolation

		    mslope = ( wcumR[kwave_file+1] - wcumR[kwave_file] )
			       / (  wave[kwave_file+1]  - wave[kwave_file] );

			spcalib->wcum_resp_spec[kwave] = wcumR[kwave_file]  +  mslope * (spcalib->wavelength_nm[kwave] - wave[kwave_file] );



		    mslope = ( prevR[kwave_file+1] - prevR[kwave_file] )
			       / (  wave[kwave_file+1]  - wave[kwave_file] );

			spcalib->prev_resp_spec[kwave] = prevR[kwave_file]  +  mslope * (spcalib->wavelength_nm[kwave] - wave[kwave_file] );



		    mslope = ( cummR[kwave_file+1] - cummR[kwave_file] )
			       / (  wave[kwave_file+1]  - wave[kwave_file] );

			spcalib->cumm_resp_spec[kwave] = cummR[kwave_file]  +  mslope * (spcalib->wavelength_nm[kwave] - wave[kwave_file] );



		    mslope = ( modlR[kwave_file+1] - modlR[kwave_file] )
			       / (  wave[kwave_file+1]  - wave[kwave_file] );

			spcalib->modl_resp_spec[kwave] = modlR[kwave_file]  +  mslope * (spcalib->wavelength_nm[kwave] - wave[kwave_file] );



		    mslope = ( prevE[kwave_file+1] - prevE[kwave_file] )
			       / (  wave[kwave_file+1]  - wave[kwave_file] );

			spcalib->prev_extn_spec[kwave] = prevE[kwave_file]  +  mslope * (spcalib->wavelength_nm[kwave] - wave[kwave_file] );



		    mslope = ( cummE[kwave_file+1] - cummE[kwave_file] )
			       / (  wave[kwave_file+1]  - wave[kwave_file] );

			spcalib->cumm_extn_spec[kwave] = cummE[kwave_file]  +  mslope * (spcalib->wavelength_nm[kwave] - wave[kwave_file] );



	  }  //... end of output wavelength loop

	  

	  //========== Set the 1st order responsivity spectra from the modeled smoothed response
	  //           Set the 2nd order responsivity to zero (only working with 1st order)
	  //              NOTE: one could override these with lab measurements

	  for (int kwave = 0; kwave < spcalib->nwavelengths; kwave++) {
		  spcalib->ord1_resp_spec[kwave] = spcalib->modl_resp_spec[kwave];
		  spcalib->ord2_resp_spec[kwave] = 0.0;
	  }

	  

	  //========== Set an initial model for extinction for airmass of unity
	  //               90 deg elev, earth's surface height of 6378.16 km

	  ExtinctionModel(spcalib->nwavelengths, spcalib->wavelength_nm, spcalib->modl_extn_spec, 90.0, 6378.16, spcalib);

	  spcalib->fading_coef = spconfig->fading_coef;


	  //========== Set the smoothing exponential weighting kernel off emmision line center

	  for (kdelta = 0; kdelta < spcalib->nwavelengths; kdelta++) {
		  term = (double)kdelta * spconfig->del_wavelength_nm / spconfig->smooth_bandwidth_nm;
		  spcalib->exp_kernel[kdelta] = exp(-0.5 * term * term);
	  }

	  

	  //========== Free temporary memory

	  free( wave );
	  free( wcumR );
	  free( prevR );
	  free( cummR );
	  free( modlR );
	  free( prevE );
	  free( cummE );

	  return(0);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Writes the responsivity and extinction file contents at the user defined wavelengths.
//   The SPCAL*.txt file is dated to the nearest hour using the user input JDT (recommend 
//   using local midnight of the night of the processed collection).
//-----------------------------------------------------------------------------------------

int   WriteSpectralCALfile( char                   *SPCALfolder, 
							double                  jdt,
                            struct speccalibration *spcalib  )
{
FILE  *spcalibfile;
char   filename[256];
long   kthcam, kwave;
long   year, month, day, hour, minute, second, milliseconds;
char   fitstring[10];
double max_resp_value;


      //---------- Find the maximum responsivity value for normalization

      max_resp_value = 0.0;
	  for( kwave=0; kwave< spcalib->nwavelengths; kwave++ )  {
		  if( max_resp_value < spcalib->modl_resp_spec[kwave] )  {
			  max_resp_value = spcalib->modl_resp_spec[kwave];
		  }
	  }


      //------------- Build the spectral calibration filename for the given JDT to the nearest hour

      CalendarDateAndTime( jdt, &year, &month, &day, &hour, &minute, &second, &milliseconds );
	  minute       = 0;
	  second       = 0;
	  milliseconds = 0;


	  sprintf( filename, "%sSPCAL_%04ld%02ld%02ld_%02ld0000.txt", SPCALfolder, year, month, day, hour );


      //------------- Ppen calibration file SPCAL_YYYYMMDD_HH0000.txt
            
      if( ( spcalibfile = fopen( filename, "wt" )) == NULL )  {
           printf(" Cannot open spectral calibration file %s for writing \n", filename );
           return(1);
      }


	  //------------- Write some header info, then loop over each camera

      fprintf( spcalibfile, "Version = %5.2lf\n", SPversion );
      fprintf( spcalibfile, "NumCams = %d\n", spcalib->ncams_spcalib );
	  fprintf( spcalibfile, "Cam     Date       Time        Roll     Pitch       Yaw\n");
	  fprintf( spcalibfile, "---  ----------  --------  --------  --------  --------\n");

	  for( kthcam=0; kthcam< spcalib->ncams_spcalib; kthcam++ )  {
     	  
           fprintf( spcalibfile, "%03ld  %02ld/%02ld/%04ld  %02ld:%02ld:%02ld  %8.2lf  %8.2lf  %8.2lf\n",  
			        spcalib->gratinfo[kthcam].camnum,
	                spcalib->gratinfo[kthcam].month, spcalib->gratinfo[kthcam].day,    spcalib->gratinfo[kthcam].year,
	                spcalib->gratinfo[kthcam].hour,  spcalib->gratinfo[kthcam].minute, spcalib->gratinfo[kthcam].second,
                    spcalib->gratinfo[kthcam].grating_roll  * 180.0 / 3.141592654,
                    spcalib->gratinfo[kthcam].grating_pitch * 180.0 / 3.141592654, 
                    spcalib->gratinfo[kthcam].grating_yaw   * 180.0 / 3.141592654 );
	  }


	  //----------- Write out extinction/responsivity headers

	  if( spcalib->extn_fittype == 1 )  strcpy( fitstring, "Same-Day " );
	  else                            strcpy( fitstring, "Multi-Day" );

	  fprintf( spcalibfile, "-------------------------------------------------------\n");
	  fprintf( spcalibfile, "Extinction Fit: %s = %d\n", fitstring, spcalib->extn_fittype  );
	  fprintf( spcalibfile, "Extinction Coef Scale     = %lf\n",    spcalib->extn_scale    );
	  fprintf( spcalibfile, "Extinction Coef Rayleigh  = %lf\n",    spcalib->extn_rayleigh );
	  fprintf( spcalibfile, "Extinction Coef Oxygen2   = %lf\n",    spcalib->extn_oxygen2  );
	  fprintf( spcalibfile, "Extinction Coef Ozone     = %lf\n",    spcalib->extn_ozone    );
	  fprintf( spcalibfile, "Extinction Coef Dust      = %lf\n",    spcalib->extn_dust     );
	  fprintf( spcalibfile, "Extinction Coef Water     = %lf\n",    spcalib->extn_water    );
	  fprintf( spcalibfile, "-------------------------------------------------------\n");
	  fprintf( spcalibfile, "Number of wavelengths      = %d\n",    spcalib->nwavelengths );
	  fprintf( spcalibfile, "Responsivity normalization = %lf\n",   max_resp_value );
	  fprintf( spcalibfile, "-------------------------------------------------------\n");


	  //----------- Write out extinction/responsivity values at each wavelength

	  fprintf( spcalibfile,  "Wave(nm)    WeightSum   DayResp     CumResp     SmoothResp  DayExtn     CumExtn\n" );

	  for( kwave=0; kwave< spcalib->nwavelengths; kwave++ )  {

		  fprintf( spcalibfile, "%10.4lf  %10.3lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf\n", 
			  spcalib->wavelength_nm[kwave],
			  spcalib->wcum_resp_spec[kwave],
			  spcalib->prev_resp_spec[kwave] / max_resp_value,
			  spcalib->cumm_resp_spec[kwave] / max_resp_value,
			  spcalib->modl_resp_spec[kwave] / max_resp_value,
			  spcalib->prev_extn_spec[kwave],
			  spcalib->cumm_extn_spec[kwave] );
	  }

	  fclose( spcalibfile );


	  return(0);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void   FreeMemorySpectralCALfile( struct speccalibration *spcalib ) 
{
	free( spcalib->gratinfo );

	free(spcalib->wavelength_nm);
	free(spcalib->ref_star_spec);

	free(spcalib->wcum_resp_spec);
	free(spcalib->prev_resp_spec);
	free(spcalib->cumm_resp_spec);
	free(spcalib->modl_resp_spec);
	free(spcalib->prev_extn_spec);
	free(spcalib->cumm_extn_spec);

	free(spcalib->ord1_resp_spec);
	free(spcalib->ord2_resp_spec);
	free(spcalib->modl_extn_spec);

	free(spcalib->esti_resp_spec);
	free(spcalib->wsum_resp_spec);
	free(spcalib->aver_resp_spec);
	free(spcalib->esti_extn_spec);
	free(spcalib->wsum_extn_spec);
	free(spcalib->aver_extn_spec);
	free(spcalib->exp_kernel    );


	spcalib->gratinfo = NULL;

	spcalib->wavelength_nm  = NULL;
	spcalib->ref_star_spec  = NULL;

	spcalib->wcum_resp_spec = NULL;
	spcalib->prev_resp_spec = NULL;
	spcalib->cumm_resp_spec = NULL;
	spcalib->modl_resp_spec = NULL;
	spcalib->prev_extn_spec = NULL;
	spcalib->cumm_extn_spec = NULL;

	spcalib->ord1_resp_spec = NULL;
	spcalib->ord2_resp_spec = NULL;
	spcalib->modl_extn_spec = NULL;

	spcalib->esti_resp_spec = NULL;
	spcalib->wsum_resp_spec = NULL;
	spcalib->aver_resp_spec = NULL;
	spcalib->esti_extn_spec = NULL;
	spcalib->wsum_extn_spec = NULL;
	spcalib->aver_extn_spec = NULL;
	spcalib->exp_kernel     = NULL;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Returns the filename of the most recent SPCAL file that is not later in time than
//   the user specified JDT (e.g. dawn of the to be processed night's imagery).
//-----------------------------------------------------------------------------------------

int     GetMostRecentSPCAL( char *folder_pathname, char *file_firstpart, double jdtlast, char *latest_filename )
{
double   jdtcal, jdtdiff;

char     fullfilename[512];
char     text[128];
long     year, month, day, hour;


	   //========== Loop over the list of spectral calibration filenames for this camera
	   //             to find the latest in time up to the date specified by the user.
	   //             Branch for Windows versus Linux implementations.


       #ifdef _WIN32 /*#######################  WINDOWS  #########################*/


	   HANDLE              hfindFile;
       WIN32_FIND_DATA     findFileStruct;
	   char                genericpathname[512];


	   strcpy(genericpathname, folder_pathname);

	   strcat(genericpathname, file_firstpart);

	   strcat(genericpathname, "*.txt");

	   hfindFile = FindFirstFile( genericpathname, &findFileStruct );

       if( hfindFile == INVALID_HANDLE_VALUE )  {
	       printf(" FindFirstFile: Invalid handle returned searching for %s files\n", file_firstpart );
		   printf(" No %s* file found with valid date/time\n", file_firstpart );
           return(1);
       }

	   strcpy( fullfilename, findFileStruct.cFileName );

       sscanf( fullfilename, "%[^_]_%4ld%2ld%2ld_%02ld0000.txt", text, &year, &month, &day, &hour );

	   jdtcal = JulianDateAndTime( year, month, day, hour, 0, 0, 0 ); 

	   strcpy( latest_filename, findFileStruct.cFileName );
      
	   jdtdiff = jdtlast - jdtcal;


	   //========== Pick the closest response file in time no later than the 
	   //              current date (in case the system changed at a future date) 

	   while( FindNextFile( hfindFile, &findFileStruct ) )  {

		   strcpy( fullfilename, findFileStruct.cFileName );

		   //-------- parse the calibration filename for date and time
		   //              SPCAL###_YYYYMMDD_HH0000.txt

           sscanf( fullfilename, "%[^_]_%4ld%2ld%2ld_%02ld0000.txt", text, &year, &month, &day, &hour );

	       jdtcal = JulianDateAndTime( year, month, day, hour, 0, 0, 0 ); 



		   if( ( fabs( jdtlast - jdtcal ) <= fabs( jdtdiff ) ) &&
			   (                 jdtcal   <= jdtlast         )     )  {

			   strcpy( latest_filename, findFileStruct.cFileName );

		       jdtdiff = jdtlast - jdtcal;
		   }

	   } //... end of loop searching for 

       FindClose( hfindFile );



       #else  /*#######################  LINUX  #######################*/


	   int     icount;
	   struct  dirent *entry;
	   DIR    *directory_ptr;


	   icount = 0;

	   directory_ptr = opendir(folder_pathname);

	   if (directory_ptr != NULL) {

		   while (entry = readdir(directory_ptr))  {

			   strcpy(fullfilename, entry->d_name);

			   //--------- Only process *.txt files

			   if ((strstr(fullfilename, file_firstpart) != NULL) && (strstr(fullfilename, ".txt") != NULL)) {

				   //........ Parse the calibration filename for date and time  SPCAL###_YYYYMMDD_HH0000.txt

				   sscanf(fullfilename, "%[^_]_%4ld%2ld%2ld_%02ld0000.txt", text, &year, &month, &day, &hour);

				   jdtcal = JulianDateAndTime(year, month, day, hour, 0, 0, 0);

				   if ((icount == 0) || ((fabs(jdtlast - jdtcal) <= fabs(jdtdiff)) && (jdtcal <= jdtlast))) {
					   strcpy(latest_filename, fullfilename);
					   jdtdiff = jdtlast - jdtcal;
				   }

				   icount += 1;

			   } //... end of test if file has a .txt extension

		   } //... end of while loop looking for any files in the folder

       }  //... end if directory found

	   else  printf(" No SPCAL files found in folder %s  \n", folder_pathname);

	   closedir(directory_ptr);


       #endif /*#######################################################*/



       //========== Ensure we have a null character at the end of the filename string.
       //           Spectral Cal file convention is 25 characters  SPCAL_YYYYMMDD_HH0000.txt 

       latest_filename[25] = '\0';  //SPCAL

       //printf(" Latest Spectral Calibration filename is %s \n", latest_filename );
	   
	   return(0);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Calculate the responsivity and extinction estimates given a spectral measurement.
//   Validity of the responsivity estimate should be assessed after return from this 
//       function by testing the number of airmasses returned and whether it violates  
//       the simple extinction approximation used herein - which is typically good
//       for < 1.5 airmasses or 42 degrees altitude angle above horizon.
//   Validity of the extinction estimate is assumed OK for all airmasses since it uses
//       a multi-day responsivity measure which should be relatively stable day-to-day.
//-----------------------------------------------------------------------------------------

double  ResponsivityExtinctionEstimate( double                  *meas_spectrum,
	                                    double                  *star_spectrum,
	                                    double                   altitude_deg,
	                                    double                   height_observer_km,
	                                    struct speccalibration  *spcalib )
{
	int     kwave;
	double  extinction, Xairmass;


	//======== Compute the airmass

	Xairmass = Airmass(altitude_deg, 6378.16 + height_observer_km);


	//======== Responsivity assumes extinction can be approximated (small enough airmass)
	//         Note responsivity is computed no matter what the airmass is, and the user 
	//            decides to use it or not based on airmass and data quality.

	for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

		extinction = ExtinctionApprox(altitude_deg, height_observer_km, spcalib->wavelength_nm[kwave]);

		if (star_spectrum[kwave] > 0.0 && extinction > 0.0)
			spcalib->esti_resp_spec[kwave] = meas_spectrum[kwave] / star_spectrum[kwave] / extinction;
		else
			spcalib->esti_resp_spec[kwave] = 0.0;

	}


	//======== Extinction uses the cummulative smoothed multi-day first order responsivity (stored in ord1)

	for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

		if (star_spectrum[kwave] > 0.0 && spcalib->ord1_resp_spec[kwave] > 0.0)
			extinction = meas_spectrum[kwave] / star_spectrum[kwave] / spcalib->ord1_resp_spec[kwave];
		else
			extinction = 1.0e-6;
		
		extinction = pow(extinction, 1.0 / Xairmass);

		if (extinction > 1.0) extinction = 1.0;

		spcalib->esti_extn_spec[kwave] = extinction;
	}


	return(Xairmass);


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Initialize, update the daily, or the multi-day responsivity.
//   ResponsivityExtinctionEstimate should be called prior to this function call so
//       that the estimated responsivity is infilled into spcalib.esti_resp_spec.
//   User provides the weight which could be unity, or a function of Vmagnitude, or
//       dependent on some other spectrum quality measure. It is used only for REDAILY.
//
//   REflag = 0 = RECLEAR   Clear the daily running weight and weighted signal vectors
//   REflag = 1 = REDAILY   Add to the daily running sums
//   REflag = 2 = REMULTI   Final daily average and update the multi-day running sums
//-----------------------------------------------------------------------------------------

void    ResponsivityAccumulate(int REflag, double weight, struct speccalibration  *spcalib)
{
	int     kwave, jwave;
	double  sum, numer, denom, kernel;


	//======== Clear the running weight and weighted average vectors

	if (REflag == RECLEAR) {

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {
			spcalib->wsum_resp_spec[kwave] = 0.0;
			spcalib->aver_resp_spec[kwave] = 0.0;
		}

	}


	//======== Add latest responsivity estimate to the running sum vectors
	//         Note this form of expression computes the average on the fly

	else if (REflag == REDAILY) {

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

			sum = spcalib->wsum_resp_spec[kwave] * spcalib->aver_resp_spec[kwave]
				+ spcalib->esti_resp_spec[kwave] * weight;

			spcalib->wsum_resp_spec[kwave] += weight;

			spcalib->aver_resp_spec[kwave] = sum / spcalib->wsum_resp_spec[kwave];
		}

	}


	//======== Merge daily average responsivity into the multi-day responsivity
	//           and save the daily average in prev_resp_spec and a smoothed
	//           responsivity in modl_resp_spec.

	else if (REflag == REMULTI) {

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

			sum = spcalib->wcum_resp_spec[kwave] * spcalib->cumm_resp_spec[kwave]
				+ spcalib->wsum_resp_spec[kwave] * spcalib->aver_resp_spec[kwave];

			spcalib->wcum_resp_spec[kwave] += spcalib->wsum_resp_spec[kwave];

			spcalib->cumm_resp_spec[kwave] = sum / spcalib->wcum_resp_spec[kwave];

			spcalib->prev_resp_spec[kwave] = spcalib->aver_resp_spec[kwave];
		}


		//-------- Compute a smoothed "model" responsivity

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

			numer = 0.0;
			denom = 0.0;

			for (jwave = 0; jwave < spcalib->nwavelengths; jwave++) {

				kernel = spcalib->exp_kernel[labs(kwave - jwave)];

				numer += kernel * spcalib->cumm_resp_spec[jwave];

			    denom += kernel;

			}

			spcalib->modl_resp_spec[kwave] = numer / denom;

			if (spcalib->modl_resp_spec[kwave] < 0.0)  spcalib->modl_resp_spec[kwave] = 0.0;

		}

	}




	//======== REflag not implemented

	else {
		printf("ERROR===> REflag = %d not implemented in ResponsivityAccumulate\n", REflag);
		Delay_msec(15000);
		exit(1);
	}

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Initialize, update the daily, or the multi-day extinction.
//   ResponsivityExtinctionEstimate should be called prior to this function call so
//       that the estimated extinction is infilled into spcalib.esti_extn_spec.
//   User provides the weight which could be unity, or a function of Vmagnitude, or
//       dependent on some other spectrum quality measure. It is used only for REDAILY.
//
//   REflag = 0 = RECLEAR   Clear the daily running weight and weighted signal vectors
//   REflag = 1 = REDAILY   Add to the daily running sums
//   REflag = 2 = REMULTI   Final daily average and update the multi-day running sums
//-----------------------------------------------------------------------------------------

void      ExtinctionAccumulate(int REflag, double weight, struct speccalibration  *spcalib)
{
	int     kwave;
	double  sum;


	//======== Clear the running weight and weighted average vectors

	if (REflag == RECLEAR) {

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {
			spcalib->wsum_extn_spec[kwave] = 0.0;
			spcalib->aver_extn_spec[kwave] = 0.0;
		}

	}


	//======== Add latest extinction estimate to the running sum vectors
	//         Note this form of expression computes the average on the fly

	else if (REflag == REDAILY) {

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

			sum = spcalib->wsum_extn_spec[kwave] * spcalib->aver_extn_spec[kwave]
				+ spcalib->esti_extn_spec[kwave] * weight;

			spcalib->wsum_extn_spec[kwave] += weight;

			spcalib->aver_extn_spec[kwave] = sum / spcalib->wsum_extn_spec[kwave];
		}

	}


	//======== Merge daily average extinction into a fading memory multi-day extinction
	//           cumm_extn_spec and save the daily average in prev_extn_spec.

	else if (REflag == REMULTI) {

		for (kwave = 0; kwave < spcalib->nwavelengths; kwave++) {

			spcalib->cumm_extn_spec[kwave] = (1.0 - 1.0 / spcalib->fading_coef) * spcalib->cumm_extn_spec[kwave]
				                                 + (1.0 / spcalib->fading_coef) * spcalib->aver_extn_spec[kwave];

			spcalib->prev_extn_spec[kwave] = spcalib->aver_extn_spec[kwave];
		}

		FitExtinctionCoefs(spcalib->nwavelengths, spcalib->wavelength_nm, spcalib->cumm_extn_spec, spcalib);

	}


	//======== REflag not implemented

	else {
		printf("ERROR===> REflag = %d not implemented in ExtinctionAccumulate\n", REflag);
		Delay_msec(15000);
		exit(1);
	}



}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Weighting schema for responsivity and extinction averaging. Currently has unity
//      weights or V-magnitude based weights.
//-----------------------------------------------------------------------------------------

double     Weighting(int WEIGHTSflag, double Vmagnitude)
{

	if (WEIGHTSflag == UNITYWEIGHTS)  return(1.0);

	else if (WEIGHTSflag == VMAGWEIGHTS) return(pow(10.0, -Vmagnitude) );

	else {
		printf("ERROR===> WEIGHTSflag = %d not implemented in Weighting\n", WEIGHTSflag);
		Delay_msec(15000);
		exit(1);
	}


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Read all the available elements data from the ElementsInfo.txt file located in the
//   "elements_data_folder" into the elements_data structure. Also allocates memory for
//   later spectra generation computations per element to be done at the user defined 
//   input wavelengths.
//-----------------------------------------------------------------------------------------

void  ReadElementsData( int     nwavelengths,   // number of user defined wavelengths
	                    double *wavelength_nm,  // user defined wavelengths for spectra calcs
	                    char   *elements_data_folder,
	                    struct  elements_data  *elemdata)
{
	int     kwave, kelem, jelem, nelem, kneutral, kline, nstring, kpf, Npf, Tpf[MAXTPF];
	double  wvac, wair, ephoton, log10gf, avalue, code, specvalue;
	double  Elower, Jlower, Eupper, Jupper, gf;
	char    text[256];
	char    elements_data_pathname[256];
	struct  element_lines_spectrum  *els;
	FILE   *elemfile;


	//========== open the element compendium file

	strcpy(elements_data_pathname, elements_data_folder);

	strcat(elements_data_pathname, "ElementsInfo.txt");

	if ((elemfile = fopen(elements_data_pathname, "rt")) == NULL) {
		printf(" Cannot open file %s\n", elements_data_pathname);
		Delay_msec(15000);
		exit(1);
	}

	//========== Determine number of elements total (neutral + ionized + molecules)

	fscanf(elemfile, "%[^\n]\n", text);  // read the 5 header lines
	fscanf(elemfile, "%[^\n]\n", text);
	fscanf(elemfile, "%[^\n]\n", text);
	fscanf(elemfile, "%[^\n]\n", text);
	fscanf(elemfile, "%[^\n]\n", text);

	nelem = 0;
	while (feof(elemfile) == 0) {
		fscanf(elemfile, "%[^\n]\n", text);
		// printf("%s\n",text);
		nelem++;
	}

	printf("Number of elements read:\n");
	printf("%d\n", nelem);
	printf("****\n");

	elemdata->nelements = nelem;

	//========== Allocate memory for the main structure file

	elemdata->els = (struct element_lines_spectrum *) malloc(nelem * sizeof(struct element_lines_spectrum));

	if (elemdata->els == NULL) {
		printf(" Memory not allocated for els structure in ReadElementsData\n");
		Delay_msec(15000);
		exit(2);
	}


	elemdata->fit_element_index = (int*)malloc(nelem * sizeof(int));

	if (elemdata->fit_element_index == NULL) {
		printf(" Memory not allocated for fit_element_index in ReadElementsData\n");
		Delay_msec(15000);
		exit(2);
	}

	elemdata->neutral_index = (int*)malloc(nelem * sizeof(int));

	if (elemdata->neutral_index == NULL) {
		printf(" Memory not allocated for neutral_index in ReadElementsData\n");
		Delay_msec(15000);
		exit(2);
	}


	//========== Loop over each element (neutrals, ions, molecules), 
	//              which are on separate text lines in the file

	rewind(elemfile);

	fscanf(elemfile, "%[^\n]\n", text);          //... read the first dashed header line
	fscanf(elemfile, "%[^=]= %d\n", text, &Npf); //... read the second header line
	fscanf(elemfile, "%[^\n]\n", text);          //... read the third header line

	if (Npf > MAXTPF) {
		printf(" Increase MAXTPF since %d > %d in Spectral_Functions.h\n", Npf, MAXTPF);
		Delay_msec(15000);
		exit(3);
	}

	//for (kpf = 0; kpf < 4; kpf++)  fscanf(elemfile, "%s", text);
	fscanf(elemfile, "%[^=]=", text);
	for (kpf = 0; kpf < Npf + 1; kpf++)  fscanf(elemfile, "%dK", &Tpf[kpf]);
	fscanf(elemfile, "%[^\n]\n", text);  //... read to the end of the fourth header line

	fscanf(elemfile, "%[^\n]\n", text);  //... read last dashed header line



	//========== Allocate memory to assign wavelengths

	elemdata->nwave = nwavelengths;

	elemdata->wave    = (double*)malloc(nwavelengths * sizeof(double));
	elemdata->lockfit = (double*)malloc(nwavelengths * sizeof(double));
	elemdata->corrfit = (double*)malloc(nwavelengths * sizeof(double));
	elemdata->dispfit = (double*)malloc(nwavelengths * sizeof(double));

	if (elemdata->wave == NULL || elemdata->lockfit == NULL || elemdata->corrfit == NULL || elemdata->dispfit == NULL) {
		printf(" Memory not allocated in ReadElementsData\n");
		Delay_msec(15000);
		exit(4);
	}

	for (kwave = 0; kwave < elemdata->nwave; kwave++)  elemdata->wave[kwave] = wavelength_nm[kwave];
	for (kwave = 0; kwave < elemdata->nwave; kwave++)  elemdata->lockfit[kwave] = 0.0;
	for (kwave = 0; kwave < elemdata->nwave; kwave++)  elemdata->corrfit[kwave] = 0.0;
	for (kwave = 0; kwave < elemdata->nwave; kwave++)  elemdata->dispfit[kwave] = 0.0;

	//========== Allocate memory for output spectra of each element

	kneutral = 0;



	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		els = &elemdata->els[kelem];

		fscanf(elemfile, "%d %s %lf %d %lf %lf %lf %lf",
			&els->init_fitflag,
			els->element_filename,
			&els->elementcode,
			&els->warmhot,
			&els->init_abundance,
			&els->ionenergy,
			&els->Vo,
			&els->c);

		if (els->elementcode == 12.0)  els->init_fitflag = FITTING;  //... always fit with Mg to start

		els->init_fitflag = FITLESS;  //... override all elements to not fitting at the start

		for (kpf = 0; kpf < Npf; kpf++)  fscanf(elemfile, "%lf", &els->partfunc[kpf]);

		fscanf(elemfile, "%[^\n]\n", text);  //... read to the end of the line

		for (kpf = 0; kpf < Npf; kpf++)  els->Tpf[kpf] = (double)Tpf[kpf];

		els->Npf = Npf;


		els->user_fitflag = els->init_fitflag;

		nstring = strstr(els->element_filename, ".txt") - els->element_filename;

		strncpy(text, els->element_filename, nstring);

		text[nstring + 0] = '\0';
		text[nstring + 1] = '\0';

		nstring = strstr(els->element_filename, "-") - els->element_filename;

		if (nstring == 1)  strcpy(els->element_string, " ");  // add a leading blank if element acronym is one letter long
		else                strcpy(els->element_string, "");

		strcat(els->element_string, text);

		els->ioncode = 1 + ((int)(els->elementcode * 100.0) % 100);  // 1, 2, ... = Neutral, first ionization ...
																     // 51 = molecule
		// printf("*** Ion code = %d\n", els->ioncode);
		//-------- Add this element to list of only neutrals if it is a neutral

		if (els->ioncode == 1 || els->ioncode == 51) {
			elemdata->neutral_index[kneutral] = kelem;
			kneutral++;
		}


		//-------- Allocate memory for output spectra, assign wavelengths

		els->speclo = (double*)malloc(nwavelengths * sizeof(double));
		els->spechi = (double*)malloc(nwavelengths * sizeof(double));

		if (els->speclo == NULL && els->spechi == NULL) {
			printf(" Memory not allocated for spec wrt element %s in ReadElementsData\n", els->element_string);
			Delay_msec(15000);
			exit(4);
		}

		for (kwave = 0; kwave < elemdata->nwave; kwave++)  els->speclo[kwave] = 0.0;
		for (kwave = 0; kwave < elemdata->nwave; kwave++)  els->spechi[kwave] = 0.0;

	}  //...end of element loop to read basic parameters

	elemdata->nneutrals = kneutral;

	fclose(elemfile);

	setvbuf (stdout, NULL, _IONBF, 0); 

	//========== Check the file was ordered correctly 

	for (kelem = 1; kelem < elemdata->nelements; kelem++) {

		/*
		if( elemdata->els[kelem].elementcode <= elemdata->els[kelem-1].elementcode )  {
			printf(" Error with file %s as elements not in increasing elementcode order\n", elements_data_pathname );
		    Delay_msec(15000);
			exit(9);
		}
		*/

		if (elemdata->els[kelem].ioncode == 2 && elemdata->els[kelem].elementcode - elemdata->els[kelem - 1].elementcode > 0.015) {
			printf(" Error with file %s having singly ionized element %lf not following a neutral\n", elements_data_pathname, elemdata->els[kelem].elementcode);
			Delay_msec(15000);
			exit(9);
		}

		if (elemdata->els[kelem].ioncode > 2 && elemdata->els[kelem].ioncode < 51) {
			printf(" Error with file %s having doubly ionized entry %lf\n", elements_data_pathname, elemdata->els[kelem].elementcode);
			Delay_msec(15000);
			exit(9);
		}

	}


	//========== Loop over each element's neutrals, ions, and molecules to read the spectral multiple line data

	printf("%d nelements\n", elemdata->nelements);

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		els = &elemdata->els[kelem];

		//printf("%d kelem\n", kelem);

		//-------- Open the element spectral line file

		strcpy(elements_data_pathname, elements_data_folder);

		strcat(elements_data_pathname, els->element_filename);

		if ((elemfile = fopen(elements_data_pathname, "rt")) == NULL) {
			printf(" Cannot open file %s\n", elements_data_pathname);
			Delay_msec(15000);
			exit(5);
		}

		//-------- Read the number of lines and allocate memory

		fscanf(elemfile, "Number of lines in output :%d\n", &els->nlines);

		if (els->nlines > 0) {

			els->wave_vacuum = (double*)malloc(els->nlines * sizeof(double));
			els->gf013divw3  = (double*)malloc(els->nlines * sizeof(double));
			els->Eupdiv8     = (double*)malloc(els->nlines * sizeof(double));

			if (els->wave_vacuum == NULL || els->gf013divw3 == NULL || els->Eupdiv8 == NULL) {
				printf(" Memory not allocated for line parameters wrt element %s in ReadElementsData\n", els->element_string);
				Delay_msec(15000);
				exit(7);
			}
			

			//-------- Read the individual line parameters

			fscanf(elemfile, "%[^\n]\n", text); // read the four header lines
			fscanf(elemfile, "%[^\n]\n", text);
			fscanf(elemfile, "%[^\n]\n", text);
			fscanf(elemfile, "%[^\n]\n", text);


			if (els->ioncode < 51) {  //--- read element and ion data line content

				

				for (kline = 0; kline < els->nlines; kline++) {

					fscanf(elemfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %[^\n]\n",
						&wvac, &wair, &ephoton, &log10gf, &avalue, &code, &Elower, &Jlower, &Eupper, &Jupper, text);

					gf = pow(10.0, log10gf);

					els->wave_vacuum[kline] = wvac;
					els->gf013divw3[kline] = gf * 0.013248 / pow(wvac, 3.0);
					els->Eupdiv8[kline] = -Eupper / 8.617343e-5;

				};


			}


			else {  //--- read molecule data line content (assumes each spectral value is a line of emmision)

				for (kline = 0; kline < els->nlines; kline++) {

					fscanf(elemfile, "%lf %lf\n",
						&wvac, &specvalue);

					els->wave_vacuum[kline] = wvac;
					els->gf013divw3[kline] = specvalue;
					els->Eupdiv8[kline] = 0.0;

				};


			}


			//--------- Check wavelengths are monotonically increasing in the array

			for (kline = 1; kline < els->nlines; kline++) {

				if (els->wave_vacuum[kline] < els->wave_vacuum[kline - 1]) {

					printf(" Element file %s does not have monotonically increasing wavelength \n\n", els->element_filename);
					printf(" Entry %d of %lf nm is less than previous line of %lf nm\n\n", kline, els->wave_vacuum[kline], els->wave_vacuum[kline - 1]);
					Delay_msec(15000);
					exit(8);

				}

			};


		} //... end of if any lines in file


		else { //------ No lines in the file so do not fit to this element

			els->init_fitflag = FITLESS;
			els->user_fitflag = FITLESS;

		}


		//--------- Close file (finished reading the element spectral line file)

		fclose(elemfile);


	}  //... end of element line reading loop



	//========== Set the ion indices that correspond to a neutral and vice versa

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {
		// printf("Setting ion indices...\n");
		els = &elemdata->els[kelem];

		els->ionindex = -1; // indicator that there may be no ion in the table for this element
		els->neuindex = -1; // indicator that there may be no neutral in the table for this element

		//-----------------------------------------------------------------------
		//     Checking for Neutral(1), 1st-Ionization(2), .... Molecule(51)
		//-----------------------------------------------------------------------

		//........ If this is a neutral, find its corresponding ion

		if (els->ioncode == 1) {
			// printf("This is neutral. Finding corresponding ion...\n");
			for (jelem = 0; jelem < elemdata->nelements; jelem++) {

				if (fabs(elemdata->els[jelem].elementcode - (els->elementcode + 0.01)) < 1.0e-3) {

					els->ionindex = jelem;
					// printf("%d\n", els->ionindex);

					break;
				}
			}

		}


		//........ If this is an ion, find its corresponding neutral

		if (els->ioncode == 2) {
			// printf("This is an ion. Finding corresponding neutral...\n");
			for (jelem = 0; jelem < elemdata->nelements; jelem++) {

				if (fabs(elemdata->els[jelem].elementcode - (els->elementcode - 0.01)) < 1.0e-3) {

					els->neuindex = jelem;
					// printf("%d\n", els->neuindex);
					break;
				}
			}

		}


		//........ If this is a molecule, there is no corresponding ion at this time

		if (els->ioncode == 51) {

			els->ionindex = -1;

		}


	} //... end of loop over all elements and ions



}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Frees all the memory concerning elements content for the elements_data structure
//-----------------------------------------------------------------------------------------

void  FreeElementsMemory(struct  elements_data  *elemdata)
{
	int  kelem;

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {
		if (elemdata->els[kelem].nlines > 0) {
			free(elemdata->els[kelem].wave_vacuum);
			free(elemdata->els[kelem].gf013divw3);
			free(elemdata->els[kelem].Eupdiv8);
		}
		free(elemdata->els[kelem].speclo);
		free(elemdata->els[kelem].spechi);
	}

	free(elemdata->wave);
	free(elemdata->lockfit);
	free(elemdata->corrfit);
	free(elemdata->dispfit);
	free(elemdata->fit_element_index);
	free(elemdata->neutral_index);

	free(elemdata->els);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Read reference star spectra file (which allocates memory as needed)
//-----------------------------------------------------------------------------------------

void  ReadStarSpectra(char *pathname, struct StarSpectraInfo *starspectra)
{
	int     kw, ks, nstars, nwave_star;
	double  star_specscale;
	char    textline[128];
	FILE   *starfile;


	//======== open the star spectra file for reading

	starfile = fopen(pathname, "rt");

	if (starfile == NULL) {
		printf(" Cannot open star spectra file %s for reading\n", pathname);
		exit(30);
	}


	//======== Read the header lines

	fscanf(starfile, "%[^=]= %d\n", textline, &nwave_star);

	fscanf(starfile, "%[^=]= %d\n", textline, &nstars);

	fscanf(starfile, "%[^=]= %lf\n", textline, &star_specscale);

	fscanf(starfile, "%[^\n]\n", textline);
	fscanf(starfile, "%[^\n]\n", textline);
	fscanf(starfile, "%[^\n]\n", textline);

	starspectra->nstars = nstars;

	starspectra->nwave = nwave_star;


	//======== Allocate memory for the file wavelengths, stars, and their spectra

	starspectra->wave = (double*)malloc(starspectra->nwave * sizeof(double));

	starspectra->star = (struct StarInfo*) malloc(starspectra->nstars * sizeof(struct StarInfo));

	for (ks = 0; ks < starspectra->nstars; ks++) {

		starspectra->star[ks].specval = (double*)malloc(starspectra->nwave * sizeof(double));

	}


	//======== Read the wavelengths line 

	fscanf(starfile, "%[^:]:", textline);

	for (kw = 0; kw < nwave_star; kw++)  fscanf(starfile, "%lf", &starspectra->wave[kw]);

	fscanf(starfile, "%[^\n]\n", textline);

	fscanf(starfile, "%[^\n]\n", textline);


	//======== Read each star line for its metadata and spectrum 

	for (ks = 0; ks < starspectra->nstars; ks++) {

		fscanf(starfile, "%ld %5s %lf %lf %lf %lf %lf %lf",
			&starspectra->star[ks].hip,
			&starspectra->star[ks].specType,
			&starspectra->star[ks].ra_deg,
			&starspectra->star[ks].dec_deg,
			&starspectra->star[ks].vmag,
			&starspectra->star[ks].bv,
			&starspectra->star[ks].vi,
			&starspectra->star[ks].bv0);

		starspectra->star[ks].userID = ks;  //... default until user reassigns

		for (kw = 0; kw < nwave_star; kw++)  fscanf(starfile, "%lf", &starspectra->star[ks].specval[kw]);

		fscanf(starfile, "\n");


		//-------- Scale the spectrum and compute the min and max spectral values for this star

		starspectra->star[ks].specmax = 0.0;
		starspectra->star[ks].specmin = 1.0e+30;

		for (kw = 0; kw < starspectra->nwave; kw++) {

			starspectra->star[ks].specval[kw] *= star_specscale;

			if (starspectra->star[ks].specmax < starspectra->star[ks].specval[kw]) {
				starspectra->star[ks].specmax = starspectra->star[ks].specval[kw];
			}

			if (starspectra->star[ks].specval[kw] > 0.0         &&
				starspectra->star[ks].specmin > starspectra->star[ks].specval[kw]) {
				starspectra->star[ks].specmin = starspectra->star[ks].specval[kw];
			}

		} //... end of wavelength loop

	} //... end of star lines reading loop


	//======== Close the star spectra files

	fclose(starfile);

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Get the star index into the StarSpectra structure based on a star's spectral type. 
//  If the spectral type is not found, returns -1.
//-----------------------------------------------------------------------------------------

int  GetStarIndexFromSPType(struct StarSpectraInfo *starspectra, char *specType)
{
	int   ks;

	for (ks = 0; ks < starspectra->nstars; ks++) {

		if (strstr(starspectra->star[ks].specType, specType) != NULL)  return(ks);

	}

	return(-1);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Get the star index into the StarSpectra structure based on a user's ID#. The user must
//  filled the StarSpectraInfo structure with the ID prior to being able to use this
//  function to pull a spectrum. Alternatively use GetStarIndexFromHIP to obtain the
//  star index given the Hipparcos catalog number (provided in the star spectra file so
//  prefill required)
//  If the userID is not found, returns -1.
//-----------------------------------------------------------------------------------------

int  GetStarIndexFromUserID(struct StarSpectraInfo *starspectra, int userID )
{
	int   ks;

	for (ks = 0; ks < starspectra->nstars; ks++) {

		if ( starspectra->star[ks].userID == userID )  return(ks);

	}

	return(-1);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Get the star index into the StarSpectra structure based on a star's Hipparcos number. 
//  This number is part of the star spectra file and is prefilled when reading the spectra.
//  If the Hipparcos number is not found, returns -1.
//-----------------------------------------------------------------------------------------

long  GetStarIndexFromHIP(struct StarSpectraInfo *starspectra, long hip)
{
	long   ks;

	for (ks = 0; ks < starspectra->nstars; ks++) {

		if (starspectra->star[ks].hip == hip)  return(ks);

	}

	return(-1);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Interpolate a specific "kstar" indexed reference star spectrum to the 
//  user's desired wavelengths. Assumes the reference star wavelengths are in  
//  monotonically increasing order.
//-----------------------------------------------------------------------------------------

void  InterpolateStarSpectrum(struct StarSpectraInfo *starspectra, long kstar, int nwave_desired, double *wavelengths_desired, double *spectrum_desired)
{
	int     kwave, kcurr;
	double  mslope;


	//============ Interpolate the spectrum to the desired output wavelengths

	kcurr = 0;

	for (kwave = 0; kwave < nwave_desired; kwave++) {

		//======== Are we below or above the star spectral wavelength band

		if (wavelengths_desired[kwave] < starspectra->wave[0] ||
			wavelengths_desired[kwave] > starspectra->wave[starspectra->nwave - 1]) {

			spectrum_desired[kwave] = starspectra->star[kstar].specmin;

			continue;
		}


		//======== Find the wavelength bin in the star spectrum vector 
		//           (assumed monotonically increasing in wavelength)

		while (kcurr < starspectra->nwave - 2) {

			if (wavelengths_desired[kwave] <= starspectra->wave[kcurr + 1]) break;

			kcurr++;
		}


		//======== Compute spectral slope in the wavelength bin and do linear interpolation

		mslope = (starspectra->star[kstar].specval[kcurr + 1] - starspectra->star[kstar].specval[kcurr])
			   / (starspectra->wave[kcurr + 1] - starspectra->wave[kcurr]);

		spectrum_desired[kwave] = starspectra->star[kstar].specval[kcurr] + mslope * (wavelengths_desired[kwave] - starspectra->wave[kcurr]);


	} //... End of interpolated star spectrum wavelength loop


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Free the memory for the StarSpectra structure
//-----------------------------------------------------------------------------------------

void  FreeMemoryStarSpectra(struct StarSpectraInfo *starspectra)
{
	int   ks;

	for (ks = 0; ks < starspectra->nstars; ks++)  free(starspectra->star[ks].specval);

	free(starspectra->star);

	free(starspectra->wave);

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Set user adjustable values in the elemdata structure to their starting defaults as
//      given in the input configuration file. These include parameters such as
//      sigma, warm and hot temperatures, electron density, hot-to-warm, airmass factor
//-----------------------------------------------------------------------------------------

void   AdjustableParametersDefaults(struct specconfiguration *spconfig, struct elements_data *elemdata)
{

	elemdata->sigma0   = spconfig->nominal_sigma0;
	elemdata->Tlo      = spconfig->nominal_lo_exc_temp;
	elemdata->Thi      = spconfig->nominal_hi_exc_temp;
	elemdata->hot2warm = spconfig->default_hot2warm;
	elemdata->ne_jones = spconfig->default_ne;
	elemdata->ne_iter  = spconfig->default_ne;
	elemdata->Xfactor  = 1.0;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Extinction Approximation
//
//  Computes an approximate extinction based on zenith angle (input as altitude elevation
//  above the horizon in degrees "elev_deg") and the site's height "hobs" above the
//  WGS84 surface in km.
//
//  Inputs:   elev_deg, hobs
//  Outputs:  extinction
//-----------------------------------------------------------------------------------------

double  ExtinctionApprox( double elev_deg, double hobs_km, double wavelength_nm )
{
double  k1, extinct, cosz;


    cosz = cos( ( 90.0 - elev_deg ) * 3.141592654 / 180.0 );

    k1 = 1.87e+11 * ( 0.1451 * exp(-hobs_km /7.996) + 0.120 * exp(-hobs_km /1.5) + 0.016 ) * pow( wavelength_nm, -4.18 );

    extinct = exp ( -k1 / ( cosz + 0.025 * exp( -11.0 * cosz) ) );

	return( extinct );

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Airmass
//
//  Computes the airmass through the US standard atmospheric density given the observer's
//  look direction altitude (elevation) above the horizon in degrees "elev_deg" and the
//  site's radius "rhobs" in km equal to WGS84 height + WGS84 Earth radius.
//
//  Inputs:   elev_deg, rhobs
//  Outputs:  airmass
//-----------------------------------------------------------------------------------------

double  Airmass( double elev_deg, double rhobs )
{
double  x;

		x = RhoIntegral( elev_deg, rhobs ) / RhoIntegral( 90.0, rhobs );
        
        return( x );
}     



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  RhoIntegral
//
//  Computes the line integral through the US standard atmospheric density given the 
//  observer's look direction altitude (elevation) above the horizon in degrees "elev_deg"
//   and the site's radius "rhobs" in km equal to WGS84 height + WGS84 Earth radius.   
//
//  Inputs:   elev_deg, rhobs
//  Outputs:  integrated_density
//-----------------------------------------------------------------------------------------

double  RhoIntegral( double elev_deg, double rhobs )
{
double  elev, sinelev, zmax, lzmax, lz, z, rearth;
double  logrho, integrated_density;

        
        rearth = 6378.137;  //... km

        elev = elev_deg * 3.141592654 / 180.0;

        sinelev = sin(elev);
        
        zmax = 50.0;
        
        if( elev_deg == 90.0 )  lzmax = 50.0;
        else                    lzmax = (rearth + zmax) 
                                      * sin( acos( rhobs*cos(elev)
                                      / (rearth + zmax) ) - elev )
                                      / cos(elev);
                                      
        integrated_density = 0.0;

        for( lz=0.0; lz<lzmax; lz+=0.01 )  {
        
             z = sqrt( rhobs*rhobs + lz*lz + 2.0*lz*rhobs*sinelev )
               - rearth;
                            
             logrho =       0.0829530
                    - z * ( 0.0346810 
                    + z * ( 0.0014865
                    - z * ( 1.5604e-5
                    + z * ( 3.0649e-7
                    - z * ( 4.7987e-9 )))));
               
             integrated_density += 0.01 * pow( 10.0, logrho );
        }
        
        return( integrated_density );
}     


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Fit the extinction coefficients for Scale_factor, Rayleigh+O2, Ozone, Dust, and Water
//  for the given extinction spectrum. Fit coefficients are placed in the 
//  speccalibration structure.
//-----------------------------------------------------------------------------------------

void   FitExtinctionCoefs( int      nwave,
	                       double  *wavelength_nm, 
	                       double  *extinction_spectrum,
						   struct   speccalibration  *spcalib )
{
int      k, krow, kcol, kwave, kiter, khigh;
double   tauR, tauO2A, tauO2B, tauRO2, tauO3, tauD, tauW, logextn, sum, coefhigh;

#define  NUMCOEFS  5   //... Scale_factor, Rayleigh+O2, Ozone, Dust, Water

double   *tpe;
double  **tpt;  // svd of tpt = USV'
double  **US;
double  **V;
double   *Ssq;  // squared singular values
double   *work;
double   *coefs;


    //========== Allocate memory for the 1D and 2D arrays used for the SVD LMS solution

    tpe   = (double*)malloc(NUMCOEFS * sizeof(double));
    Ssq   = (double*)malloc(NUMCOEFS * sizeof(double));
	work  = (double*)malloc(NUMCOEFS * sizeof(double));
	coefs = (double*)malloc(NUMCOEFS * sizeof(double));

    tpt = (double**)malloc((NUMCOEFS * sizeof(double*)) + (NUMCOEFS * NUMCOEFS * sizeof(double)));
    for (k = 0; k < NUMCOEFS; ++k)  tpt[k] = (double*)(tpt + NUMCOEFS) + k * NUMCOEFS;

    US = (double**)malloc((NUMCOEFS * sizeof(double*)) + (NUMCOEFS * NUMCOEFS * sizeof(double)));
    for (k = 0; k < NUMCOEFS; ++k)  US[k] = (double*)(US + NUMCOEFS) + k * NUMCOEFS;

    V = (double**)malloc((NUMCOEFS * sizeof(double*)) + (NUMCOEFS * NUMCOEFS * sizeof(double)));
    for (k = 0; k < NUMCOEFS; ++k)  V[k] = (double*)(V + NUMCOEFS) + k * NUMCOEFS;


    //========== clear the T'E vector and T'T matrix

    for( krow=0; krow<NUMCOEFS; krow++ )  {

		tpe[krow] = 0.0;

		for( kcol=0; kcol<NUMCOEFS; kcol++ )  tpt[krow][kcol] = 0.0;

	}


	//========== fill the T'E vector and T'T matrix

	for( kwave=0; kwave<nwave; kwave++ )  {

		tauR   = TauOpticalDepth( RAYLEIGH, wavelength_nm[kwave] );

		tauO2A = TauOpticalDepth( O2BAND1,  wavelength_nm[kwave] );

		tauO2B = TauOpticalDepth( O2BAND2,  wavelength_nm[kwave] );

		tauO3  = TauOpticalDepth( OZONE,    wavelength_nm[kwave] );

		tauD   = TauOpticalDepth( DUST,     wavelength_nm[kwave] );

		tauW   = TauOpticalDepth( WATER,    wavelength_nm[kwave] );


		tauRO2 = tauR + tauO2A + tauO2B;


        if( extinction_spectrum[kwave] <= 0.0  ||  
			extinction_spectrum[kwave] >= 1.0     ) continue;

		logextn = log( extinction_spectrum[kwave] );

		tpe[0] += logextn;
		tpe[1] += logextn * tauRO2;
		tpe[2] += logextn * tauO3;
		tpe[3] += logextn * tauD;
		tpe[4] += logextn * tauW;

		tpt[0][0] += 1.0 * 1.0;
		tpt[0][1] += 1.0 * tauRO2;
		tpt[0][2] += 1.0 * tauO3;
		tpt[0][3] += 1.0 * tauD;
		tpt[0][4] += 1.0 * tauW;
	
		tpt[1][0] += tauRO2 * 1.0;
		tpt[1][1] += tauRO2 * tauRO2;
		tpt[1][2] += tauRO2 * tauO3;
		tpt[1][3] += tauRO2 * tauD;
		tpt[1][4] += tauRO2 * tauW;
	
		tpt[2][0] += tauO3 * 1.0;
		tpt[2][1] += tauO3 * tauRO2;
		tpt[2][2] += tauO3 * tauO3;
		tpt[2][3] += tauO3 * tauD;
		tpt[2][4] += tauO3 * tauW;
	
		tpt[3][0] += tauD * 1.0;
		tpt[3][1] += tauD * tauRO2;
		tpt[3][2] += tauD * tauO3;
		tpt[3][3] += tauD * tauD;	
		tpt[3][4] += tauD * tauW;	

		tpt[4][0] += tauW * 1.0;
		tpt[4][1] += tauW * tauRO2;
		tpt[4][2] += tauW * tauO3;
		tpt[4][3] += tauW * tauD;	
		tpt[4][4] += tauW * tauW;	

	}


	//========== Iterate through until all remaining coefs are absorption lines (except scale_factor)

	for( kiter=0; kiter<NUMCOEFS-1; kiter++ )  {

        //-------- compute the solution [T'T]^-1 T'E  =  U S V' T' E

		svd_NxN( &tpt[0][0], &US[0][0], &Ssq[0], &V[0][0], NUMCOEFS, NUMCOEFS ); 

		for( kcol=0; kcol<NUMCOEFS; kcol++ )  {

			 sum = 0.0;

		     for( krow=0; krow<NUMCOEFS; krow++ )  {
				  sum += US[krow][kcol] * tpe[krow];
			 }

			 if( Ssq[kcol] > 0.0 )  work[kcol] = sum  / Ssq[kcol];
			 else                   work[kcol] = 0.0;
		}


		for( krow=0; krow<NUMCOEFS; krow++ )  {

			 sum = 0.0;

			 for( kcol=0; kcol<NUMCOEFS; kcol++ )  {
				 sum += V[krow][kcol] * work[kcol];
			 }
			 coefs[krow] = sum;
		}

		//-------- now check for any emission coefs > 0

		khigh = -1;

		coefhigh = -1.0e+30;

		for( krow=1; krow<NUMCOEFS; krow++ )  {

			if( coefhigh < coefs[krow] )  {
				coefhigh = coefs[krow];
				khigh = krow;
			}
		}

		if( coefhigh > 0.0 )  { //... found an emission coef

			//.... adjust corresponding row and column to isolate equations

			for( krow=0; krow<NUMCOEFS; krow++ )  {
				tpt[krow][khigh] = 0.0;
				tpt[khigh][krow] = 0.0;
			}
			tpt[khigh][khigh] = 1.0;
			tpe[khigh] = 0.0;

		}
		else  khigh = -1;


		//..... end iteration if all absorption coefs

		if( khigh == -1 )  break;

	}
		
		
    //========== assign the fit extinction coefficients

	spcalib->extn_scale    = exp( coefs[0] ); 
	spcalib->extn_rayleigh = -coefs[1]; 
	spcalib->extn_oxygen2  = -coefs[1];   
	spcalib->extn_ozone    = -coefs[2];      
	spcalib->extn_dust     = -coefs[3];   
	spcalib->extn_water    = -coefs[4];


	//======== Free temporary memory

	free(tpe);
	free(Ssq);
	free(work);
	free(coefs);
	free(tpt);
	free(US);
	free(V);

						   
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Compute the optical depth for the absorption type: 
//       Rayleigh, O2 bands 1 & 2, Ozone, Dust, Water
//-----------------------------------------------------------------------------------------

double  TauOpticalDepth( int ktype, double wavelength_nm )
{
double  tau;

#define  NWATERLINES  21
int      kthline;
double   cen_wid_amp[NWATERLINES][3] = {  // water absorption center wavelength, line width, amplitude
                                 {  648.0,    6.0, 0.0001 },
                                 {  698.0,    6.0, 0.0003 },
                                 {  728.0,    6.0, 0.0007 },
                                 {  818.0,    6.0, 0.0014 },
                                 {  855.0,   10.0, 0.0010 },
                                 {  902.0,   12.0, 0.0011 },
                                 {  942.0,   21.0, 0.0040 },
                                 {  995.0,   41.0, 0.0006 },
                                 { 1064.0,   22.0, 0.0020 },
                                 { 1127.0,   18.0, 0.0059 },
                                 { 1146.0,   16.0, 0.0031 },
                                 { 1207.0,   10.0, 0.0022 },
                                 { 1222.0,   10.0, 0.0015 },
                                 { 1267.0,   10.0, 0.0220 },
                                 { 1363.0,   11.0, 0.0290 },
                                 { 1382.0,   18.0, 0.0230 },
                                 { 1400.0,   17.0, 0.0180 },
                                 { 1436.0,   15.0, 0.0250 },
                                 { 1472.0,   10.0, 0.0040 },
                                 { 1538.0,    8.0, 0.0024 },
                                 { 1575.0,    9.0, 0.0130 }
                                       };


     switch( ktype )
	 {
	     case RAYLEIGH:  
		     tau = pow( wavelength_nm, -4.0 ) / 5.7e-10;
		     break;

	     case O2BAND1:   
		     tau = AbsorptionLine( wavelength_nm, 688.0, 8.0, 0.02 );
		     break;

	     case O2BAND2: 
		     tau = AbsorptionLine( wavelength_nm, 762.0, 8.0, 0.09 );
		     break;

	     case OZONE:  
	    	 tau = AbsorptionLine( wavelength_nm, 590.0, 110.0, 0.04 );
		     break;

	     case DUST: 
		     tau = 4.348 / wavelength_nm;
		     break;

	     case WATER: 
		     tau = 0.0;
			 for( kthline=0; kthline<NWATERLINES; kthline++ )  {
				  tau += AbsorptionLine( wavelength_nm, cen_wid_amp[kthline][0], cen_wid_amp[kthline][1], cen_wid_amp[kthline][2] );
			 }
		     break;

	     default:
		     tau = 0.0; 
		     break;
	 }

	 return( tau );	 


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the extinction "extinction_model" at "nwave" user specified wavelengths 
//   "wavelength_nm" given the extinction coefficients in the spcalib structure, site radius
//   "rhobs" = WGS84 height + Earth radius, and the look altitude (elevation) in degrees.
//-----------------------------------------------------------------------------------------

void   ExtinctionModel( int      nwave,
	                    double  *wavelength_nm, 
	                    double  *extinction_model,
						double   elevation_deg,   //... 90 deg gives unity airmass
						double   rhobs,  //... site H + Rearth
						struct   speccalibration  *spcalib )
{
int      kwave;
double   E, X;


		//========== Get the airmass and compute extinction at each wavelength

        X = Airmass( elevation_deg, rhobs );

		for( kwave=0; kwave<nwave; kwave++ )  {

			E = spcalib->extn_scale    * exp(
		      - spcalib->extn_rayleigh * TauOpticalDepth( RAYLEIGH, wavelength_nm[kwave] )
			  - spcalib->extn_oxygen2  * TauOpticalDepth( O2BAND1,  wavelength_nm[kwave] )
			  - spcalib->extn_oxygen2  * TauOpticalDepth( O2BAND2,  wavelength_nm[kwave] )
			  - spcalib->extn_ozone    * TauOpticalDepth( OZONE,    wavelength_nm[kwave] )
			  - spcalib->extn_dust     * TauOpticalDepth( DUST,     wavelength_nm[kwave] )
			  - spcalib->extn_water    * TauOpticalDepth( WATER,    wavelength_nm[kwave] ) );

			extinction_model[kwave] = pow( E, X );

		}

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Compute absorption line response at a wavelength "wavelength_nm" off from the line's
//  center wavelength "center_nm", line width "width_nm" and "depth_fraction".
//  Wavelengths are in nanometers (nm).
//-----------------------------------------------------------------------------------------

double   AbsorptionLine( double  wavelength_nm, 
	                     double  center_nm,
					     double  width_nm,
					     double  depth_fraction )
{

double  tau = depth_fraction * pow( width_nm/2.0, 2.0 ) / ( pow( center_nm - wavelength_nm, 2.0 ) + pow( width_nm/2.0, 2.0 ) );

return( tau );

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    Allocate memory for the measured and fit meteor spectra and column density
//-----------------------------------------------------------------------------------------

int      AllocMemorySpectrum( int                       nwavelengths,
	                          int                       nelements, 
	                          struct meas_fit_spectra  *meteor_spec )
{

	meteor_spec->meas_spectrum  = (double*)malloc(nwavelengths * sizeof(double));
	meteor_spec->fit_spectrum   = (double*)malloc(nwavelengths * sizeof(double));
	meteor_spec->columndensity  = (double*)malloc(   nelements * sizeof(double));

	if (meteor_spec->meas_spectrum == NULL ||
		meteor_spec->fit_spectrum  == NULL ||
		meteor_spec->columndensity == NULL)  {
		printf(" Memory not allocated for meteorspectrum structure in AllocMemorySpectra\n");
		return(3);
	}

	return(0);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    Free memory for the measured and fit meteor spectra and column density
//-----------------------------------------------------------------------------------------

void      FreeMemorySpectrum(struct meas_fit_spectra    *meteor_spec)
{
	free( meteor_spec->meas_spectrum );
	free( meteor_spec->fit_spectrum  );
	free( meteor_spec->columndensity );
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the approximate plasma radius at a height above the surface in km.
//-----------------------------------------------------------------------------------------

double  PlasmaRadius(double height_km)
{
	double  plasma_radius_meters, log10diameter_meters, H100;


	H100 = height_km / 100.0;

	log10diameter_meters = + 3.1516   * pow(H100, 8)
		                   - 29.3745  * pow(H100, 7)
		                   + 111.0692 * pow(H100, 6)
		                   - 217.0351 * pow(H100, 5)
		                   + 230.3760 * pow(H100, 4)
		                   - 127.2633 * pow(H100, 3)
		                   + 30.4579  * pow(H100, 2)
		                   + 5.1829   * H100
		                   - 5.7075;

	plasma_radius_meters = 0.5 * pow(10.0, log10diameter_meters);

	////printf(" H = %lf   r = %lf \n", H_km, plasma_radius_meters );

	return(plasma_radius_meters);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the warm and hot plasma volumes any time one of these changes:
//       Height above the surface in km
//       Range to the meteor in km
//       Hot-to-Warm ratio
//       Approach angle in radians between the look angle and the velocity vector (the same
//             as the angle between the radiant and look equatorial positions).
//
//   Computed products are saved in the elements_data structure.
//-----------------------------------------------------------------------------------------

void    PlasmaVolumes(double height_km, double range_km, double approach_angle, double hot2warm, struct elements_data *elemdata)
{

	//--------- Compute the plasma radius and volumes, convert range to meters

	elemdata->plasma_radius_meters = PlasmaRadius(height_km);

	elemdata->range_meteor_meters = 1000.0 * range_km;


	elemdata->VolumeDepth = 2.0 * elemdata->plasma_radius_meters * sin(approach_angle);

	elemdata->VolumeTlo = (4.0 / 3.0) * 3.141592654 * pow(elemdata->plasma_radius_meters, 3);  //... m^3

	elemdata->VolumeThi = elemdata->VolumeTlo * hot2warm;  //... m^3


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Find the index into the table of elements (neutrals and ions) given the atomicNumber
//   of the element, the ionization state (neutral atom=0, first ionization atom=1, neutral
//   molecule=50, first ionization molecule=51). Returns -1 if not found in table.
//-----------------------------------------------------------------------------------------

int    GetElementIndex(int atomicNumber, int ionization, struct elements_data *elemdata)
{
	int kelem;
	double  elementcode;


	elementcode = (double)atomicNumber + (double)ionization / 100.0;

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {
		if (elementcode == elemdata->els[kelem].elementcode)  return(kelem);
	}

	return(-1);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Find the primary neutral element in order of Mg-I, Fe-I, then Na-I that is 
//   actively being fit.
//-----------------------------------------------------------------------------------------

int   GetPrimaryNeutralElement(struct elements_data *elemdata)
{
	int    kelem;


	//========== First test for neutral Mg as an active element

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		if (elemdata->els[kelem].elementcode == 12.0  &&  elemdata->els[kelem].user_fitflag != FITLESS)  return(kelem);

	}


	//========== Second test for Fe as an active element

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		if (elemdata->els[kelem].elementcode == 26.0  &&  elemdata->els[kelem].user_fitflag != FITLESS)  return(kelem);

	}


	//========== Third test for Na as an active element

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		if (elemdata->els[kelem].elementcode == 11.0  &&  elemdata->els[kelem].user_fitflag != FITLESS)  return(kelem);

	}


	//========== Additonal tests would go here

	printf("\a No active element for Mg-I, Fe-I, or Na-I - Densities will not be computed\n");

	return(-1);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the Jones electron density given the primary neutral element
//-----------------------------------------------------------------------------------------

double  JonesElectronDensity(struct elements_data *elemdata, int kneutral_primary)
{
	int     kneutral, kelem;
	double  ne, n_neutrals_primary, n_atoms_primary, n_atoms_cosmic, n_ions_cosmic;


	//========== Number density of neutrals for the primary element

	n_neutrals_primary = elemdata->els[kneutral_primary].N_warm / elemdata->VolumeDepth;
	printf("Number of neutral = %e\n", n_neutrals_primary);


	//========== Number density of atoms (neutrals + ions) for the primary element based on Jones fraction

	n_atoms_primary = n_neutrals_primary / (1.0 - elemdata->els[kneutral_primary].beta_jones);
	printf("Number of atoms = %e\n", n_atoms_primary);
	printf("Beta Jones = %f\n", elemdata->els[kneutral_primary].beta_jones);


	//========== Compute sum over all elements to get number density of electrons = number density of ions

	ne = 0.0;


	//========== Loop over all NEUTRAL elements (so we don't double count with the ion elements)

	for (kneutral = 0; kneutral < elemdata->nneutrals; kneutral++) {

		kelem = elemdata->neutral_index[kneutral];

		//-------- Number density of atoms total for this element based on cosmic abundance ratio re primary element

		n_atoms_cosmic = n_atoms_primary * elemdata->els[kelem].init_abundance / elemdata->els[kneutral_primary].init_abundance;
		// printf("Number of atoms based on cosmic abundance = %e\n", n_atoms_cosmic);
		//-------- Add in number of ions for this element based on Jones fraction beta

		n_ions_cosmic = n_atoms_cosmic * elemdata->els[kelem].beta_jones;
		// printf("Number of ions_cosmic = %e   beta_jones = %e\n", n_ions_cosmic, elemdata->els[kelem].beta_jones);

		ne += n_ions_cosmic;
		// printf("ne = %e\n", ne);

	}

	elemdata->ne_jones = ne;
	printf("ne_jones = %e\n", elemdata->ne_jones);

	return(ne);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Iteratively compute the electron density given the primary neutral element and
//   a procedure developed by P. Jenniskens.
//-----------------------------------------------------------------------------------------

double    IterativeElectronDensity(struct elements_data *elemdata, int kneutral_primary, double ne_guess)
{
	int     kloop, nloops, kneutral, kelem, kion;
	double  ne_sum, ne, nratio, n_neutrals_primary, n_atoms_primary, n_neutrals, n_atoms, n_ions;


	//========== Iteratively loop until the electron number density converges, uses Jones ne for initial guess
	printf("Iteratively computing electron density...\n");
	nloops = 50;

	ne = ne_guess;

	for (kloop = 0; kloop < nloops; kloop++) {
		// printf("Iteration #%d\n", kloop+1);

		//-------- Number density of warm neutrals fit as column density for the primary reference element
		// printf("N_warm: %f\n", elemdata->els[kneutral_primary].N_warm);
		// printf("Volume Depth %f\n", elemdata->VolumeDepth);
		n_neutrals_primary = elemdata->els[kneutral_primary].N_warm / elemdata->VolumeDepth;


		//-------- Number density of atoms (neutral + ions) for the primary element using Saha's equation

		kion = elemdata->els[kneutral_primary].ionindex;

		nratio = 2.0 * 2.4146874e+21 * pow(elemdata->Tlo, 1.5)
			* exp(-elemdata->els[kneutral_primary].ionenergy / 8.617343e-5 / elemdata->Tlo)
			* elemdata->els[kion].partfuncTlo / elemdata->els[kneutral_primary].partfuncTlo / ne;

		n_atoms_primary = n_neutrals_primary * (1.0 + nratio);
		// printf("nratio: %e\nn_atoms: %f\n", nratio, n_atoms_primary);

		//-------- Sum the contributing #ions = #electrons for all elements

		ne_sum = 0.0;

		// printf("%le %d %le %le\n", n_neutrals_primary, kion, nratio, n_atoms_primary);

		//-------- Loop over all NEUTRAL elements to get column densities and get number of atoms

		for (kneutral = 0; kneutral < elemdata->nneutrals; kneutral++) {

			kelem = elemdata->neutral_index[kneutral];


			//....... Get the element index of the corresponding ion for this neutral

			kion = elemdata->els[kelem].ionindex;

			if (kion < 0)  continue;  //... No corresponding ion in the table (do not include in ne)


			//....... Use the Saha equation to get ratio of ions to neutrals for this element

			nratio = 2.0 * 2.4146874e+21 * pow(elemdata->Tlo, 1.5)
				* exp(-elemdata->els[kelem].ionenergy / 8.617343e-5 / elemdata->Tlo)
				* elemdata->els[kion].partfuncTlo / elemdata->els[kelem].partfuncTlo / ne;


			//....... EITHER the element was fit and has a known neutral column density, in which
			//                case divide by the column depth for number density and scale by   
			//                the Saha ratio to get the number of ions for this element,
			//        OR use the cosmic abundance ratio and the primary element's atom number density 
			//               to obtain the total number of atoms (neutrals + ions) for this element,  
			//               and then use the Saha ratio to obtain the number of ions for this element.

			if (elemdata->els[kelem].user_fitflag != FITLESS) {

				n_neutrals = elemdata->els[kelem].N_warm / elemdata->VolumeDepth;  //... from the fit to a warm neutral element

				n_ions = n_neutrals * nratio;

			}

			else {

				n_atoms = n_atoms_primary * elemdata->els[kelem].init_abundance / elemdata->els[kneutral_primary].init_abundance;

				n_ions = n_atoms * nratio / (1.0 + nratio);

			}


			//....... Add in number of ions for this element to the running sum of electrons

			ne_sum += n_ions;


		} //... end of elements loop


		//-------- Test for convergence such that the relative change is less than 0.1%

		if (fabs((double)ne - (double)ne_sum) / ((double)ne + 1.0e-30) < 0.001)  break;


		//-------- Next guess for ne is based on an average of the previous and current estimate

		ne = (ne + ne_sum) / 2;

	} //... end of iterative loop


	elemdata->ne_iter = ne;

	elemdata->numiter = kloop;

	return(ne);

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the column densities and number of atoms from a given fit.
//-----------------------------------------------------------------------------------------

void    ColumnDensities_NumberAtoms(struct elements_data *elemdata, double ne)
{
	int     kneutral, kelem, kion, kneu;
	double  N_hot_total, number_density;


	//========== Make sure any previous fit values have been cleared

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		elemdata->els[kelem].N_warm_total = 0.0;  // warm total column density (Neutrals + Ions)
		elemdata->els[kelem].number_atoms = 0.0;  // warm number of atoms

	}


	//========== Loop over all NEUTRAL elements to get column densities and get number of atoms

	for (kneutral = 0; kneutral < elemdata->nneutrals; kneutral++) {

		kelem = elemdata->neutral_index[kneutral];

		//-------- Only perform these operations on fit elements

		if (elemdata->els[kelem].user_fitflag != FITLESS) {

			//------- Compute column densities ONLY if a warm element

			if (elemdata->els[kelem].warmhot == WARM) {

				//....... Get the element index of the corresponding neutral and ion

				kneu = kelem;
				kion = elemdata->els[kelem].ionindex;


				//....... If an element has a neutral in the table...

				if (kneu >= 0) {

					//.... Compute the ratio of ions to neutrals via Saha's equation for the WARM component
					//        and via Jone's Beta for the HOT component. If there is no ion for this neutral
					//        then set these ratios to zero. Store in both the neutral and ion data structure
					//        indices kneu and kion respectively.

					if (kion >= 0) {

						elemdata->els[kneu].ions2neut_warm = 2.0 * 2.4146874e+21 * pow(elemdata->Tlo, 1.5)
							* exp(-elemdata->els[kneu].ionenergy / 8.617343e-5 / elemdata->Tlo)
							* elemdata->els[kion].partfuncTlo / elemdata->els[kneu].partfuncTlo
							/ ne;

						elemdata->els[kneu].ions2neut_hot = elemdata->els[kneu].beta_jones / (1.0 - elemdata->els[kneu].beta_jones);

					}

					else {

						elemdata->els[kneu].ions2neut_warm = 0.0;

						elemdata->els[kneu].ions2neut_hot = 0.0;

					}


					//########################################################################################
					//    Note: These sections assume we have a neutral warm column density N_warm
					//########################################################################################

					//========================  WARM Plasma Component  =======================================

					//  Given the neutral warm column density compute both the ion column density and the total
					//    atom (ion + neutral) column density

					elemdata->els[kion].N_warm = elemdata->els[kneu].N_warm * elemdata->els[kneu].ions2neut_warm;

					elemdata->els[kneu].N_warm_total = elemdata->els[kneu].N_warm + elemdata->els[kion].N_warm;


					//.... Get number density and the total number of atoms for the warm component

					number_density = elemdata->els[kneu].N_warm_total / elemdata->VolumeDepth;

					elemdata->els[kneu].number_atoms = number_density * elemdata->VolumeTlo;


					//========================  HOT Plasma Component  =======================================

					//.... Compute the element's hot plasma neutral and ion column densities by assuming
					//        the abundance is the same for warm and hot plasma so the total atoms column
					//        density between hot and warm scale only by the hot-to-warm volume ratio.

					N_hot_total = elemdata->els[kneu].N_warm_total * elemdata->hot2warm;

					elemdata->els[kneu].N_hot = N_hot_total / (1.0 + elemdata->els[kneu].ions2neut_hot);

					elemdata->els[kion].N_hot = elemdata->els[kneu].N_hot * elemdata->els[kneu].ions2neut_hot;

					//.... Print statements for de-bubbing purposes

					//printf("neu_ion %d %d %le %le %le %le %le %le %le\n", kneu, kion, 
					//elemdata->els[kion].N_warm,
					//elemdata->els[kneu].N_warm_total,
					//number_density,
					//elemdata->els[kneu].number_atoms,
					//N_hot_total,
					//elemdata->els[kneu].N_hot,
					//elemdata->els[kion].N_hot);

					//.... End of de-bugging print statements


				} //... end IF there exists both a neutral and ion

			} //... end of IF a warm neutral element


			//------- Compute fluxes for neutral/ion and warm/hot ONLY if a neutral hot element (e.g. Hydrogen)

			if (elemdata->els[kelem].ioncode != 2 && elemdata->els[kelem].warmhot == HOT) {
				printf(" ===> ERROR: ColDensities_NumberAtoms has a HOT neutral %s not implemented\n", elemdata->els[kelem].element_string);
			}


		} //... end of IF an actively fit element

	} //... end of elements loop


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Clear a specific element (index kelem) to no column density and zero atoms.
//-----------------------------------------------------------------------------------------

void  ResetOneElementAbundance( int kelem, double Vinfinity, struct elements_data *elemdata )
{
	double  cvterm;

    //========== Initialize specific element to zero column densities and no atoms

	elemdata->els[kelem].N_warm       = 0.0;  // warm column density (Neutrals + Ions)
	elemdata->els[kelem].N_warm_total = 0.0;  // warm total column density (Neutrals + Ions)
	elemdata->els[kelem].N_hot        = 0.0;  // hot column density (Neutrals + Ions)
	elemdata->els[kelem].number_atoms = 0.0;  // warm number of atoms


	//========== Initialize the Jones' beta value

	cvterm = elemdata->els[kelem].c * pow(Vinfinity - elemdata->els[kelem].Vo, 2) * pow(Vinfinity, 0.8);
	//printf("CVTERM = %e\n", cvterm);
	elemdata->els[kelem].beta_jones = cvterm / (1.0 + cvterm);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Set all elements to no column density, zero atoms, compute Jones' beta.
//-----------------------------------------------------------------------------------------

void  ResetAllElementAbundances(double Vinfinity, struct elements_data *elemdata)
{
	int     kelem;

	for (kelem = 0; kelem < elemdata->nelements; kelem++) {

		elemdata->els[kelem].user_fitflag = FITLESS;

		ResetOneElementAbundance(kelem, Vinfinity, elemdata);

	}

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the relative abundances for all the elements
//-----------------------------------------------------------------------------------------

void  ComputeRelativeAbundance( struct  elements_data  *elemdata )
{
int     kneutral, kelem;
struct  element_lines_spectrum  *els_neu, *els_ref;


    //========== Get the reference neutral elements total warm column density

	els_ref = &elemdata->els[elemdata->kelem_ref];


    //========== Loop over all NEUTRAL elements to get abundances

	for( kneutral=0; kneutral<elemdata->nneutrals; kneutral++ )  {

		 kelem = elemdata->neutral_index[kneutral];

		 els_neu = &elemdata->els[kelem];


		 //-------- Check for no reference fit, in which case set the abundance to zero

		 els_neu->abundance = 0.0;

		 if( els_ref->user_fitflag == FITLESS )  continue;
			 

	     //-------- Check if element NOT fit for neutral, then set to solar adundance, 
		 //           otherwise set to ratio of total warm column densities

		 if( els_neu->user_fitflag == FITLESS )  els_neu->abundance = els_neu->init_abundance / els_ref->init_abundance;
		 else                                    els_neu->abundance = els_neu->N_warm_total   / els_ref->N_warm_total; 

	} //... end of element loop


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the model spectrum due to Fe only.
//   Needs to be called if sigma0 or T (Tlo) changes from last Fe only spectrum calculation.
//-----------------------------------------------------------------------------------------

void  IronOnlySpectrumT( double T, 
                         struct  elements_data  *elemdata,
						 double                 *responsivity_1st, 
						 double                 *responsivity_2nd, 
					     double                  grating_area_scaling,
						 double                 *extinction )
{
int     kelem, kline, kwave, kwave_start, kline_start, kline_finis, kpf;
int     kmin, kmax, kwave2nd;
double  invsigmasq_1st, invsigma2pi_1st;
double  invsigmasq_2nd, invsigma2pi_2nd;
double  sigma, sigmahalf, dlimit, dlambda, expterm, pf;
double  maxspec, sumterm, RespExtn, solid_angle;
struct  element_lines_spectrum  *els;


    //========== Compute broadening terms

    sigma = elemdata->sigma0 / 2.3548;

	sigmahalf = sigma / 2.0;

    invsigmasq_1st = -0.5 /  ( sigma * sigma );

	invsigma2pi_1st = +1.0 / ( sigma * sqrt( 2.0 * 3.141592654 ) );

    invsigmasq_2nd = -0.5 /  ( sigmahalf * sigmahalf );

	invsigma2pi_2nd = +1.0 / ( sigmahalf * sqrt( 2.0 * 3.141592654 ) );


	//========== Compute lambda cutoff range when contribution = 10^-6

	dlimit = sigma * sqrt( -2.0 * log( 1.0e-6 ) );

	kmin = (int)( elemdata->wave[0] / ( elemdata->wave[1] - elemdata->wave[0] ) );


    //========== Loop through each element and build spectra

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

		//......... only do the Fe lines

		if( elemdata->els[kelem].elementcode != 26.0 )  continue;

 		 els = &elemdata->els[kelem];

		 //-------- clear the spectrum and check if user active

		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  els->speclo[kwave] = 1.0e-60;

		 if( els->user_fitflag == FITLESS )  continue;

		 //printf("------------------------------------------  %s\n", els->element_string );

		 //-------- Compute partition function at T using piecewise linear interpolation

	     kpf = 0;

	     while( T > els->Tpf[kpf+1]  &&  kpf < els->Npf-2 )  kpf++;

		 pf = els->partfunc[kpf]  +  ( els->partfunc[kpf+1] - els->partfunc[kpf] ) * ( T - els->Tpf[kpf] ) / ( els->Tpf[kpf+1] - els->Tpf[kpf] );


		 //-------- Loop over the emission lines for this element
		 //         NOTE: line wavelengths are assumed sorted in ascending order
		 //               and this was checked when they were read in. This
		 //               allows the use of an increasing wavelength sliding 
		 //               window of output wavelengths that bracket the line 
		 //               wavelengths as we progress through the lines.

		 //-------- get the line start and stop indices

		 kline_start = 0;

		 while( elemdata->wave[0] - els->wave_vacuum[kline_start] > dlimit  &&  kline_start < els->nlines-1 )  kline_start++;

		 kline_finis = els->nlines-1;

		 while( els->wave_vacuum[kline_finis] - elemdata->wave[elemdata->nwave-1] > dlimit  &&  kline_finis > 0 )  kline_finis--;


		 kwave_start = 0;

		 for( kline=kline_start; kline<=kline_finis; kline++ )  {

			  //- - - - - Shift the minimum output wavelength upwards, if it is > dlimit below the line wavelength

			  while( els->wave_vacuum[kline] - elemdata->wave[kwave_start] > dlimit  &&  kwave_start < elemdata->nwave-1 )  kwave_start++;


			  //- - - - - Sum the spectra over the narrow contributing wavelength band

			  for( kwave=kwave_start; kwave<elemdata->nwave; kwave++ )  {

				   dlambda = elemdata->wave[kwave] - els->wave_vacuum[kline];

				   if( dlambda > dlimit )  break;


				   //..... Compute this for 1st order

				   expterm = invsigmasq_1st * dlambda * dlambda + els->Eupdiv8[kline] / T;

				   sumterm = els->gf013divw3[kline] * exp( expterm );

				   RespExtn = responsivity_1st[kwave] * extinction[kwave];

				   els->speclo[kwave] += invsigma2pi_1st * RespExtn * sumterm;


				   //..... Compute this for 2nd order adding half the term because of 2:1 wavelength mapping 1st to 2nd order

				   kwave2nd = ( kwave - kmin ) / 2;

				   expterm = invsigmasq_2nd * dlambda * dlambda + els->Eupdiv8[kline] / T;

				   sumterm = els->gf013divw3[kline] * exp( expterm );

				   RespExtn = responsivity_2nd[kwave] * extinction[kwave];

				   if( kwave2nd >= 0  &&  kwave2nd < elemdata->nwave )  els->speclo[kwave2nd] += invsigma2pi_2nd * RespExtn * sumterm;

			  }

		 }

		 //-------- Scale the final answer

		 solid_angle = 3.141592654 * elemdata->plasma_radius_meters * elemdata->plasma_radius_meters / elemdata->range_meteor_meters / elemdata->range_meteor_meters;

		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  els->speclo[kwave] *= solid_angle * grating_area_scaling / pf;


		 //--------- Find the peak spectral value

		 kmax = 0;

		 maxspec = 0.0;
		 
		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  {
			  
			  if( maxspec < els->speclo[kwave])  {
			      maxspec = els->speclo[kwave];
				  kmax = kwave;
			  }
		 }

		 els->maxspeclo = maxspec;
         

	} //... end of element loop
				   

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the model spectrum at the Tlo temperature for the active elements to be used
//   in a fit. Needs to be called if sigma0, Thi or any user_fitflag changes from the last 
//   model spectrum calculation.
//-----------------------------------------------------------------------------------------

void  ActiveElementsSpectraTlo( struct  elements_data  *elemdata,
								double                 *responsivity_1st, 
								double                 *responsivity_2nd, 
							    double                  grating_area_scaling,
							    double                 *extinction )
{
int     kelem, kline, kwave, kmin, kwave2nd, kneu, kion;
int     kwave_start, kline_start, kline_finis, kpf, kmax;
double  invsigmasq_1st, invsigma2pi_1st;
double  invsigmasq_2nd, invsigma2pi_2nd;
double  sigma, sigmahalf, dlimit, dlambda, expterm, pf;
double  maxspec, sumterm, RespExtn, solid_angle;
struct  element_lines_spectrum  *els;


    //========== Compute broadening terms

    sigma = elemdata->sigma0 / 2.3548;
    // printf("Sigma... %f\n", elemdata->sigma0);

	sigmahalf = sigma / 2.0;
	// printf("Sigmahalf... %f\n", sigmahalf);

    invsigmasq_1st = -0.5 /  ( sigma * sigma );

	invsigma2pi_1st = +1.0 / ( sigma * sqrt( 2.0 * 3.141592654 ) );

    invsigmasq_2nd = -0.5 /  ( sigmahalf * sigmahalf );

	invsigma2pi_2nd = +1.0 / ( sigmahalf * sqrt( 2.0 * 3.141592654 ) );


	//========== Compute lambda cutoff range when contribution = 10^-6

	dlimit = sigma * sqrt( -2.0 * log( 1.0e-6 ) );
	// printf("dlimit... %f\n", dlimit);

	kmin = (int)( elemdata->wave[0] / ( elemdata->wave[1] - elemdata->wave[0] ) );


	//========== Set the compute model flag for low temperature

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  elemdata->els[kelem].compute_lo = 0;

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

		 if( elemdata->els[kelem].user_fitflag == FITLESS )  continue;    //... do not compute unused elements

		 if( elemdata->els[kelem].ioncode == 2 )  {  //... This is an ion so get the neutral index
			 kneu = elemdata->els[kelem].neuindex;   //       
			 kion = kelem;
		 }

		 else  {                          //... This is a neutral or molecule so get the ion index
			 kneu = kelem;
			 kion = elemdata->els[kelem].ionindex;
		 }

		 elemdata->els[kneu].compute_lo = 1;
		 elemdata->els[kion].compute_lo = 1;

	}


    //========== Loop through each element and build spectra

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

 		 els = &elemdata->els[kelem];

		 //-------- clear the spectrum

		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  els->speclo[kwave] = 1.0e-60;


		 //-------- Compute partition function at T using piecewise linear interpolation

	     kpf = 0;

	     while( elemdata->Tlo > els->Tpf[kpf+1]  &&  kpf < els->Npf-2 )  kpf++;

		 pf = els->partfunc[kpf]  +  ( els->partfunc[kpf+1] - els->partfunc[kpf] ) * ( elemdata->Tlo - els->Tpf[kpf] ) / ( els->Tpf[kpf+1] - els->Tpf[kpf] );

		 els->partfuncTlo = pf;


		 //-------- check if user has made this an active element (either neutral or ion)

		 if( els->compute_lo == 0 )  continue;


		 //printf("------------------------------------------  %s\n", els->element_string );

		 //-------- Loop over the emmission lines for this element
		 //         NOTE: line wavelengths are assumed sorted in ascending order
		 //               and this was checked when they were read in. This
		 //               allows the use of an increasing wavelength sliding 
		 //               window of output wavelengths that bracket the line 
		 //               wavelengths as we progress through the lines.

		 //-------- get the line start and stop indices

		 if( els->nlines > 0 )  {

		     kline_start = 0;

		     while( elemdata->wave[0] - els->wave_vacuum[kline_start] > dlimit  &&  kline_start < els->nlines-1 )  kline_start++;

		     kline_finis = els->nlines-1;

		     while( els->wave_vacuum[kline_finis] - elemdata->wave[elemdata->nwave-1] > dlimit  &&  kline_finis > 0 )  kline_finis--;


		     kwave_start = 0;

		     for( kline=kline_start; kline<=kline_finis; kline++ )  {

			      //- - - - - Shift the minimum output wavelength upwards, if it is > dlimit below the line wavelength

			      while( els->wave_vacuum[kline] - elemdata->wave[kwave_start] > dlimit  &&  kwave_start < elemdata->nwave-1 )  kwave_start++;


			      //- - - - - Sum this emmision line's contributing spectra over the narrow wavelength band defined by PSF (sigma)

			      for( kwave=kwave_start; kwave<elemdata->nwave; kwave++ )  {

				       dlambda = elemdata->wave[kwave] - els->wave_vacuum[kline];

				       if( dlambda > dlimit )  break;


				       //..... Compute this for 1st order

				       expterm = invsigmasq_1st * dlambda * dlambda + els->Eupdiv8[kline] / elemdata->Tlo;

				       sumterm = els->gf013divw3[kline] * exp( expterm );

				       RespExtn = responsivity_1st[kwave] * extinction[kwave];

				       els->speclo[kwave] += invsigma2pi_1st * RespExtn * sumterm;


				       //..... Compute this for 2nd order adding half the term because of 2:1 wavelength mapping 1st to 2nd order

				       kwave2nd = ( kwave - kmin ) / 2;

				       expterm = invsigmasq_2nd * dlambda * dlambda + els->Eupdiv8[kline] / elemdata->Tlo;

				       sumterm = els->gf013divw3[kline] * exp( expterm );

				       RespExtn = responsivity_2nd[kwave] * extinction[kwave];

				       if( kwave2nd >= 0  &&  kwave2nd < elemdata->nwave )  els->speclo[kwave2nd] += invsigma2pi_2nd * RespExtn * sumterm;

			      } //... end of narrow wavelength band loop

		     } //... end of emission line loop

		 } //... end if there any lines for the element


		 //-------- Scale the final answer 

		 solid_angle = 3.141592654 * elemdata->plasma_radius_meters * elemdata->plasma_radius_meters / elemdata->range_meteor_meters / elemdata->range_meteor_meters;

		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  els->spechi[kwave] *= solid_angle * grating_area_scaling / pf;

		 /*
		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  {  // looking for O777 to O844 lines and ratio
			 if( elemdata->wave[kwave] > 773.0 && elemdata->wave[kwave] < 781.0 )  printf("%lf  %lg\n", elemdata->wave[kwave], els->speclo[kwave] );
		     if( elemdata->wave[kwave] > 840.0 && elemdata->wave[kwave] < 848.0 )  printf("%lf  %lg\n", elemdata->wave[kwave], els->speclo[kwave] );
		 }
		 */

		 //-------- Find the peak spectral value

		 kmax = 0;

		 maxspec = 0.0;
		 
		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  {
			  
			  if( maxspec < els->speclo[kwave])  {
			      maxspec = els->speclo[kwave];
				  kmax = kwave;
			  }
		 }

		 els->maxspeclo = maxspec;

	} //... end of element loop
				   

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the model spectrum at the Thi temperature for the active elements to be used
//   in a fit. Needs to be called if sigma0, Thi or any user_fitflag changes from the last 
//   model spectrum calculation.
//-----------------------------------------------------------------------------------------

void  ActiveElementsSpectraThi( struct  elements_data  *elemdata,
								double                 *responsivity_1st, 
								double                 *responsivity_2nd, 
							    double                  grating_area_scaling,
							    double                 *extinction )
{
int     kelem, kline, kwave, kwave_start, kline_start, kline_finis, kpf;
int     kmin, kmax, kwave2nd, kneu, kion;
double  invsigmasq_1st, invsigma2pi_1st;
double  invsigmasq_2nd, invsigma2pi_2nd;
double  sigma, sigmahalf, dlimit, dlambda, expterm, pf;
double  maxspec, sumterm, RespExtn, solid_angle;
struct  element_lines_spectrum  *els;



    //========== Compute broadening terms

    sigma = elemdata->sigma0 / 2.3548;

	sigmahalf = sigma / 2.0;

    invsigmasq_1st = -0.5 /  ( sigma * sigma );

	invsigma2pi_1st = +1.0 / ( sigma * sqrt( 2.0 * 3.141592654 ) );

    invsigmasq_2nd = -0.5 /  ( sigmahalf * sigmahalf );

	invsigma2pi_2nd = +1.0 / ( sigmahalf * sqrt( 2.0 * 3.141592654 ) );


	//========== Compute lambda cutoff range when contribution = 10^-6

	dlimit = sigma * sqrt( -2.0 * log( 1.0e-6 ) );

	kmin = (int)( elemdata->wave[0] / ( elemdata->wave[1] - elemdata->wave[0] ) );


	//========== Set the compute model flag for high temperature

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  elemdata->els[kelem].compute_hi = 0;

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

		 if( elemdata->els[kelem].user_fitflag == FITLESS )  continue;  //... do not compute unused elements

		 if( elemdata->els[kelem].ioncode > 50 )  continue;   //... do not include molecular spectrum for Thi

		 if( elemdata->els[kelem].ioncode == 2 )  {  //... This is an ion so get the neutral index
			 kneu = elemdata->els[kelem].neuindex;   //       
			 kion = kelem;
		 }

		 else  {                          //... This is a neutral or molecule so get the ion index
			 kneu = kelem;
			 kion = elemdata->els[kelem].ionindex;
		 }

		 elemdata->els[kneu].compute_hi = 1;
		 elemdata->els[kion].compute_hi = 1;

	}


	//========== Loop through each element and build spectra

	for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

 		 els = &elemdata->els[kelem];

		 //........ clear the spectrum

		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  els->spechi[kwave] = 1.0e-60;


		 //-------- Compute partition function at T using piecewise linear interpolation

	     kpf = 0;

	     while( elemdata->Thi > els->Tpf[kpf+1]  &&  kpf < els->Npf-2 )  kpf++;

		 pf = els->partfunc[kpf]  +  ( els->partfunc[kpf+1] - els->partfunc[kpf] ) * ( elemdata->Thi - els->Tpf[kpf] ) / ( els->Tpf[kpf+1] - els->Tpf[kpf] );

		 els->partfuncThi = pf;


		 //-------- check if user made this an active element

		 if( els->compute_hi == 0 )  continue;


		 //printf("------------------------------------------  %s\n", els->element_string );


		 //-------- Loop over the emmission lines for this element
		 //         NOTE: line wavelengths are assumed sorted in ascending order
		 //               and this was checked when they were read in. This
		 //               allows the use of an increasing wavelength sliding 
		 //               window of output wavelengths that bracket the line 
		 //               wavelengths as we progress through the lines.

		 //-------- get the line start and stop indices

		 if( els->nlines > 0 )  {

		     kline_start = 0;

		     while( elemdata->wave[0] - els->wave_vacuum[kline_start] > dlimit  &&  kline_start < els->nlines-1 )  kline_start++;

		     kline_finis = els->nlines-1;

		     while( els->wave_vacuum[kline_finis] - elemdata->wave[elemdata->nwave-1] > dlimit  &&  kline_finis > 0 )  kline_finis--;


		     kwave_start = 0;

		     for( kline=kline_start; kline<=kline_finis; kline++ )  {

			      //- - - - - Shift the minimum output wavelength upwards, if it is > dlimit below the line wavelength

			      while( els->wave_vacuum[kline] - elemdata->wave[kwave_start] > dlimit  &&  kwave_start < elemdata->nwave-1 )  kwave_start++;


			      //- - - - - Sum the spectra over the narrow contributing wavelength band

			      for( kwave=kwave_start; kwave<elemdata->nwave; kwave++ )  {

				       dlambda = elemdata->wave[kwave] - els->wave_vacuum[kline];

				       if( dlambda > dlimit )  break;


				       //..... Compute this for 1st order

				       expterm = invsigmasq_1st * dlambda * dlambda + els->Eupdiv8[kline] / elemdata->Thi;

				       sumterm = els->gf013divw3[kline] * exp( expterm );

				       RespExtn = responsivity_1st[kwave] * extinction[kwave];

				       els->spechi[kwave] += invsigma2pi_1st * RespExtn * sumterm;


				       //..... Compute this for 2nd order adding half ??? the term because of 2:1 wavelength mapping 1st to 2nd order

				       kwave2nd = ( kwave - kmin ) / 2;

				       expterm = invsigmasq_2nd * dlambda * dlambda + els->Eupdiv8[kline] / elemdata->Thi;

				       sumterm = els->gf013divw3[kline] * exp( expterm );

				       RespExtn = responsivity_2nd[kwave2nd] * extinction[kwave2nd];

				       if( kwave2nd >= 0  &&  kwave2nd < elemdata->nwave )  els->spechi[kwave2nd] += invsigma2pi_2nd * RespExtn * sumterm;

			      } //... end of narrow wavelength band loop 

			 } //... end of emission line loop

		 } //... end if there any lines for the element


		 //-------- Scale the final answer (include the hot-to-warm ratio in the hot model spectra)

		 solid_angle = 3.141592654 * elemdata->plasma_radius_meters * elemdata->plasma_radius_meters / elemdata->range_meteor_meters / elemdata->range_meteor_meters;

		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  els->spechi[kwave] *= solid_angle * grating_area_scaling / pf;


		 //------- Find the peak spectral value

		 kmax = 0;

		 maxspec = 0.0;
		 
		 for( kwave=0; kwave<elemdata->nwave; kwave++ )  {
			  
			  if( maxspec < els->spechi[kwave])  {
			      maxspec = els->spechi[kwave];
				  kmax = kwave;
			  }
		 }

		 els->maxspechi = maxspec;
   

	} //... end of element loop
				   

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Performs an LMS fit of the active elements to the measured "spectrum".
//-----------------------------------------------------------------------------------------

void  FitSpectralCoefficients( double                     nwavelengths,
	                           double                    *wavelength_nm,
							   double                    *meas_spectrum, 
                               struct  elements_data     *elemdata,
							   struct  specconfiguration *spconfig )
{
int      nfit, kfit, jfit, kelem, jelem, kwave, k, kwlo, kwhi, allzero;
double   sum, spec_kelem, spec_jelem; 
double   maxSpectrum, maxspec_kelem, maxspec_jelem;
double   Qts1, QtQ1;

double   *Qts;
double   *Ssq;  // squared singular values
double   *work;
double  **QtQ;  // svd of QtQ = USV'
double  **US;
double  **V;


   //========== Determine number of active elements to be used in the fit 
   //              by checking for the FITTING flag but not fit locked (FITLOCK)
   printf("Determining the number of active elements...\n");
   printf("%d\n", elemdata->nelements);

   kfit = 0;

   for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

   		printf("User fit flag: %d\n", elemdata->els[kelem].user_fitflag);
   		// printf("%d\n", kelem);
   		// printf("Init fit flag: %d\n", elemdata->els[kelem].init_fitflag);

	    if( elemdata->els[kelem].user_fitflag == FITTING )  {

			elemdata->fit_element_index[kfit] = kelem;
			
			kfit++;
			printf("kfit: %d\n", kfit);

		}
   }

   nfit = kfit;

   elemdata->nfit = nfit;

   printf("Number of active elements = %d\n", nfit);


   //========== Find limits and largest value in the spectrum and check wavelengths match
   
   maxSpectrum = 0.0;

   allzero = 1;  // assume all zeros

   for( kwave=0; kwave<nwavelengths; kwave++ )  {

	    if( wavelength_nm[kwave] <= spconfig->min_fit_wavelength_nm )  kwlo = kwave;
	    if( wavelength_nm[kwave] <= spconfig->max_fit_wavelength_nm )  kwhi = kwave;

	    if( wavelength_nm[kwave] < spconfig->min_fit_wavelength_nm  ||
			wavelength_nm[kwave] > spconfig->max_fit_wavelength_nm     )  continue;

		if(meas_spectrum[kwave] != 0.0 )  allzero = 0;  // non-zero spectral value found

		if( maxSpectrum < meas_spectrum[kwave] )
			maxSpectrum = meas_spectrum[kwave];

		if( fabs( wavelength_nm[kwave] - elemdata->wave[kwave] ) > 1.0e-6 )  {
			printf("FitSpectralCoefficients called with wavelength range different from ReadElementsData \n");
			Delay_msec(15000);
			exit(20);
		}
   }

  

   //========== Don't fit spectrum if all the values are zero

   if( allzero == 1  &&  nfit > 0 )  {

	    for( kfit=0; kfit<nfit; kfit++ )  {

			kelem = elemdata->fit_element_index[kfit];

			elemdata->els[kelem].N_warm         = 0.0;
			elemdata->els[kelem].N_warm_default = 0.0;
		}
   }


   //========== Fit spectrum for no elements

   else if( nfit == 0 )  {

       for( kelem=0; kelem<elemdata->nelements; kelem++ )  {
		    elemdata->els[kelem].N_warm         = 0.0;
		    elemdata->els[kelem].N_warm_default = 0.0;
	   }

   }


   //========== Fit spectrum for one element

   else if( nfit == 1 )  {

       //-------- Compute Q-transpose * spectrum and Q-transpose * Q
	    
       Qts1 = 0.0;
       QtQ1 = 0.0;

	   kelem = elemdata->fit_element_index[0];

       for( kwave=kwlo; kwave<=kwhi; kwave++ )  {

			if( elemdata->els[kelem].warmhot == WARM )  spec_kelem = elemdata->els[kelem].speclo[kwave];
			if( elemdata->els[kelem].warmhot == HOT  )  spec_kelem = elemdata->els[kelem].spechi[kwave];

		    Qts1 += spec_kelem * meas_spectrum[kwave];

		    QtQ1 += spec_kelem * spec_kelem;

       }

       elemdata->els[kelem].N_warm = Qts1 / QtQ1;

	   if( elemdata->els[kelem].N_warm < 0.0 )   elemdata->els[kelem].N_warm = 0.0;

	   elemdata->els[kelem].N_warm_default = elemdata->els[kelem].N_warm;

   }  //... end of fit spectrum for one element


   //========== Fit spectrum for more than one element

   else if( nfit > 1 )  {

       //-------- Allocate memory for 1D and 2D arrays where Q
       //         respresents the user desired element spectra to fit.
       //                  Qts[nfit]
       //                  Ssq[nfit];
       //                  work[nfit];
       //                  QtQ[nfit][nfit]
       //                  US[nfit][nfit]
       //                  V[nfit][nfit]

       Qts  = (double*) malloc( nfit * sizeof(double) );
       Ssq  = (double*) malloc( nfit * sizeof(double) );
       work = (double*) malloc( nfit * sizeof(double) );

       QtQ = (double**) malloc( (nfit * sizeof(double*)) + (nfit * nfit * sizeof(double)) );
       for( k=0; k<nfit; ++k )  QtQ[k] = (double*)(QtQ + nfit) + k*nfit;

       US = (double**) malloc( (nfit * sizeof(double*)) + (nfit * nfit * sizeof(double)) );
       for( k=0; k<nfit; ++k )  US[k] = (double*)(US + nfit) + k*nfit;

       V = (double**) malloc( (nfit * sizeof(double*)) + (nfit * nfit * sizeof(double)) );
       for( k=0; k<nfit; ++k )  V[k] = (double*)(V + nfit) + k*nfit;


       //-------- Fill the Q-transpose * spectrum vector

       for( kfit=0; kfit<nfit; kfit++ )  {

	        kelem = elemdata->fit_element_index[kfit];

	        sum = 0.0;

		    for( kwave=kwlo; kwave<=kwhi; kwave++ )  {

			     if( elemdata->els[kelem].warmhot == WARM )  spec_kelem = elemdata->els[kelem].speclo[kwave];
			     if( elemdata->els[kelem].warmhot == HOT  )  spec_kelem = elemdata->els[kelem].spechi[kwave];

			     sum += spec_kelem * meas_spectrum[kwave];

		    }

		    maxspec_kelem = elemdata->els[kelem].maxspeclo + elemdata->els[kelem].maxspechi;

		    Qts[kfit] = sum / maxspec_kelem / maxSpectrum;

	        ////printf( "Fitting %d Qts=%lf\n", kfit, Qts[kfit] );

       }


       //-------- Fill the QtQ matrix

       for( jfit=0; jfit<nfit; jfit++ )  {

	        ////printf( "Fitting  %d  QtQ=", kfit );
	    
		    jelem = elemdata->fit_element_index[jfit];

            for( kfit=jfit; kfit<nfit; kfit++ )  {

	             kelem = elemdata->fit_element_index[kfit];

	             sum = 0.0;

		         for( kwave=kwlo; kwave<=kwhi; kwave++ )  {

			          if( elemdata->els[kelem].warmhot == WARM )  spec_kelem = elemdata->els[kelem].speclo[kwave];
			          if( elemdata->els[kelem].warmhot == HOT  )  spec_kelem = elemdata->els[kelem].spechi[kwave];

			          if( elemdata->els[kelem].warmhot == WARM )  spec_jelem = elemdata->els[jelem].speclo[kwave];
			          if( elemdata->els[kelem].warmhot == HOT  )  spec_jelem = elemdata->els[jelem].spechi[kwave];

			          sum += spec_kelem * spec_jelem;

			     }

		         maxspec_kelem = elemdata->els[kelem].maxspeclo + elemdata->els[kelem].maxspechi;

		         maxspec_jelem = elemdata->els[jelem].maxspeclo + elemdata->els[jelem].maxspechi;

			     sum = sum / maxspec_kelem / maxspec_jelem;

			     QtQ[jfit][kfit] = sum;
			     QtQ[kfit][jfit] = sum;

	             ////printf( "  %lf", QtQ[jfit][kfit] );
		    }

		    ////printf("\n");

       }


       //-------- Compute SVD of Q-transpose * Q = USV'

       svd_NxN( *QtQ, *US, Ssq, *V, nfit, nfit );


       //-------- Obtain solution flux = V * US' * (Qt * spectrum)

       for( kfit=0; kfit<nfit; kfit++ )  {
	   
	        sum = 0.0;

		    for( jfit=0; jfit<nfit; jfit++ )  sum += US[jfit][kfit] * Qts[jfit];
		
		    if( Ssq[kfit] > 0.0 )  work[kfit] = sum  / Ssq[kfit];
		    else                   work[kfit] = 0.0;
       }


       for( kfit=0; kfit<nfit; kfit++ )  {

		    kelem = elemdata->fit_element_index[kfit];
		 
	        sum = 0.0;

		    for( jfit=0; jfit<nfit; jfit++ )  sum += V[kfit][jfit] * work[jfit];

		    maxspec_kelem = elemdata->els[kelem].maxspeclo + elemdata->els[kelem].maxspechi;

		    elemdata->els[kelem].N_warm = sum * maxSpectrum / maxspec_kelem;

		    ////printf( "Fitting  %d  flux=%lf\n", kfit, elemdata->els[kelem].flux );

			if( elemdata->els[kelem].N_warm < 0.0 )   elemdata->els[kelem].N_warm = 0.0;

	        elemdata->els[kelem].N_warm_default = elemdata->els[kelem].N_warm;

       }
                     

       //-------- Free memory

       free( Qts  );
       free( Ssq  );
       free( work );
       free( QtQ  );
       free( US   );
       free( V    );

   }  //... end of fit more than one element
  
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the spectrum for the "locked" elements and associated coefficients
//-----------------------------------------------------------------------------------------

void  SpectrumGivenOnlyLockedCoefs( struct  elements_data     *elemdata,
							        double                    *fitspectrum )
{
int      kneutral, kelem, kwave, kneu, kion;
double   sum; 


   //========== loop over all output wavelengths

   for( kwave=0; kwave<elemdata->nwave; kwave++ )  {

	    sum = 1.0e-99;

        //-------- Loop over all NEUTRAL elements to get column densities

	    for( kneutral=0; kneutral<elemdata->nneutrals; kneutral++ )  {

		     kelem = elemdata->neutral_index[kneutral];
 

			 //........ Accumulate only if element was already fit and locked

	         if( elemdata->els[kelem].user_fitflag == FITLOCK )  {	  

				 kneu = kelem;
			     kion = elemdata->els[kelem].ionindex;

			     sum += elemdata->els[kneu].N_warm * elemdata->els[kneu].speclo[kwave]
					 +  elemdata->els[kion].N_warm * elemdata->els[kion].speclo[kwave]
			         +  elemdata->els[kneu].N_hot  * elemdata->els[kneu].spechi[kwave]
					 +  elemdata->els[kion].N_hot  * elemdata->els[kion].spechi[kwave];

			 }

		}

		fitspectrum[kwave] = sum;

	}

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the spectrum for the "fitting" elements and associated coefficients
//-----------------------------------------------------------------------------------------

void  SpectrumGivenOnlyFittingCoefs( struct  elements_data     *elemdata,
						        	 double                    *fitspectrum )
{
int      kneutral, kelem, kwave, kneu, kion;
double   sum; 


   //========== loop over all output wavelengths

   for( kwave=0; kwave<elemdata->nwave; kwave++ )  {

	    sum = 1.0e-99;

        //-------- Loop over all NEUTRAL elements to get column densities

	    for( kneutral=0; kneutral<elemdata->nneutrals; kneutral++ )  {

		     kelem = elemdata->neutral_index[kneutral];
 

			 //........ Accumulate only if element was fit but not locked

	         if( elemdata->els[kelem].user_fitflag == FITTING )  {	  

				 kneu = kelem;
			     kion = elemdata->els[kelem].ionindex;

			     sum += elemdata->els[kneu].N_warm * elemdata->els[kneu].speclo[kwave]
					 +  elemdata->els[kion].N_warm * elemdata->els[kion].speclo[kwave]
			         +  elemdata->els[kneu].N_hot  * elemdata->els[kneu].spechi[kwave]
					 +  elemdata->els[kion].N_hot  * elemdata->els[kion].spechi[kwave];

			 }

		}

		fitspectrum[kwave] = sum;

	}

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the spectrum for ALL the coefficients (fitting and locked)
//-----------------------------------------------------------------------------------------

void  SpectrumGivenAllCoefs( struct  elements_data     *elemdata,
							 double                    *fitspectrum )
{
int      kneutral, kelem, kwave, kneu, kion;
double   sum; 
	printf("elemdata\n");
	printf("fitspectrum\n");

   //========== loop over all output wavelengths

   for( kwave=0; kwave<elemdata->nwave; kwave++ )  {

	    sum = 1.0e-99;

        //-------- Loop over all NEUTRAL elements to get column densities

	    for( kneutral=0; kneutral<elemdata->nneutrals; kneutral++ )  {

		     kelem = elemdata->neutral_index[kneutral];
 

			 //....... Accumulate only if element was actively fit. Note that all the column
			 //        densities needed for the calculation are stored in the neutral "kneu"
			 //        or "kion" index of the element_lines_spectrum structure. 

	         if( elemdata->els[kelem].user_fitflag != FITLESS )  {	  

				 kneu = kelem;
			     kion = elemdata->els[kelem].ionindex;

			     sum += elemdata->els[kneu].N_warm * elemdata->els[kneu].speclo[kwave]
					 +  elemdata->els[kion].N_warm * elemdata->els[kion].speclo[kwave]
			         +  elemdata->els[kneu].N_hot  * elemdata->els[kneu].spechi[kwave]
					 +  elemdata->els[kion].N_hot  * elemdata->els[kion].spechi[kwave];

			 }

		}

		fitspectrum[kwave] = sum;

	}

	printf("Successfully fit. I think...");


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Compute the spectrum for ALL the coefficients (fit and locked) for specific
//   settings of neutral or ion, and warm or hot.
//   That is a spectrum is computed for one the following four combinations:
//         Neutral - Warm
//         Neutral - Hot
//         Ion     - Warm
//         Ion     - Hot
//-----------------------------------------------------------------------------------------

void  SpectrumGivenAllCoefs_Subset( struct  elements_data     *elemdata,
							        double                    *fitspectrum,
									int                        neutral_ion_display,
									int                        warm_hot_display )
{
int      kneutral, kelem, kwave, kneu, kion;
double   sum; 


   //========== neutral_ion_display:  0 = show all  
   //                                 1 = show only neutrals and molecules
   //                                 2 = show only ions


   //========== loop over all output wavelengths

   for( kwave=0; kwave<elemdata->nwave; kwave++ )  {

	    sum = 1.0e-99;

        //-------- Loop over all NEUTRAL elements to get column densities

	    for( kneutral=0; kneutral<elemdata->nneutrals; kneutral++ )  {

		     kelem = elemdata->neutral_index[kneutral];
 
			 
			 //........ Accumulate only if element was actively fit

	         if( elemdata->els[kelem].user_fitflag != FITLESS )  {	 //... accumulate if element was fit   

				 kneu = kelem; 
			     kion = elemdata->els[kelem].ionindex;

			     //.... Add in the NEUTRAL WARM contribution for this element

			     if( neutral_ion_display != 2  &&  warm_hot_display != 2 )  {

					 sum += elemdata->els[kneu].N_warm * elemdata->els[kneu].speclo[kwave];
				 }

				 //.... Add in the ION WARM contribution for this element

			     if( neutral_ion_display != 1  &&  warm_hot_display != 2 )  {

				     sum += elemdata->els[kion].N_warm * elemdata->els[kion].speclo[kwave];
				 }

			     //.... Add in the NEUTRAL HOT contribution for this element

			     if( neutral_ion_display != 2  &&  warm_hot_display != 1 )  {

			         sum += elemdata->els[kneu].N_hot  * elemdata->els[kneu].spechi[kwave];
				 }

			     //.... Add in the ION HOT contribution for this element

			     if( neutral_ion_display != 1  &&  warm_hot_display != 1 )  {

				     sum += elemdata->els[kion].N_hot  * elemdata->els[kion].spechi[kwave];
				 }

			 }

		}

		fitspectrum[kwave] = sum;

	}


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    Performs an LMS fit for Fe only to the measured "spectrum".
//-----------------------------------------------------------------------------------------

void    FitFeSpectrum( double                     nwavelengths,
	                   double                    *wavelength_nm,
					   double                    *spectrum, 
                       struct  elements_data     *elemdata,
					   struct  specconfiguration *spconfig,
					   double                    *Fe_coef,
					   double                    *Fe_residual )
{
int      kelem, kwave, kwlo = 0, kwhi = 0;
double   sum, spec_kelem, maxspec_kelem;
double   Qts, QtQ;



   //========== Determine the element index for Iron

   for( kelem=0; kelem<elemdata->nelements; kelem++ )  {

	    if( elemdata->els[kelem].elementcode == 26.0 )  break;

   }

   if( kelem >= elemdata->nelements )  {
	   printf("FitFeSpectrumCoefficient could not identify element 26 (Fe) in elemdata structure \n");
	   Delay_msec(15000);
	   exit(20);
   }


   //========== Check if wavelengths match
   
   for( kwave=0; kwave<nwavelengths; kwave++ )  {

	    if( wavelength_nm[kwave] <= spconfig->min_fit_wavelength_nm )  kwlo = kwave;
	    if( wavelength_nm[kwave] <= spconfig->max_fit_wavelength_nm )  kwhi = kwave;

	    if( wavelength_nm[kwave] < spconfig->min_fit_wavelength_nm  ||
			wavelength_nm[kwave] > spconfig->max_fit_wavelength_nm     )  continue;

	    if( fabs( wavelength_nm[kwave] - elemdata->wave[kwave] ) > 1.0e-6 )  {
			printf("FitFeSpectrumCoefficient called with wavelength range different from ReadElementsData \n");
			Delay_msec(15000);
			exit(20);
		}
   }


   //========== Compute Q-transpose * spectrum and Q-transpose * Q
	    
   Qts = 0.0;
   QtQ = 0.0;

   maxspec_kelem = 0.0;

   for( kwave=kwlo; kwave<=kwhi; kwave++ )  {

		spec_kelem = elemdata->els[kelem].speclo[kwave] + elemdata->els[kelem].spechi[kwave];

		if( maxspec_kelem < spec_kelem )  maxspec_kelem = spec_kelem;

		Qts += spec_kelem * spectrum[kwave];

		QtQ += spec_kelem * spec_kelem;

   }


   //========== Obtain solution flux = inverse(Qt * Q) * (Qt * spectrum)

   *Fe_coef = Qts / QtQ;


   //------------- Compute fit Fe spectrum to actual measured spectrum residual

   sum = 0.0;

   for( kwave=kwlo; kwave<=kwhi; kwave++ )  {

		spec_kelem = elemdata->els[kelem].speclo[kwave] + elemdata->els[kelem].spechi[kwave];

		if( spec_kelem > 1.0e-6 * maxspec_kelem )  sum += ( *Fe_coef * spec_kelem - spectrum[kwave] )
			                                            * ( *Fe_coef * spec_kelem - spectrum[kwave] );

   }

   *Fe_residual = sqrt( sum );


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
 * Perform a singular value decomposition on A = USV' of a nxn square matrix.
 *
 * Call this with  svd_NxN( &a[0][0], &us[0][0], &s2[0], &v[0][0], n, maxn );
 *
 * This routine was originally adapted from a Pascal implementation
 * (c) 1988 J. C. Nash, "Compact numerical methods for computers", Hilger 1990.
 * The A, US, V matrices must be pre-allocated with minimum n rows and n columns. 
 * Upon calling, the matrix to be decomposed is contained in the first n rows 
 * of A. On return, the US matrix contains the product U * sqrt(S2) and the V
 * matrix contains V (not V'). The S2 vector returns the square of the singular
 * values. Note that maxn is the max column dimension of A, US, and V.
*/
//-------------------------------------------------------------------------------

void svd_NxN( double *A, double *US, double *S2, double *V, int n, int maxn )
{
  int     i, j, k, EstColRank = n, RotCount = n; 
  int     SweepCount = 0, slimit = (n<120) ? 30 : n/4;
  double  eps = 1e-15, e2 = 10.0*n*eps*eps, tol = 0.1*eps;
  double  vt, p, x0, y0, q, r, c0, s0, d1, d2;

  for (i=0; i<n; i++) { for (j=0; j<n; j++) *(US+i*maxn+j) = *(A+i*maxn+j); }

  for (i=0; i<n; i++) { for (j=0; j<n; j++) *(V+i*maxn+j) = 0.0; *(V+i*maxn+i) = 1.0; }
  while (RotCount != 0 && SweepCount++ <= slimit) {
    RotCount = EstColRank*(EstColRank-1)/2;
    for (j=0; j<EstColRank-1; j++) 
      for (k=j+1; k<EstColRank; k++) {
        p = q = r = 0.0;
        for (i=0; i<n; i++) {
          x0 = *(US+i*maxn+j); y0 = *(US+i*maxn+k);
          p += x0*y0; q += x0*x0; r += y0*y0;
        }
        S2[j] = q; S2[k] = r;
        if (q >= r) {
          if (q<=e2*S2[0] || fabs(p)<=tol*q)
            RotCount--;
          else {
            p /= q; r = 1.0-r/q; vt = sqrt(4.0*p*p+r*r);
            c0 = sqrt(0.5*(1.0+r/vt)); s0 = p/(vt*c0);
            for (i=0; i<n; i++) {
              d1 = *(US+i*maxn+j); d2 = *(US+i*maxn+k);
              *(US+i*maxn+j) = d1*c0+d2*s0; *(US+i*maxn+k) = -d1*s0+d2*c0;
            }
            for (i=0; i<n; i++) {
              d1 = *(V+i*maxn+j); d2 = *(V+i*maxn+k);
              *(V+i*maxn+j) = d1*c0+d2*s0; *(V+i*maxn+k) = -d1*s0+d2*c0;
            }
          }
        } else {
          p /= r; q = q/r-1.0; vt = sqrt(4.0*p*p+q*q);
          s0 = sqrt(0.5*(1.0-q/vt));
          if (p<0.0) s0 = -s0;
          c0 = p/(vt*s0);
          for (i=0; i<n; i++) {
            d1 = *(US+i*maxn+j); d2 = *(US+i*maxn+k);
            *(US+i*maxn+j) = d1*c0+d2*s0; *(US+i*maxn+k) = -d1*s0+d2*c0;
          }
          for (i=0; i<n; i++) {
            d1 = *(V+i*maxn+j); d2 = *(V+i*maxn+k);
            *(V+i*maxn+j) = d1*c0+d2*s0; *(V+i*maxn+k) = -d1*s0+d2*c0;
          }
        }
      }
    while (EstColRank>2 && S2[EstColRank-1]<=S2[0]*tol+tol*tol) EstColRank--;
  }
  if (SweepCount > slimit)
    printf("Warning: Reached maximum number of sweeps (%d) in SVD routine...\n"
	   ,slimit);
}
//                                                             svd_NxN
//====================================================================


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Converts year, month, day, ... into a Julian date and time
//-----------------------------------------------------------------------------------------

double  JulianDateAndTime(long year, long month, long day,
	                      long hour, long minute, long second, long milliseconds)
{
	long    aterm, bterm, jdate;
	double  jdt;


	//========== compute julian date terms

	if ((month == 1L) || (month == 2L)) {
		year = year - 1L;
		month = month + 12L;
	}
	aterm = (long)(year / 100L);
	bterm = 2L - aterm + (long)(aterm / 4L);

	//========== julian day at 12hr UT

	jdate = (long)(365.25 * (double)(year + 4716L))
		+ (long)(30.6001 * (double)(month + 1L))
		+ day + bterm - 1524L;

	//========== add on the actual UT

	jdt = (double)jdate - 0.5
		+ (double)hour / (24.0)
		+ (double)minute / (24.0 * 60.0)
		+ (double)second / (24.0 * 60.0 * 60.0)
		+ (double)milliseconds / (24.0 * 60.0 * 60.0 * 1000.0);

	return(jdt);

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Converts Julian date and time into year, month, ....
//-----------------------------------------------------------------------------------------

void  CalendarDateAndTime( double jdt,
	                       long *year, long *month, long *day,
	                       long *hour, long *minute, long *second, long *milliseconds)
{
	long        jyear, jmonth, jday, jhour, jminute, jsecond, jmilliseconds;
	long        aterm, bterm, cterm, dterm, eterm, alpha, intpart;
	double      fracpart;

	//========== JD integer and fractional part

	intpart = (long)(jdt + 0.5);
	fracpart = (jdt + 0.5) - (double)intpart;

	alpha = (long)(((double)intpart - 1867216.25) / 36524.25);
	aterm = intpart + 1 + alpha - (long)(alpha / 4L);
	bterm = aterm + 1524L;
	cterm = (long)(((double)bterm - 122.1) / 365.25);
	dterm = (long)(365.25 * (double)cterm);
	eterm = (long)(((double)bterm - (double)dterm) / 30.6001);

	//========== date calculation

	jday = bterm - dterm - (long)(30.6001 * (double)eterm);

	if (eterm < 14L)  jmonth = (long)eterm - 1;
	else              jmonth = (long)eterm - 13;

	if (jmonth > 2L)  jyear = (long)cterm - 4716;
	else              jyear = (long)cterm - 4715;


	//========== time calculation

	fracpart += 0.0001 / (24.0 * 60.0 * 60.0);  // Add 0.1 msec for rounding
	fracpart *= 24.0;
	jhour = (long)fracpart;
	fracpart -= (double)jhour;
	fracpart *= 60.0;
	jminute = (long)fracpart;
	fracpart -= (double)jminute;
	fracpart *= 60.0;
	jsecond = (long)fracpart;
	fracpart -= (double)jsecond;
	fracpart *= 1000.0;
	jmilliseconds = (long)fracpart;


	//========== put into output pointers

	*year         = jyear;
	*month        = jmonth;
	*day          = jday;
	*hour         = jhour;
	*minute       = jminute;
	*second       = jsecond;
	*milliseconds = jmilliseconds;


}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   Writes some kind of spectrum data to a file
//   FIXME: someone please describe how this is actually used by the library
//-----------------------------------------------------------------------------------------
void  WriteSpectrum(const char *pathname, int nwavelengths, double *wavelength, double *spectrum1, double *spectrum2)
{
	FILE *outfile;

	if ((outfile = fopen(pathname, "wt")) == NULL) {
		printf(" Cannot open output spectra file %s for writing\n", pathname);
		return;
	}

	fprintf(outfile, "%d\n", nwavelengths);

	for (int kwave = 0; kwave < nwavelengths; kwave++) {

		if( spectrum2 == NULL ) fprintf(outfile, "%8.3lf  %15.6le\n", wavelength[kwave], spectrum1[kwave] );
		else                    fprintf(outfile, "%8.3lf  %15.6le  %15.6le\n", wavelength[kwave], spectrum1[kwave], spectrum2[kwave]);

	}

	fclose(outfile);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

