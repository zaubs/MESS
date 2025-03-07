

//====================================================================
//        Structures
//====================================================================

struct  StarInfo
{
	int      hip;             // Hipparcos star identifier
	double   ra_deg;          // RA in degrees (J2000)
	double   dec_deg;         // DEC in degrees (J2000)
	double   vmag;            // V magnitude 
	double   bv;              // B-V magnitude
	double   vi;              // V-I magnitude
	double   bv0;             // B-V (0) magnitude
	char     name[32];        // Star name
	double  *specval;         // Vector of spectral values (units?)
	double   specmax;         // Max spectral value
	double   specmin;         // Min spectral value > 0
};

struct  StarSpectra
{
	int      nwave;           // Number of wavelengths in the data file
	int      nstars;          // Number of stars in the data file
	double  *wave;            // Vector of spectral wavelengths (nm)
	struct   StarInfo  *star; // Structure to the star information
};



//====================================================================
//        Prototypes
//====================================================================

void  ReadStarSpectra( char *pathname, struct StarSpectra *starspect );

void  InterpolateStarSpectrum( struct StarSpectra *starspect, int kstar, int nwave, double *wave_spec, double *star_spec );

void  FreeMemoryStarSpectra( struct StarSpectra *starspect );



//====================================================================
//        Functions
//====================================================================

void  ReadStarSpectra( char *pathname, struct StarSpectra *starspect )
{
int     kw, ks, nstars, nwave_star;
double  star_specscale;
char    textline[128];
char    star_filename[256];
FILE   *starfile;

    
    //------ build the filenames for stars and planets

    strcpy( star_filename, pathname );

	strcat( star_filename, "StarSpectra.txt" );
      

    //------ open the stars and planets spectra files for reading

	starfile = fopen( star_filename, "rt" );

	if( starfile == NULL )  {
		printf(" Cannot open star spectra file %s for reading\n", star_filename );
		exit(30);
	}


	//------ Read the header lines of both files

	fscanf( starfile, "%[^\n]\n", textline );

	fscanf( starfile, "%d%d%lf", &nwave_star, &nstars, &star_specscale );


	starspect->nstars = nstars;

	starspect->nwave = nwave_star;


    //------ Allocate memory for wavelengths, stars, and their spectra

	starspect->wave = (double*) malloc( starspect->nwave * sizeof(double) );

	starspect->star = (struct StarInfo*) malloc( starspect->nstars * sizeof( struct StarInfo ) );

	for( ks=0; ks<starspect->nstars; ks++ )  {

		starspect->star[ks].specval = (double*) malloc( starspect->nwave * sizeof(double) );

	}


	//------ Read the star file's Hipparcos #, RA, Dec, Vmag, B-V, V-I, BV0, and star name 

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%d", &starspect->star[ks].hip );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%lf", &starspect->star[ks].ra_deg );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%lf", &starspect->star[ks].dec_deg );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%lf", &starspect->star[ks].vmag );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%lf", &starspect->star[ks].bv );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%lf", &starspect->star[ks].vi );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%lf", &starspect->star[ks].bv0 );

	fscanf( starfile, "\n" );

	for( ks=0; ks<nstars; ks++ )  fscanf( starfile, "%s", &starspect->star[ks].name );

	fscanf( starfile, "\n" );


	//------ Clear the min and max spectral value per star and planet

	for( ks=0; ks<starspect->nstars; ks++ )  starspect->star[ks].specmax = 0.0;
	for( ks=0; ks<starspect->nstars; ks++ )  starspect->star[ks].specmin = 1.0e+30;


	//------ Read the star wavelengths and spectra

	for( kw=0; kw<starspect->nwave; kw++ )  {

		fscanf( starfile, "%lf", &starspect->wave[kw] );

		for( ks=0; ks<starspect->nstars; ks++ )  {

			fscanf( starfile, "%lf", &starspect->star[ks].specval[kw] );
			
			starspect->star[ks].specval[kw] *= star_specscale;

			if( starspect->star[ks].specmax < starspect->star[ks].specval[kw] )  {
				starspect->star[ks].specmax = starspect->star[ks].specval[kw];
			}

			if( starspect->star[ks].specval[kw] > 0.0         &&
			    starspect->star[ks].specmin > starspect->star[ks].specval[kw] )  {
				starspect->star[ks].specmin = starspect->star[ks].specval[kw];
			}

		}

		fscanf( starfile, "\n" );
	}


	//------ Close the star and planet spectral files

	fclose(starfile);


}


//---------------------------------------------------------------------------------------

void  InterpolateStarSpectrum( struct StarSpectra *starspect, int kstar, int nwave, double *wave_spec, double *star_spec )
{
int     kwave, kcurr;
double  mslope;


    kcurr = 0;

    for( kwave=0; kwave<nwave; kwave++ )  {

		//-------- Are we below or above the star spectral wavelength band

		if( wave_spec[kwave] < starspect->wave[0]  ||
			wave_spec[kwave] > starspect->wave[starspect->nwave-1] )  {

			star_spec[kwave] = starspect->star[kstar].specmin;

			continue;
		}


		//-------- Find the wavelength bin in the star spectrum vector 
		//           (assumed monotonically increasing in wavelength)

        while( kcurr < starspect->nwave-2 )  {

			if( wave_spec[kwave] <= starspect->wave[kcurr+1] ) break;

			kcurr++;
		}


		//-------- Compute spectral slope in the wavelength bin and do linear interpolation

		mslope = ( starspect->star[kstar].specval[kcurr+1] - starspect->star[kstar].specval[kcurr] )
			   / ( starspect->wave[kcurr+1]                - starspect->wave[kcurr]    );

		star_spec[kwave] = starspect->star[kstar].specval[kcurr]  +  mslope * ( wave_spec[kwave] - starspect->wave[kcurr] );


	} //... End of interpolated star spectrum wavelength loop


}


//--------------------------------------------------------------------------------------

void  FreeMemoryStarSpectra( struct StarSpectra *starspect )
{
int   ks;

	for( ks=0; ks<starspect->nstars; ks++ )  free( starspect->star[ks].specval );

	free( starspect->star );

	free( starspect->wave );

}


//---------------------------------------------------------------------------------------

