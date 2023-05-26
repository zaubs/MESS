
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Spectral functions  -  Header file
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Input functions for reading configs, element emission lines, partition functions.
// IO and processing functions for extinction and responsivity.
// Processing functions for computing an emission line spectrum. 
//
//    Date     Ver   Author        Description
// ----------  ----  ------------  ---------------------------------------------------------
// 2021-05-12  1.00  Gural	       Derived from CAMS-Spectral code base
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef _H_GUARD_SPECTRAL_FUNCTIONS
#define _H_GUARD_SPECTRAL_FUNCTIONS


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%    Mnemonics and Hard Dimensioning     %%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define  MAXGRATINGS   32    // max number of cameras with gratings
#define  MAXTPF        20    // max # of temperatures in partition functions file

#define  WARM           0    // warmhot flag
#define  HOT            1
#define  BOTH           2

#define  FITLESS        0    // user_fitflag
#define  FITTING        1
#define  FITLOCK        2

#define  RAYLEIGH       1    // Extinction components
#define  O2BAND1        2
#define  O2BAND2        3
#define  OZONE          4
#define  DUST           5
#define  WATER          6

#define  RECLEAR        0   // REflag = Responsivity and extinction control flag
#define  REDAILY        1
#define  REMULTI        2

#define  UNITYWEIGHTS   1   // Weighting schema for responsivity and extinction averaging
#define  VMAGWEIGHTS    2


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%       Structure definitions       %%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//===============  Spectral processing configuration file input parameters

struct specconfiguration     
{
	double   version;                 // Config file version number
	long     order4spcalib;             // Order to use for spectral calibration
	long     rowdelta;                // Plus/minus #rows for integrating spectra in pixels
	double   min_cal_wavelength_nm;   // Minimum wavelength to calibrate response in nm
	double   max_cal_wavelength_nm;   // Maximum wavelength to calibrate response in nm
	double   del_wavelength_nm;       // Stepping wavelength to calibrate response in nm
	double   min_fit_wavelength_nm;   // Minimum fitting wavelength in nm
	double   max_fit_wavelength_nm;   // Maximum fitting wavelength in nm
    double   minspan_wavelength_nm;   // Min bandwidth the order4spcalib spectrum must span to use
	double   smooth_bandwidth_nm;     // 1 sigma width of Gaussian smoothing filter
	double   faintest_star_vmag;      // Faintest star to process in mV
    double   airmass_limit;           // Airmass value for inclusion in response or extinction
	double   fading_coef;             // Fading memory coefficient
	double   coin_time_tolerance;     // Time tolerance for coincidence in seconds
	double   min_lo_exc_temp;         // Minimum low excitation temperature in deg K
	double   max_lo_exc_temp;         // Maximum low excitation temperature in deg K
	double   step_lo_exc_temp;        // Step low excitation temperature in deg K
	double   nominal_lo_exc_temp;     // Nominal low excitation temperature in deg K
	double   nominal_hi_exc_temp;     // Nominal high excitation temperature in deg K
	double   nominal_sigma0;          // Nominal broadening in nm
	double   grating_offnormal_deg;   // Grating normal to look direction angle in deg
	double   default_roll_deg;        // Default grating roll angle in deg
	double   default_pitch_deg;       // Default grating picth angle in deg
	double   default_yaw_deg;         // Default grating yaw angle in deg
	double   default_ne;              // Default density of electrons
	double   default_hot2warm;        // Default hot to warm plasma contribution ratio
	long     ncams_grating;           // Number of cameras with gratings
	long     camnum[MAXGRATINGS];     // Camera numbers for each grating camera
	double   linespermm[MAXGRATINGS]; // Lines per mm for each grating camera

};


//===============  Grating parameters associated with a given date/time

struct gratinginfo
{                     
	long     camnum;              // Camera number

    long     year;                //------- Date and UT Time of Calibration (local midnight)
    long     month;
    long     day;
    long     hour;
    long     minute;
    long     second;
	                              //------- Grating information ---------------------------
 	double   grating_linespermm;  // Lines per millimeter
    double   grating_roll;        // Grating rotation around lens axis in radians
    double   grating_pitch;       // Grating normal off lens axis perpendicular to ruled grooves (rad)
    double   grating_yaw;         // Grating normal off lens axis parallel to ruled grooves (rad)
	double   grating_area_scale;  // Grating normal off look direction scaling factor

};


//===============  Extinction and responsivity fit coefficients

struct speccalibration
{
	int      ncams_spcalib;         // Number of cameras with gratings for spectral cal

    struct gratinginfo   *gratinfo;

	int      extn_fittype;        // Extinction fit type (1=SameDay, 2=Multi-Day)
	double   extn_scale;          // Extinction fit leading scale factor
	double   extn_rayleigh;       // Extinction fit coefficient for Rayleigh scattering
	double   extn_oxygen2;        // Extinction fit coefficient for O2 absorption
	double   extn_ozone;          // Extinction fit coefficient for Ozone absorption
	double   extn_dust;           // Extinction fit coefficient for Dust scattering
	double   extn_water;          // Extinction fit coefficient for Water absorption

	double   resp_normalization;  // Responsivity normalization factor in eng units
	double   fading_coef;         // Fading memory coefficient for extinction

								  // ------------- SPCAL file contents
	int       nwavelengths;       // Number of user defined wavelengths
	double   *wavelength_nm;      // user defined wavelengths in nanometers

	double   *wcum_resp_spec;     // multi-day cummulative weight sum
	double   *prev_resp_spec;     // last processed day's responsivity
	double   *cumm_resp_spec;     // multi-day cummulative responsivity
	double   *modl_resp_spec;     // sliding average (smoothed) responsivity
	double   *prev_extn_spec;     // last processed day's extinction
	double   *cumm_extn_spec;     // multi-day fading memory cummulative extinction

								  // ------------- Filled when reading SPCAL
	double   *ord1_resp_spec;     // 1st order responsivity spectrum
	double   *ord2_resp_spec;     // 2nd order responsivity spectrum
	double   *modl_extn_spec;     // model fit of last processed day's extinction

								  // ------------- Work vectors
	double   *esti_resp_spec;     // estimated responsivity ratio
	double   *wsum_resp_spec;     // daily weight sum for responsivity
	double   *aver_resp_spec;     // daily average weighted responsivity
	double   *esti_extn_spec;     // estimated extinction ratio
	double   *wsum_extn_spec;     // daily weight sum for extinction
	double   *aver_extn_spec;     // daily average weighted extinction
	double   *exp_kernel;         // exponential wavelength weighting kernel

	double   *ref_star_spec;      // reference star spectrum from catalog file


};


//===============  Fitting of elements and associated data components

struct  elements_data
{
	int       nelements;          // Number of all meteor element species (neutrals and ions)
	int       kelem_adj;          // Element index to adjust scalable coef
	int       kelem_ref;          // Element index of the abundance reference
	double    sigma0;             // Broadening in nm
	double    Tlo;                // Low excitation Temperature in K
	double    Thi;                // High excitation Temperature in K
	double    plasma_radius_meters;  // Plasma radius at Tlo
	double    range_meteor_meters;   // Range to meteor
	double    VolumeTlo;          // Plasma volume at Tlo in meters^3
	double    VolumeThi;          // Plasma volume at Thi in meters^3
	double    VolumeDepth;        // Plasma volume depth at Tlo in meters
	double    ne_jones;           // Electron density (Jones estimation method)
	double    ne_iter;            // Electron density (Iterative method)
	int       numiter;            // #iterations to reach ne convergence 0.1%
	double    hot2warm;           // Hot to warm plasma contribution ratio
	double    Xfactor;            // Elevation angle factor to scale airmass
	int       nwave;              // Number of desired wavelengths
	double   *wave;               // Desired wavelength in nm
	double   *lockfit;            // Combined element post-fit spectra for locked elements
	double   *corrfit;            // Combined element post-fit spectral flux corrected for saha
	double   *dispfit;            // Combined element post-fit spectral flux display subset

	int       nneutrals;          // Number of neutral only elements
	int      *neutral_index;      // Index to the neutral elements

	int       nfit;               // Number of fitted elements
	int      *fit_element_index;  // Indices of fitted elements

	struct    element_lines_spectrum   *els;  // Per element parameters
};


//===============  Element parameters, emission lines and computed spectrum

struct  element_lines_spectrum
{
    char         element_string[20];    // string describing element species (e.g. Fe-II)
	char         element_filename[20];  // filename containing all emmission lines for this element
	double       elementcode;      // unique number for this element species = atomic#.ionization
	int          ioncode;          // 1=neutral, 2=ion, 51=molecule, 52=molecular ion
	int          ionindex;         // Index of the ion if this is a neutral
	int          neuindex;         // Index of the neutral if this is an ion
	int          warmhot;          // flag for warm=0, hot=1, both=2 contributor
	int          init_fitflag;     // starting flag setting for using element in fit
	int          user_fitflag;     // current flag setting for using element in fit
	int          Npf;              
	double       Tpf[MAXTPF];      // temperatures for partition function
	double       partfunc[MAXTPF]; // partition function values
	double       partfuncTlo;      // partition function at Tlo
	double       partfuncThi;      // partition function at Thi
	double       beta_jones;       // Fraction ionized n+ / (n+ + no) from Jones 1997
	double       init_abundance;   // Starting abundance setting
	double       abundance;        // Actual abundance relative to first fit element
	double       ionenergy;        // Eion
	double       Vo;               // V reference from Jones for ion/neutral ratio
	double       c;                // Coef from Jones for ion/neutral ratio

	double       z_re_Fe;          // Abundance relative to Fe
	double       z_re_Mg;          // Abundance relative to Mg
	double       z_re_Na;          // Abundance relative to Na

	//---------- Parameter vectors across an element's emission lines
	int          nlines;           // Number of emmission lines for this element species
	double      *wave_vacuum;      // Vacuum wavelength in nm
	double      *gf013divw3;       // gf * 0.013248 / w^3
	double      *Eupdiv8;          // -Eupper / 8.617343e-5  (Eupper in eV)

	//---------- Computed spectrum at sigma0 and T
	int          compute_lo;       // Flag to compute Tlo model
	int          compute_hi;       // Flag to compute Thi model
	double      *speclo;           // Spectrum for element at sigma0, Tlo
	double      *spechi;           // Spectrum for element at sigma0, Thi
	double       maxspeclo;        // max spectral value across Tlo wavelength band
	double       maxspechi;        // max spectral value across Thi wavelength band

	//---------- Column densities and ratios for neutral OR ion per "kelem"
	double       N_warm;           // Warm column density (stored in neutral OR ion)
	double       N_hot;            // Hot column density  (stored in neutral OR ion)

	//---------- Column densities =====> STORED ONLY IN THE NEUTRAL ELEMENT'S STRUCTURE
	double       N_warm_default;   // Default warm neutral column density
	double       N_warm_total;     // Warm column density for neutrals + ions
	double       ions2neut_warm;   // # ions / #neutrals ratio warm using Saha
	double       ions2neut_hot;    // # ions / #neutrals ratio hot using Jones
	double       number_atoms;     // # of warm contribution atoms = neutrals + ions
	double       f2fsum_natoms;    // Frame-by-frame running sum of number of atoms

	                               // N_warm(ion)     = N_warm(neutral) * ions2neut_warm
	                               // N_hot (neutral) = N_warm(neutral) * hot2warm * ( 1 + ions2neut_warm ) / ( 1 + ions2neut_hot )
	                               // N_hot (ion)     = N_hot (neutral) * ions2neut_hot
};


//=============== Meteor spectrum, fit spectrum and column density for given frame (or multi-frame combo)

struct  meas_fit_spectra      // Measured and fit are at the user defined wavelength_nm in struct speccalibration
{
	double    height;         // Height of the meteor in km for specified frame
	double   *meas_spectrum;  // Measured, integrated, mean-removed spectrum x #wavelengths
	double   *fit_spectrum;   // Fit spectrum scaled by responsivity and extinction x #wavelengths
	double   *columndensity;  // Final coldensity x #elements (neutrals or ions)
};


//================ Star spectra

struct  StarInfo
{
	long      hip;             // Hipparcos star identifier
	int      userID;          // Reserved for user specific identifier
	char     specType[8];     // Assigned spectral type from Gunn-Stryker catalog
	double   ra_deg;          // RA in degrees (J2000)
	double   dec_deg;         // DEC in degrees (J2000)
	double   vmag;            // V magnitude 
	double   bv;              // B-V magnitude
	double   vi;              // V-I magnitude
	double   bv0;             // B-V (0) magnitude
	double  *specval;         // Vector of spectral values (units?)
	double   specmax;         // Max spectral value
	double   specmin;         // Min spectral value > 0
};

struct  StarSpectraInfo
{
	int      nwave;           // Number of wavelengths in the spectral data file
	int      nstars;          // Number of stars in the data file
	double  *wave;            // Vector of input file spectral wavelengths (nm)
	struct   StarInfo  *star; // Structure to the star information
};



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%        Function Prototypes        %%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//=========================================================================================
//                       Configuration file input function
//=========================================================================================

void     ReadSpectralConfigfile( char                     *pathname,
                                 struct specconfiguration *spconfig  );


//=========================================================================================
//           Spectral calibration file IO functions (responsivity, extinction)
//=========================================================================================

int AllocMemorySpectralCalWavelengths( struct speccalibration   *spcalib,
	                             int                       ncams_spcalib,
	                             int                       nwavelengths, 
	                             double                    startwavelength_nm, 
	                             double                    deltawavelength_nm);

int         ReadSpectralCALfile( char                     *pathname, 
	                             struct specconfiguration *spconfig,
                                 struct speccalibration   *spcalib          );


int        WriteSpectralCALfile( char                     *SPCALfolder, 
								 double                    jdt_midnight, 
                                 struct speccalibration   *spcalib          );

void  FreeMemorySpectralCALfile( struct speccalibration   *spcalib          );


int          GetMostRecentSPCAL( char                     *folder_pathname, 
	                             char                     *file_firstpart, 
								 double                    jdtlast, 
								 char                     *latest_filename );


//=========================================================================================
//                   General purpose spectral file IO functions
//=========================================================================================
void  WriteSpectrum(const char *pathname,
					int nwavelengths,
					double *wavelength,
					double *spectrum1,
					double *spectrum2);


//=========================================================================================
//                      Star Spectra functions
//=========================================================================================

void           ReadStarSpectra( char                      *pathname,  
	                            struct StarSpectraInfo    *starspectra     );

int     GetStarIndexFromSPType( struct StarSpectraInfo    *starspectra,
	                            char                      *specType        );

int     GetStarIndexFromUserID( struct StarSpectraInfo    *starspectra,
	                            int                        userID);

long        GetStarIndexFromHIP( struct StarSpectraInfo    *starspectra, 
	                            long                        HIP             );


void   InterpolateStarSpectrum( struct StarSpectraInfo    *starspectra,
	                            long                       kstar, 
	                            int                        nwave_desired,             
	                            double                    *wavelengths_desired, 
	                            double                    *spectrum_desired  );

void     FreeMemoryStarSpectra( struct StarSpectraInfo    *starspectra     );


//=========================================================================================
//                      Responsivity and Extinction functions
//=========================================================================================

double  ResponsivityExtinctionEstimate( double                  *meas_spectrum,
	                                    double                  *star_spectrum,
	                                    double                   altitude_deg,
	                                    double                   height_observer,
	                                    struct speccalibration  *spcalib );

void            ResponsivityAccumulate( int                      REflag,
	                                    double                   weight,
	                                    struct speccalibration  *spcalib  );

void              ExtinctionAccumulate( int                      REflag,
	                                    double                   weight,
	                                    struct speccalibration  *spcalib  );

double                       Weighting( int                      WEIGHTSflag,
	                                    double                   Vmagnitude );


//=========================================================================================
//                     Elements and Molecules input functions
//=========================================================================================

void              ReadElementsData( int     nwavelengths,
	                                double *wavelength_nm,
						            char   *elements_data_folder,
						            struct  elements_data  *elemdata );

void            FreeElementsMemory( struct  elements_data  *elemdata );

void  AdjustableParametersDefaults( struct  specconfiguration *spconfig,
	                                struct  elements_data     *elemdata);



//=========================================================================================
//                       Support functions for extinction fitting
//=========================================================================================

double         ExtinctionApprox( double                    elev_deg, 
	                             double                    hobs_km, 
								 double                    wavelength_nm );

double                  Airmass( double                    elev_deg, 
	                             double                    rhobs         );

double              RhoIntegral( double                    elev_deg, 
	                             double                    rhobs         );

void         FitExtinctionCoefs( int                       nwave,
	                             double                   *wavelength_nm, 
	                             double                   *extinction_spectrum,
						         struct speccalibration   *spcalibinfo     );
							   
double          TauOpticalDepth( int                       ktype,
	                             double                    wavelength_nm );

void            ExtinctionModel( int                       nwave,
	                             double                   *wavelength_nm, 
	                             double                   *extinction_model,
								 double                    elevation_deg,   //... 90 deg gives unity airmass
								 double                    rhobs,           //... site Height + Rearth
						         struct speccalibration   *spcalibinfo     );

double           AbsorptionLine( double                    wavelength_nm,
	                             double                    center_nm,
					             double                    width_nm,
					             double                    depth_fraction );


//=========================================================================================
//                Measured and fit spectra memory allocation
//=========================================================================================

int         AllocMemorySpectrum( int                       nwavelengths,
	                             int                       nelements,
	                             struct meas_fit_spectra  *meteor_spec );

void         FreeMemorySpectrum( struct meas_fit_spectra  *meteor_spec );



//=========================================================================================
//                      Electron and column density functions
//=========================================================================================

double  PlasmaRadius(double height_km);

void    PlasmaVolumes(double height_km, double range_km, double approach_angle, double hot2warm, struct elements_data *elemdata);

int     GetElementIndex(int atomicNumber, int ionization, struct elements_data *elemdata);

int     GetPrimaryNeutralElement(struct elements_data *elemdata);

double  JonesElectronDensity(struct elements_data *elemdata, int kneutral_primary);

double  IterativeElectronDensity(struct elements_data *elemdata, int kneutral_primary, double ne_guess);

void    ColumnDensities_NumberAtoms(struct elements_data *elemdata, double ne);


//=========================================================================================
//                        Abundance calculation functions
//=========================================================================================

void      ResetOneElementAbundance( int                     kelem,
	                                double                  Vinfinity,
                                    struct  elements_data  *elemdata );

void     ResetAllElementAbundances( double                  Vinfinity, 
	                                struct elements_data   *elemdata );

void      ComputeRelativeAbundance( struct  elements_data  *elemdata );


//=========================================================================================
//                     Spectra calculation and support functions
//=========================================================================================

void             IronOnlySpectrumT( double  T,
									struct  elements_data  *elemdata,
								    double                 *responsivity_1st, 
								    double                 *responsivity_2nd,
									double                  grating_area_scaling,
							        double                 *extinction );

void      ActiveElementsSpectraTlo( struct  elements_data  *elemdata,
								    double                 *responsivity_1st, 
								    double                 *responsivity_2nd, 
									double                  grating_area_scaling,
							        double                 *extinction );

void      ActiveElementsSpectraThi( struct  elements_data  *elemdata,
								    double                 *responsivity_1st, 
								    double                 *responsivity_2nd, 
									double                  grating_area_scaling,
							        double                 *extinction );

void       FitSpectralCoefficients( double                     nwavelengths,
	                                double                    *wavelength_nm,
							        double                    *meas_spectrum,
                                    struct  elements_data     *elemdata,
									struct  specconfiguration *spconfig );

void   SpectrumGivenOnlyLockedCoefs( struct  elements_data     *elemdata,
									double                    *fitspectrum );

void SpectrumGivenOnlyFittingCoefs( struct  elements_data     *elemdata,
									double                    *fitspectrum );

void         SpectrumGivenAllCoefs( struct  elements_data     *elemdata,
									double                    *fitspectrum );
									
void  SpectrumGivenAllCoefs_Subset( struct  elements_data     *elemdata,
									double                    *fitspectrum,
									int                        neutral_ion_display,
									int                        warm_hot_display );
									
void                 FitFeSpectrum( double                     nwavelengths,
	                                double                    *wavelength_nm,
					                double                    *spectrum, 
                                    struct  elements_data     *elemdata,
									struct  specconfiguration *spconfig,
					                double                    *Fe_coef,
					                double                    *Fe_residual );

void                       svd_NxN( double   *A, 
	                                double  *US, 
								    double  *S2, 
								    double   *V, 
								    int       n, 
								    int    maxn );


//=========================================================================================
//                        Julian Date and Calendar functions
//=========================================================================================

double           JulianDateAndTime( long    year, 
	                                long    month, 
	                                long    day,
	                                long    hour, 
	                                long    minute, 
	                                long    second, 
	                                long    milliseconds );

void           CalendarDateAndTime( double  jdt,
		                            long    *year, 
	                                long    *month, 
	                                long    *day,
		                            long    *hour, 
	                                long    *minute, 
	                                long    *second, 
	                                long    *milliseconds );



#endif /* _H_GUARD_SPECTRAL_FUNCTIONS */