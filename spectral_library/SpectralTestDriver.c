//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Test driver for generating line spectra and exercising functionality
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#pragma warning(disable: 4996)  // disable warning on strcpy, fopen, ...


#ifdef _WIN32 /************* WINDOWS ******************************/

	#include <windows.h>    // string fncts, malloc/free, system

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Spectral_Functions.h"
#include "System_FileFunctions.h"


void  WriteSpectrumFile(const char *pathname, int nwavelengths, double *wavelength, double *spectrum1, double *spectrum2);


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[])
{


	struct specconfiguration  spconfig;
	struct speccalibration    spcalib;
	struct elements_data      elemdata;
	struct meas_fit_spectra   spectra, subspectra;
	struct StarSpectraInfo    starspectra;


	//======== Read the spectral configuration file parameters

	ReadSpectralConfigfile("SpectralConfig.txt", &spconfig);

	/*
	printf("config camera number %ld\n", spconfig.camnum[0]);
	printf("config airmass limit %lf\n", spconfig.airmass_limit);
	printf("config hot to warm   %lf\n", spconfig.default_hot2warm);
	printf("config version       %lf\n", spconfig.version);
	*/


	//======== Allocate memory for reading the spectral calibration SPCAL file
	//             which contains the latest responsivity and extinction
	//         Set the number of cameras with gratings
	//         Set the user defined number of processing wavelengths and their start/step in nm

	int     ncameras = 1;

	int     nwave = 1 + (int)((spconfig.max_cal_wavelength_nm - spconfig.min_cal_wavelength_nm) / spconfig.del_wavelength_nm);

	double  startwave = spconfig.min_cal_wavelength_nm;

	double  deltawave = spconfig.del_wavelength_nm;


	AllocMemorySpectralCalWavelengths(&spcalib, ncameras, nwave, startwave, deltawave);

	WriteSpectrumFile("Wave_1.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.wavelength_nm, spcalib.wavelength_nm);



	int camos_camera_index = 0;

	spcalib.gratinfo[camos_camera_index].camnum = 101;


	ReadSpectralCALfile( "SPCAL\\SPCAL_20130421_080000.txt", &spconfig, &spcalib );

	WriteSpectrumFile("Wave_2.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.wavelength_nm, spcalib.wavelength_nm);


	/*
	printf("nwave %ld\n", spcalib.nwavelengths);
	for (int k = 0; k < 501; k++) {
		printf("%lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
			spcalib.wavelength_nm[k],
			spcalib.wcum_resp_spec[k],
			spcalib.prev_resp_spec[k],
			spcalib.cumm_resp_spec[k],
			spcalib.modl_resp_spec[k],
			spcalib.prev_extn_spec[k],
			spcalib.cumm_extn_spec[k]);
	}
	*/


	//========== Read the elements and molecules files

	ReadElementsData(spcalib.nwavelengths, spcalib.wavelength_nm, "ElementsData\\", &elemdata);


    //========== Allocate memeory for measured and fit spectra, column densities
	//           NOTE:  You will need the number of wavelengths and the number
	//                  of elements from ReadElementsData.

	AllocMemorySpectrum(spcalib.nwavelengths, elemdata.nelements, &spectra);

	AllocMemorySpectrum(spcalib.nwavelengths, elemdata.nelements, &subspectra);


	//============================================================================
	//   At this stage, you should extract an integrated spectrum from the 
	//      imagery (for a frame or aggregated frames) that is to be fit and
	//      place the spectrum in the vector "spectrum.integ_spec" at the 
	//      corresponding vector of wavelengths of spcal.wavelength_nm.
	//      NOTE: For the purposes of this example we will infill integ_spec
	//            later with a model spectrum plus noise to be fit.
	//
	//   You will also need the corresponding metadata of the event such as
	//      heights above earth's surface, range to the meteor, the approach 
	//      angle (look direction angle off the radiant), entry velocity, 
	//      altitude angle (elevation above the horizon).   
	//============================================================================

	double  vinfinity_kmsec = 40.0;     //... km/sec

	double  approach_angle_radians = 55.0 * 3.141592654 / 180.0;

	double  earth_radius_km = 6378.16;  //... WGS84

	double  site_height_km = 0.2;       //... Above WGS84 Earth surface

	double  meteor_height_km = 85.0;    //... Above WGS84 Earth surface

	double  altitude_deg = 45.0;        //... elevation angle above horizon

	double  Rsin = earth_radius_km * sin(altitude_deg * 3.141592654 / 180.0);

	double  meteor_range_km = sqrt( Rsin * Rsin + 2.0 * earth_radius_km * meteor_height_km + meteor_height_km * meteor_height_km) - Rsin;


	//========== Set the grating camera with its scaling factor based off the config file input. 
	//           CAMO-S has one grating camera so only index [0] is set. 
	//           On other systems this scaling factor varies with each event due the different angle
	//              of incidence of the light path onto the grating.

	spcalib.gratinfo[camos_camera_index].grating_area_scale = cos(spconfig.grating_offnormal_deg * 3.141592654 / 180.0);


	//========== Set user adjustable values in the elemdata structure to their starting defaults
	//              such as sigma, temperatures, electron density, hot-to-warm, airmass factor

	AdjustableParametersDefaults(&spconfig, &elemdata);


	//========== If any of the input arguments in the call below to PlasmaVolumes changes
	//              at a later time, you must call PlasmaVolumes again to infill the  
	//              elemdata structure values for height, range, and volumes. 
	//              For example, a different frame (height and range) is to be 
	//              processed, or the user adjusts the hot-to-warm ratio for 
	//              fitting purposes.

	PlasmaVolumes(meteor_height_km, meteor_range_km, approach_angle_radians, spconfig.default_hot2warm, &elemdata);


	//========== Compute the model for extinction and the airmass given the event's metadata.
	//           You may want to update for each altitude change or keep it fixed for all frames.

	ExtinctionModel(spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.modl_extn_spec, elemdata.Xfactor * altitude_deg, earth_radius_km + site_height_km, &spcalib);


	//========== Zero all the elemental abundances and #atoms
	//           Set all element fitting flags to not selected for fitting = FITLESS 
	//           Compute Jones 1997 fraction of ionized atoms Beta = n+ / ( n+ + no ) = function( Vinf )

	ResetAllElementAbundances(vinfinity_kmsec, &elemdata);


	//========== Get some element indices into the elements table

	int kelem_Na = GetElementIndex(11, 0, &elemdata);  //... neutral Sodium

	int kelem_Mg = GetElementIndex(12, 0, &elemdata);  //... neutral Magnesium 

	int kelem_Fe = GetElementIndex(26, 0, &elemdata);  //... neutral Iron

	int kelem_C  = GetElementIndex( 6, 0, &elemdata);  //... neutral Carbon


	//========== Obtain model spectra for Iron for warm and hot temperature plasmas

	elemdata.els[kelem_Fe].user_fitflag = FITTING;    //... turn ON iron


	ActiveElementsSpectraTlo( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	ActiveElementsSpectraThi( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec); 
	
	WriteSpectrumFile("FeModelSpectra_LoHi.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Fe].speclo, elemdata.els[kelem_Fe].spechi);


	//========== Obtain new warm temperature model spectra for Iron because Tlo changed

	elemdata.Tlo += 1000.0;  //... add 1000 deg K

	ActiveElementsSpectraTlo( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      


	//========== Obtain new hot temperature model spectra for Iron because Thi changed

	elemdata.Thi -= 2000.0;  //... subtract 2000 deg K

	ActiveElementsSpectraThi( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	WriteSpectrumFile("FeModelSpectra_TempChg.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Fe].speclo, elemdata.els[kelem_Fe].spechi);


	//========== Obtain new plasma volumes and hot temperature model spectra for Iron 
	//           because the hot2warm ratio changed.

	elemdata.hot2warm *= 2.0;  //... doubled the hot-to-warm ratio

	PlasmaVolumes(meteor_height_km, meteor_range_km, approach_angle_radians, elemdata.hot2warm, &elemdata);


	ActiveElementsSpectraThi( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	WriteSpectrumFile("FeModelSpectra_Warm2HotChg.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Fe].speclo, elemdata.els[kelem_Fe].spechi);


	//========== Obtain new model spectra for Iron because the broadening doubled

	elemdata.sigma0 *= 2.0;   //... impacts both Tlo and Thi

	ActiveElementsSpectraTlo( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	ActiveElementsSpectraThi( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	WriteSpectrumFile("FeModelSpectra_SigmaChg.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Fe].speclo, elemdata.els[kelem_Fe].spechi);


	//========== Obtain new model spectra for Iron becasue the airmass factor changed and this 
	//           impacted the extinction model

	elemdata.Xfactor *= 1.2;  //... increased the airmass factor by 20%

	ExtinctionModel(spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.modl_extn_spec, elemdata.Xfactor * altitude_deg, earth_radius_km + site_height_km, &spcalib);

	ActiveElementsSpectraTlo( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	ActiveElementsSpectraThi( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	WriteSpectrumFile("FeModelSpectra_XfactorChg.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Fe].speclo, elemdata.els[kelem_Fe].spechi);


	//========== Obtain new model spectra for both Magnesium and Sodium because the fitting
	//              elements are being changed with Iron removed from the active elements.

	elemdata.els[kelem_Fe].user_fitflag = FITLESS;    //... turn OFF iron

	elemdata.els[kelem_Mg].user_fitflag = FITTING;    //... turn ON magnesium

	elemdata.els[kelem_Na].user_fitflag = FITTING;    //... turn ON sodium

	elemdata.Tlo = spconfig.nominal_lo_exc_temp;
	elemdata.Thi = spconfig.nominal_hi_exc_temp;

	ActiveElementsSpectraTlo( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);		                      

	ActiveElementsSpectraThi( &elemdata,
		                      spcalib.ord1_resp_spec, 
		                      spcalib.ord2_resp_spec, 
		                      spcalib.gratinfo[camos_camera_index].grating_area_scale,
                              spcalib.modl_extn_spec);

	WriteSpectrumFile("MgModelSpectra_LoHi.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Mg].speclo, elemdata.els[kelem_Mg].spechi);

	WriteSpectrumFile("NaModelSpectra_LoHi.txt", spcalib.nwavelengths, spcalib.wavelength_nm, elemdata.els[kelem_Na].speclo, elemdata.els[kelem_Na].spechi);


	//========== Set the reference element for relative abundances

	elemdata.kelem_ref = GetPrimaryNeutralElement(&elemdata);   //... Selects first "active" fitting element in the order of Mg, Fe, Na

	elemdata.kelem_ref = SetPrimaryNeutralElement(12, &elemdata);;  //... OR user specific selection (Mg = atomic number 12)

	//  if kelem_ref returned is -1, then either there is no element in the data base or element is not active (FITTING or FITLOCK)



	//========== Compute the full spectrum which requires column densities for warm, hot, neutrals, ions
	//              Artificially set the neutral warm column densities of the 2 elements
	//              Get a first guess at the electron density from Jones
	//              Iteratively refine the electron density via Jenniskens method
	//              Get the remaining column densities for warm-ions, hot-neutrals, hot-ions
	//              Compute spectrum for all active elements

	elemdata.els[kelem_Mg].N_warm = 3.0e+12;  //... normally get col density from a fit (see later)

	elemdata.els[kelem_Na].N_warm = 2.0e+09;  //... normally get col density from a fit (see later)

	double ne_guess = JonesElectronDensity(&elemdata, elemdata.kelem_ref);          //... sets ne_jones

	printf("ne guess jones  %le\n", ne_guess);

	double ne = IterativeElectronDensity(&elemdata, elemdata.kelem_ref, ne_guess);  //... sets ne_iter

	printf("ne   %le\n", ne);

	ColumnDensities_NumberAtoms(&elemdata, ne, spconfig.beta_flag);

	SpectrumGivenAllCoefs(&elemdata, spectra.fit_spectrum);

	WriteSpectrumFile("MgNaComboSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spectra.fit_spectrum, NULL);

	printf("ne  %le\n", ne);



	//========== Create a measured spectrum from the Mg, Na model with 
	//              added random noise scaled by 10% of the peak spectral magnitude

	double spectrum_max = 0.0;

	for (int kwave = 0; kwave < spcalib.nwavelengths; kwave++) {
		if (spectrum_max < spectra.fit_spectrum[kwave])
			spectrum_max = spectra.fit_spectrum[kwave];
	}

	for (int kwave = 0; kwave < spcalib.nwavelengths; kwave++) {
		double noise_component = 0.00 * spectrum_max * (2.0 * (double)rand() / (double)RAND_MAX - 1.0);
		spectra.meas_spectrum[kwave] = spectra.fit_spectrum[kwave] + noise_component;
	}


	//========== Now fit the measured spectrum to the active elements of Mg and Na.
	//           Note that the active elements were previously set to fitting and
	//               the active element's model spectra have been calculated.
	//           First computes the coefficients for the warm column densities.
	//           Next computes the remaining column densities.

	FitSpectralCoefficients(spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, &elemdata, &spconfig);

	printf("FIT warm column densities:  Mg = %le   Na = %le\n", elemdata.els[kelem_Mg].N_warm, elemdata.els[kelem_Na].N_warm);

	ColumnDensities_NumberAtoms(&elemdata, ne, spconfig.beta_flag);

	SpectrumGivenAllCoefs(&elemdata, spectra.fit_spectrum);

	WriteSpectrumFile("MgNaFitSpectra.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, spectra.fit_spectrum);


	//========== Scale up the warm column density for Mg by 20%

	elemdata.els[kelem_Mg].N_warm *= 1.20;

	ColumnDensities_NumberAtoms(&elemdata, ne, spconfig.beta_flag);

	SpectrumGivenAllCoefs(&elemdata, spectra.fit_spectrum);

	WriteSpectrumFile("MgNaFitSpectra_ChangedMgNwarm.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, spectra.fit_spectrum);


	//========== Spectral fit with FITLOCK on Na but free FITTING on Mg.
	//           Here the Na warm col dens is held fixed, while the Mg warm col dens
	//              is getting fit. You can have multiple elements held fixed and 
	//              multiple elements fit - though a user typically fixes everything 
	//              except a single fitting element.

	elemdata.els[kelem_Na].user_fitflag = FITLOCK; //... no change to Nwarm col dens

	elemdata.els[kelem_Mg].user_fitflag = FITTING; //... fits this element's Nwarm col dens

	FitSpectralCoefficients(spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, &elemdata, &spconfig);

	printf("FIT warm column densities:  Mg = %le   Na = %le\n", elemdata.els[kelem_Mg].N_warm, elemdata.els[kelem_Na].N_warm);

	ColumnDensities_NumberAtoms(&elemdata, ne, spconfig.beta_flag);

	SpectrumGivenAllCoefs(&elemdata, spectra.fit_spectrum);

	WriteSpectrumFile("MgNaFitSpectra_NaFixed.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, spectra.fit_spectrum);


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
	int kelem_xx = GetElementIndex(28, 0, &elemdata);  //... neutral Nickel

	for (int kwave = 0; kwave < spcalib.nwavelengths; kwave++) {
		spcalib.ord1_resp_spec[kwave] = 1.0;
		spcalib.ord2_resp_spec[kwave] = 1.0;
		spcalib.modl_extn_spec[kwave] = 1.0;
	}
*/

	//========== Obtain new model spectra for Iron, Magnesium and Sodium 

	elemdata.els[kelem_Fe].user_fitflag = FITTING;    //... turn ON iron

	elemdata.els[kelem_Mg].user_fitflag = FITTING;    //... turn ON magnesium

	elemdata.els[kelem_Na].user_fitflag = FITTING;    //... turn ON sodium

	elemdata.Tlo = spconfig.nominal_lo_exc_temp;
	elemdata.Thi = spconfig.nominal_hi_exc_temp;

	ActiveElementsSpectraTlo(&elemdata,
		spcalib.ord1_resp_spec,
		spcalib.ord2_resp_spec,
		spcalib.gratinfo[camos_camera_index].grating_area_scale,
		spcalib.modl_extn_spec);

	ActiveElementsSpectraThi(&elemdata,
		spcalib.ord1_resp_spec,
		spcalib.ord2_resp_spec,
		spcalib.gratinfo[camos_camera_index].grating_area_scale,
		spcalib.modl_extn_spec);


	//========== Set the reference element for relative abundances

	elemdata.kelem_ref = GetPrimaryNeutralElement(&elemdata);   //... Selects first active fitting element in the order of Mg, Na, Fe

	elemdata.kelem_ref = SetPrimaryNeutralElement(12, &elemdata);;  //... OR user specific selection (Mg = atomic number 12)


	//========== Compute the full spectrum which requires column densities for warm, hot, neutrals, ions
	//              Artificially set the neutral warm column densities of the 2 elements
	//              Get a first guess at the electron density from Jones
	//              Iteratively refine the electron density via Jenniskens method
	//              Get the remaining column densities for warm-ions, hot-neutrals, hot-ions
	//              Compute spectrum for all active elements

	elemdata.els[kelem_Fe].N_warm = 4.0e+12;  //... normally get col density from a fit (see later)

	elemdata.els[kelem_Mg].N_warm = 3.0e+12;  //... normally get col density from a fit (see later)

	elemdata.els[kelem_Na].N_warm = 2.0e+09;  //... normally get col density from a fit (see later)

	ne_guess = JonesElectronDensity(&elemdata, elemdata.kelem_ref);          //... sets ne_jones

	printf("ne guess jones  %le\n", ne_guess);

	ne = IterativeElectronDensity(&elemdata, elemdata.kelem_ref, ne_guess);  //... sets ne_iter

	printf("ne   %le\n", ne);

	ColumnDensities_NumberAtoms(&elemdata, ne, spconfig.beta_flag);

	SpectrumGivenAllCoefs(&elemdata, spectra.fit_spectrum);

	WriteSpectrumFile("FeMgNaComboSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spectra.fit_spectrum, NULL);

	printf("ne  %le\n", ne);


	//========== Get component by component subspectra for Na 

	elemdata.els[kelem_Fe].user_fitflag = FITTING;    //... turn ON iron
	elemdata.els[kelem_Mg].user_fitflag = FITLESS;    //... turn OFF magnesium
	elemdata.els[kelem_Na].user_fitflag = FITLESS;    //... turn OFF sodium

	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.fit_spectrum, 1, 1);  // Neutral warm
	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.meas_spectrum, 1, 2);  // Neutral hot

	WriteSpectrumFile("FeNeuWarmHotSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, subspectra.fit_spectrum, subspectra.meas_spectrum);

	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.fit_spectrum, 2, 1);  // Ion warm
	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.meas_spectrum, 2, 2);  // Ion hot

	WriteSpectrumFile("FeIonWarmHotSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, subspectra.fit_spectrum, subspectra.meas_spectrum);


	//========== Get component by component subspectra for Na 

	elemdata.els[kelem_Fe].user_fitflag = FITLESS;    //... turn OFF iron
	elemdata.els[kelem_Mg].user_fitflag = FITLESS;    //... turn OFF magnesium
	elemdata.els[kelem_Na].user_fitflag = FITTING;    //... turn ON sodium

	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.fit_spectrum, 1, 1);  // Neutral warm
	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.meas_spectrum, 1, 2);  // Neutral hot

	WriteSpectrumFile("NaNeuWarmHotSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, subspectra.fit_spectrum, subspectra.meas_spectrum);

	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.fit_spectrum, 2, 1);  // Ion warm
	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.meas_spectrum, 2, 2);  // Ion hot

	WriteSpectrumFile("NaIonWarmHotSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, subspectra.fit_spectrum, subspectra.meas_spectrum);


	//========== Get component by component subspectra for Mg  (Mg and Na were only two turned on at this point)

	elemdata.els[kelem_Fe].user_fitflag = FITLESS;    //... turn OFF iron
	elemdata.els[kelem_Mg].user_fitflag = FITTING;    //... turn ON magnesium
	elemdata.els[kelem_Na].user_fitflag = FITLESS;    //... turn OFF sodium

	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.fit_spectrum, 1, 1);  // Neutral warm
	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.meas_spectrum, 1, 2);  // Neutral hot

	WriteSpectrumFile("MgNeuWarmHotSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, subspectra.fit_spectrum, subspectra.meas_spectrum);

	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.fit_spectrum, 2, 1);  // Ion warm
	SpectrumGivenAllCoefs_Subset(&elemdata, subspectra.meas_spectrum, 2, 2);  // Ion hot

	WriteSpectrumFile("MgIonWarmHotSpectra_KnownColDens.txt", spcalib.nwavelengths, spcalib.wavelength_nm, subspectra.fit_spectrum, subspectra.meas_spectrum);


	elemdata.els[kelem_Fe].user_fitflag = FITTING;    //... turn ON iron
	elemdata.els[kelem_Mg].user_fitflag = FITTING;    //... turn ON magnesium
	elemdata.els[kelem_Na].user_fitflag = FITTING;    //... turn ON sodium



	//========== Create a measured spectrum from the Mg, Na model with 
	//              added random noise scaled by 10% of the peak spectral magnitude

	spectrum_max = 0.0;

	for (int kwave = 0; kwave < spcalib.nwavelengths; kwave++) {
		if (spectrum_max < spectra.fit_spectrum[kwave])
			spectrum_max = spectra.fit_spectrum[kwave];
	}

	for (int kwave = 0; kwave < spcalib.nwavelengths; kwave++) {
		double noise_component = 0.00 * spectrum_max * (2.0 * (double)rand() / (double)RAND_MAX - 1.0);
		spectra.meas_spectrum[kwave] = spectra.fit_spectrum[kwave] + noise_component;
	}


	//========== Now fit the measured spectrum to the active elements of Mg and Na.
	//           Note that the active elements were previously set to fitting and
	//               the active element's model spectra have been calculated.
	//           First computes the coefficients for the warm column densities.
	//           Next computes the remaining column densities.

	FitSpectralCoefficients(spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, &elemdata, &spconfig);

	printf("FIT warm column densities:  Fe = %le   Mg = %le   Na = %le\n", elemdata.els[kelem_Fe].N_warm, elemdata.els[kelem_Mg].N_warm, elemdata.els[kelem_Na].N_warm);

	ColumnDensities_NumberAtoms(&elemdata, ne, spconfig.beta_flag);

	SpectrumGivenAllCoefs(&elemdata, spectra.fit_spectrum);

	WriteSpectrumFile("FeMgNaFitSpectra.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spectra.meas_spectrum, spectra.fit_spectrum);




	//========== Compute and report on relative abundances

	ComputeRelativeAbundance(&elemdata);

	for (int kneutral = 0; kneutral < elemdata.nneutrals; kneutral++) {

		int kelem = elemdata.neutral_index[kneutral];

		if (elemdata.els[kelem].user_fitflag == FITLESS) continue; //... do not report non-active elements

		printf("Relative Abundance %s / %s = %lf\n", elemdata.els[kelem].element_string, elemdata.els[elemdata.kelem_ref].element_string, elemdata.els[kelem].abundance);

	}




	//=========================================================================================
	//  Example showing responsivity and extinction updates based on known star spectra
	//=========================================================================================

	//---------- Plot files for existing responsivity and extinction

	WriteSpectrumFile("ResponsivityExisting.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.prev_resp_spec, spcalib.cumm_resp_spec);

	WriteSpectrumFile("ExtinctionExisting.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.prev_extn_spec, spcalib.cumm_extn_spec);

	FitExtinctionCoefs(spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.cumm_extn_spec, &spcalib);

	ExtinctionModel(spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.modl_extn_spec, elemdata.Xfactor * altitude_deg, earth_radius_km + site_height_km, &spcalib);

	WriteSpectrumFile("RespExtnModelExisting.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.modl_resp_spec, spcalib.modl_extn_spec);


	//========== Read the star catalog, get the index for the star Canopus, then show function
	//             call to interpolate its cataloged spectrum to the user defined wavelengths.
	
	ReadStarSpectra("StarSpectra_V5.0_RA_0_360_DEC_56_90.txt", &starspectra);

	int hip = 11767;

	int star_index = (int)GetStarIndexFromHIP(&starspectra, (long)hip);

	InterpolateStarSpectrum(&starspectra, star_index, spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.ref_star_spec);

	printf("Hipparcos %ld is index %ld\n", hip, star_index);

	WriteSpectrumFile("StarSpectrum.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.ref_star_spec, spcalib.ref_star_spec);

	Delay_msec(10000);

	//========== Clear the responsivity and extinction accumulator vectors

	ResponsivityAccumulate(RECLEAR, 1.0, &spcalib);

	ExtinctionAccumulate(RECLEAR, 1.0, &spcalib);


	//========== Loop over e.g. three stars to update the responsivity and extinction for the night
	//           Find the star's index by its name (we will just use the first three indices 0, 1, 2)

	int starcounter = 0;

	while (starcounter < 3) {

		star_index = starcounter;  //... Normally get this index from the star name collected


		//-------- Interpolate this cataloged star's spectrum and place in the spcalib structural
		//         element ref_star_spec at the user specified wavelengths via linear interpolation.

		InterpolateStarSpectrum(&starspectra, star_index, spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.ref_star_spec);


		//-------- Generate an "artificial" star measurement by convolving the star with the 
		//         current 1st order responsivity, latest extinction model and added random noise,
		//         at the user specified wavelengths.
		//         NOTE: Normally you would extract the star spectrum from a star collection.

		for (int kwave = 0; kwave < spcalib.nwavelengths; kwave++) {
			double noise_component = 0.00 * starspectra.star[star_index].specmax * (2.0 * (double)rand() / (double)RAND_MAX - 1.0);
			spectra.meas_spectrum[kwave] = spcalib.ref_star_spec[kwave] * spcalib.ord1_resp_spec[kwave] * spcalib.modl_extn_spec[kwave] + noise_component;
		}


		//-------- Compute the individual star responsivity and extinction estimate and let user decide
		//         to add it to the cummulative running results or skip to next star.

		double  Xairmass = ResponsivityExtinctionEstimate(spectra.meas_spectrum, spcalib.ref_star_spec, altitude_deg, site_height_km, &spcalib);


		//-------- Compute the weighting function

		double  weight = Weighting( UNITYWEIGHTS, 0.0);  //... uses unity weighting


		//.... Display the estimated responsivity and extinction over the ord1_resp and modl_extn

		    //**** Plot the display of responsivity and extinction ****


		//.... IF airmass below the configuration file limit and ...
		//     IF the user decides to include the estimated responsivity to the daily running response

	    if(Xairmass <= spconfig.airmass_limit)  ResponsivityAccumulate( REDAILY, weight, &spcalib );


		//.... IF the user decides to include the estimated extinction to the daily running extinction tally

	 	ExtinctionAccumulate( REDAILY, weight, &spcalib );


        //-------- Go to next star spectrum collected

		starcounter++;

	} //... end of star collection loop for the night


    //========== Display the estimated responsivity and extinction (ord1_resp and modl_resp)
	//           The user should decide to save the updated results or pass on
	//               writing/updating the SPCAL file for this data if it was of poor quality.
	//           Note the SPCAL file has the date/time in the filename to make it unique
	//               and not overwrite the previous SPCAL.

		//.... IF the user decides to include the daily responsivity to the multi-day cummulative response

    ResponsivityAccumulate(REMULTI, 1.0, &spcalib);

	    //.... IF the user decides to include the daily extinction to the multi-day cummulative extinction

    ExtinctionAccumulate( REMULTI, 1.0, &spcalib);


	//========== Fit the extinction model to the 5 known contributors to update the extinction coefficients

	FitExtinctionCoefs(spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.cumm_extn_spec, &spcalib);


    //========== Write updated responsivity and extinction to a file

	double  jdt = JulianDateAndTime(2013, 4, 30, 8, 0, 0, 0);

	WriteSpectralCALfile("SPCAL\\", jdt, &spcalib);


	//---------- Plot files for updated responsivity and extinction

	WriteSpectrumFile("ResponsivityUpdated.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.prev_resp_spec, spcalib.cumm_resp_spec);

	WriteSpectrumFile("ExtinctionUpdated.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.prev_extn_spec, spcalib.cumm_extn_spec);


	//---------- Plot files for model fit extinction

	ExtinctionModel(spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.modl_extn_spec, elemdata.Xfactor * altitude_deg, earth_radius_km + site_height_km, &spcalib);

	WriteSpectrumFile("RespExtnModelFit.txt", spcalib.nwavelengths, spcalib.wavelength_nm, spcalib.modl_resp_spec, spcalib.modl_extn_spec);


	//=========================================================================================
	//   Special test section on extinction components
	//=========================================================================================
	
	/*
	FILE *outfile;
	outfile = fopen("ExtinctionComponents.txt", "wt");
	for (int kwave = 0; kwave < 1001; kwave++) {
		double wavelength_nm = 350.0 + ((double)kwave) / 2.0;
		fprintf(outfile, "%lf  %lf  %lf  %lf  %lf  %lf  %lf\n", wavelength_nm,
			TauOpticalDepth(RAYLEIGH, wavelength_nm),
			TauOpticalDepth(O2BAND1,  wavelength_nm),
			TauOpticalDepth(O2BAND2,  wavelength_nm),
			TauOpticalDepth(OZONE,    wavelength_nm),
			TauOpticalDepth(DUST,     wavelength_nm),
			TauOpticalDepth(WATER,    wavelength_nm));
	}
	fclose(outfile);
	*/

	//=========================================================================================
	//   Special test section on plasma radius versus height and range
	//=========================================================================================
	
	/*
	int nheights = 120;
	for (int kheight_km = 1; kheight_km <= nheights; kheight_km++) {
		double singamma = 6378.16 * sin((90 + 45) / 57.296) / (6378.16 + (double)kheight_km);
		double beta = 180 - 90 - 45 - 57.296*asin(singamma);
		double range = 6378.16 * sin(beta / 57.296) / singamma;

		spcalib.wavelength_nm[kheight_km] = (double)kheight_km;
		spcalib.aver_resp_spec[kheight_km] = range;
		spcalib.aver_extn_spec[kheight_km] = PlasmaRadius((double)kheight_km);
	}

	WriteSpectrumFile("PlasmaRadius.txt", nheights, spcalib.wavelength_nm, spcalib.aver_resp_spec, spcalib.aver_extn_spec);
	*/

	//=========================================================================================
	//   Cleanup - Free the allocated memory
	//=========================================================================================

	FreeMemorySpectrum(&spectra);

	FreeElementsMemory(&elemdata);

	FreeMemorySpectralCALfile(&spcalib);

	FreeMemoryStarSpectra(&starspectra);


	printf("Done\n");

	Delay_msec(30000);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void  WriteSpectrumFile(const char *pathname, int nwavelengths, double *wavelength, double *spectrum1, double *spectrum2)
{
	FILE *outfile;

	if ((outfile = fopen(pathname, "wt")) == NULL) {
		printf(" Cannot open output spectra file %s for writing\n", pathname);
		return;
	}

	fprintf(outfile, "%ld\n", nwavelengths);
	
	for (int kwave = 0; kwave < nwavelengths; kwave++) {

		if( spectrum2 == NULL ) fprintf(outfile, "%8.3lf  %15.6le\n", wavelength[kwave], spectrum1[kwave] );
		else                    fprintf(outfile, "%8.3lf  %15.6le  %15.6le\n", wavelength[kwave], spectrum1[kwave], spectrum2[kwave]);
	
	}

	fclose(outfile);

}