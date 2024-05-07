""" 
	Sels
	pectral library 
	Interfaces Pete Gural's C library with python
	Built for use with CAMO-S analysis

"""


#######################################################################################
####################################### IMPORTS #######################################
#######################################################################################


#from _typeshed import Self
import os
import sys
import time
import math
import numpy as np
import numpy.ctypeslib as npct
import ctypes as ct
import matplotlib.pyplot as plt


#######################################################################################
############################## PATHS, C TYPES, CONSTANTS ##############################
#######################################################################################

def adjustableParametersDefaults(self):
	self.spectral_lib.AdjustableParametersDefaults(self.spconfig, self.elemdata)

def loadSpectralLibrary(self):
	"""
	Loads the SPCAL folder with necessary calibration  files.
	"""

	# Path to the spectral library
	SPECTRAL_LIBRARY = os.path.join("spectral_library", "SpectralTest")

	# Load the spectral library 
	self.spectral_lib = npct.load_library(SPECTRAL_LIBRARY, os.path.dirname(__file__))


def loadElementsData(self):
	"""
	Loads ElementsData file to extract element information.
	"""

	# Load elements data
	self.spectral_lib.ReadElementsData(self.spcalib.nwavelengths, self.spcalib.wavelength_nm, "spectral_library/ElementsData/".encode('ascii'), self.elemdata)


def readSpectralConfig(self):
	"""
	Reads in the spectral configuration file.
	"""

	# Path to spectral config file
	spectral_config_path = os.path.join("spectral_library", "DriverInputFiles", "SpectralConfig.txt")

	# Read spectral configuration file
	self.spectral_lib.ReadSpectralConfigfile(os.path.join(os.path.dirname(__file__), \
			spectral_config_path).encode('ascii'), self.spconfig)

	self.ncameras = 1
	self.nwave = 1 + int(((self.spconfig.max_cal_wavelength_nm - self.spconfig.min_cal_wavelength_nm) / self.spconfig.del_wavelength_nm))
	self.startwave = self.spconfig.min_cal_wavelength_nm
	self.deltawave = self.spconfig.del_wavelength_nm


def allocMemory(self):
	"""
	Allocates memeory for the calibration wavelengths and the spectrum.
	"""

	# Allocate memory for the spectral calibration wavelengths
	self.spectral_lib.AllocMemorySpectralCalWavelengths(self.spcalib, self.ncameras, self.nwave, self.startwave, self.deltawave)

	# Allocate memory for the spectrum
	self.spectral_lib.AllocMemorySpectrum(self.spcalib.nwavelengths, self.elemdata.nelements, self.spectra) 


def readSpectralCALFile(self):
	"""
	Reads in spectral calibration file
	"""

	# Read in the spectral calibration file
	# self.spectral_lib.ReadSpectralCALfile("spectral_library/SPCAL/SPCAL_20130421_080000.txt".encode('ascii'), self.spconfig, self.spcalib)
	self.spectral_lib.ReadSpectralCALfile("spectral_library/SPCAL/SPCAL_20240503_080000.txt".encode('ascii'), self.spconfig, self.spcalib)


def readStarSpectra(self):
	"""
	Reads in star spectra file.
	"""

	# Read in the star spectra 
	self.spectral_lib.ReadStarSpectra("spectral_library/DriverInputFiles/StarSpectra_V5.0_RA_0_360_DEC_56_90.txt".encode('ascii'), self.starspectra)


# Init ctypes types
DOUBLE = ct.c_double
PDOUBLE = ct.POINTER(DOUBLE)
PPDOUBLE = ct.POINTER(PDOUBLE)
PPPDOUBLE = ct.POINTER(PPDOUBLE)
INT = ct.c_int 
LONG = ct.c_long 
VOID = ct.c_voidp
CHAR = ct.c_char
STRUCTURE = ct.Structure

### Define constants ###

# Max number of cameras with gratings
MAXGRATINGS =    32    

# Max # of temperatures in partition functions file
MAXTPF =         20    

# Warmhot flag
WARM =            0    
HOT =             1
BOTH =            2

# User_fitflag
FITLESS =         0    
FITTING =         1
FITLOCK =         2

# Extinction components
RAYLEIGH =        1    
O2BAND1 =         2
O2BAND2 =         3
OZONE =           4
DUST =            5
WATER =           6

# REflag = Responsivity and extinction control flag
RECLEAR =         0   
REDAILY =         1
REMULTI =         2

# Weighting schema for responsivity and extinction averaging
UNITYWEIGHTS =    1   
VMAGWEIGHTS =     2

# Earth's radius in km, WGS84
EARTH_RADIUS_KM = 6378.16 


#######################################################################################
################################ DEFINE PYTHON CLASSES ################################
#######################################################################################



############### Spectral processing configuration file input parameters ###############
class SpecConfiguration(STRUCTURE):
	"""
	Structure containing spectral processing configuration file, input parameters.
	Defines all associated variable times.
	"""

	_fields_ = [
	
		# Config file version number
		("version", DOUBLE),

		# Order to use for spectral calibration
		("order4spcalib", LONG),

		# Plus/minus #rows for integrating spectra in pixels
		("rowdelta", LONG),

		# Minimum wavelength to calibrate response in nm
		("min_cal_wavelength_nm", DOUBLE),

		# Maximum wavelength to calibrate response in nm
		("max_cal_wavelength_nm", DOUBLE),

		# Stepping wavelength to calibrate response in nm
		("del_wavelength_nm", DOUBLE),

		# Minimum fitting wavelength in nm
		("min_fit_wavelength_nm", DOUBLE),

		# Maximum fitting wavelength in nm
		("max_fit_wavelength_nm", DOUBLE),

		# Min bandwidth the order4spcalib spectrum must span to use
		("minspan_wavelength_nm", DOUBLE),

		# 1 sigma width of Gaussian smoothing filter
		("smooth_bandwidth_nm", DOUBLE),

		# Faintest star to process in mV
		("faintest_star_vmag", DOUBLE),

		# Airmass value for inclusion in response or extinction
		("airmass_limit", DOUBLE),

		# Fading memory coefficient
		("fading_coef", DOUBLE),

		# Time tolerance for coincidence in seconds
		("coin_time_tolerance", DOUBLE),

		# Minimum low excitation temperature in deg K
		("min_lo_exc_temp", DOUBLE),

		# Maximum low excitation temperature in deg K
		("max_lo_exc_temp", DOUBLE),

		# Step low excitation temperature in deg K
		("step_lo_exc_temp", DOUBLE),

		# Nominal low excitation temperature in deg K
		("nominal_lo_exc_temp", DOUBLE),

		# Nominal high excitation temperature in deg K
		("nominal_hi_exc_temp", DOUBLE),

		# Nominal broadening in nm
		("nominal_sigma0", DOUBLE),

		# Grating normal to look direction angle in deg
		("grating_offnormal_deg", DOUBLE),

		# Default grating roll angle in deg
		("default_roll_deg", DOUBLE),

		# Default grating picth angle in deg
		("default_pitch_deg", DOUBLE),

		# Default grating yaw angle in deg
		("default_yaw_deg", DOUBLE),

		# Default density of electrons
		("default_ne", DOUBLE),

		# Default hot to warm plasma contribution ratio
		("default_hot2warm", DOUBLE),

		# Number of cameras with gratings
		("ncams_grating", LONG),

		# Camera numbers for each grating camera
		("camnum", MAXGRATINGS*LONG),

		# Lines per mm for each grating camera
		("linespermm", MAXGRATINGS*DOUBLE)

	]


################ Grating parameters associated with a given date/time #################
class GratingInfo(STRUCTURE):
	"""
	Structure containing grating parameters. 
	Defines all associated variable types.
	"""

	_fields_ = [
		
		# Camera number
		("camnum", LONG),              

		 #------- Date and UT Time of Cam ("local" midnight)
		("year", LONG),
		("month", LONG),
		("day", LONG),
		("hour", LONG),
		("minute", LONG),
		("second", LONG),
		
		#------- Grating information ---------------------------

		# Lines per millimeter
		("grating_linespermm", DOUBLE), 

		# Grating rotation around lens axis in radians
		("grating_roll", DOUBLE),        

		# Grating normal off lens axis perpendicular to ruled grooves ("rad)
		("grating_pitch", DOUBLE),       

		# Grating normal off lens axis parallel to ruled grooves ("rad)
		("grating_yaw", DOUBLE),         

		# Grating normal off look direction scaling factor
		("grating_area_scale", DOUBLE)

	]


#################### Extinction and responsivity fit coefficients #####################
class SpecCalibration(STRUCTURE):
	"""
	Structure containing extinction and responsivity fit coeffiecients. 
	Defines all associated variable types.
	"""

	_fields_ = [

		# Number of cameras with gratings for spectral cal
		("ncams_spcalib", INT),         
		("gratinfo", ct.POINTER(GratingInfo)),

		# Extinction fit type (1=SameDay, 2=Multi-Day)
		("extn_fittype", INT),      

		# Extinction fit leading scale factor
		("extn_scale", DOUBLE),          

		# Extinction fit coefficient for Rayleigh scattering
		("extn_rayleigh", DOUBLE),       

		# Extinction fit coefficient for O2 absorption
		("extn_oxygen2", DOUBLE),        

		# Extinction fit coefficient for Ozone absorption
		("extn_ozone", DOUBLE),          

		# Extinction fit coefficient for Dust scattering
		("extn_dust", DOUBLE),           

		# Extinction fit coefficient for Water absorption
		("extn_water", DOUBLE),          

		# Responsivity normalization factor in eng units
		("resp_normalization", DOUBLE),  
	
		# Fading memory coefficient for extinction
		("fading_coef", DOUBLE),         
		

		# ------------- SPCAL file contents

		# Number of user defined wavelengt
		("nwavelengths", INT),    

		# user defined wavelengths in nanometers
		("wavelength_nm", PDOUBLE),      

		# multi-day cummulative weight sum
		("wcum_resp_spec", PDOUBLE),     

		# last processed day's responsivity
		("prev_resp_spec", PDOUBLE),     

		# multi-day cummulative responsivity
		("cumm_resp_spec", PDOUBLE),     

		# sliding average (smoothed) responsivity
		("modl_resp_spec", PDOUBLE),     

		# last processed day's extinction
		("prev_extn_spec", PDOUBLE),     

		# multi-day fading memory cummulative extinction
		("cumm_extn_spec", PDOUBLE),     

		# ------------- Filled when reading SPCAL

		# 1st order responsivity spectrum
		("ord1_resp_spec", PDOUBLE),     

		# 2nd order responsivity spectrum
		("ord2_resp_spec", PDOUBLE),     

		# model fit of last processed day's extinction
		("modl_extn_spec", PDOUBLE),     


		# ------------- Work vectors

		# estimated responsivity ratio
		("esti_resp_spec", PDOUBLE),

		# daily weight sum for responsivity     
		("wsum_resp_spec", PDOUBLE),

		# daily average weighted responsivity     
		("aver_resp_spec", PDOUBLE),

		# estimated extinction ratio     
		("esti_extn_spec", PDOUBLE),

		# daily weight sum for extinction     
		("wsum_extn_spec", PDOUBLE),

		# daily average weighted extinction     
		("aver_extn_spec", PDOUBLE),

		# exponential wavelength weighting kernel     
		("exp_kernel", PDOUBLE),    

		# reference star spectrum from catalog file
		("ref_star_spec", PDOUBLE)      

	]


############## Element parameters, emission lines and computed spectrum ###############
class ElementLinesSpectrum(STRUCTURE):
	"""
	Structure containing all element parameters with regards to emission lines/spectra.
	Defines all associated variable types.
	"""

	_fields_ = [
	
		# string describing element species (e.g. Fe-II)
		("element_string", 20*CHAR),    

		# filename containing all emmission lines for this element
		("element_filename", 20*CHAR),  

		# unique number for this element species = atomic num.ionization
		("elementcode", DOUBLE),      

		  # 1=neutral, 2=ion, 51=molecule, 52=molecular ion
		("ioncode", INT),        

		# Index of the ion if this is a neutral
		("ionindex", INT),         

		# Index of the neutral if this is an ion
		("neuindex", INT),         

		# flag for warm=0, hot=1, both=2 contributor
		("warmhot", INT),          

		# starting flag setting for using element in fit
		("init_fitflag", INT),     

		# current flag setting for using element in fit
		("user_fitflag", INT),     

		("Npf", INT), 

		# temperatures for partition function
		("Tpf", MAXTPF*DOUBLE),      

		# partition function values
		("partfunc", MAXTPF*DOUBLE), 

		# partition function at Tlo
		("partfuncTlo", DOUBLE),      

		# partition function at Thi
		("partfuncThi", DOUBLE),      

		# Fraction ionized n+ / (n+ + no) from Jones 1997
		("beta_jones", DOUBLE),       

		# Starting abundance setting
		("init_abundance", DOUBLE),   

		# Actual abundance relative to first fit element
		("abundance", DOUBLE),        

		# Eion
		("ionenergy", DOUBLE),        

		# V reference from Jones for ion/neutral ratio
		("Vo", DOUBLE),               

		# Coef from Jones for ion/neutral ratio
		("c", DOUBLE),                

		# Abundance relative to Fe
		("z_re_Fe", DOUBLE),        

		# Abundance relative to Mg
		("z_re_Mg", DOUBLE),         

		# Abundance relative to Na
		("z_re_Na", DOUBLE),          
		
		#---------- Parameter vectors across an element's emission lines

		# Number of emmission lines for this element species
		("nlines", INT),          
		
		# Vacuum wavelength in nm
		("wave_vacuum", PDOUBLE),      

		# gf * 0.013248 / w^3
		("gf013divw3", PDOUBLE),       

		# -Eupper / 8.617343e-5  (Eupper in eV)
		("Eupdiv8", PDOUBLE),          

		#---------- Computed spectrum at sigma0 and T

		# Flag to compute Tlo model
		("compute_lo", INT),       
		
		# Flag to compute Thi model
		("compute_hi", INT),       
		
		# Spectrum for element at sigma0, Tlo
		("speclo", PDOUBLE),           
		
		# Spectrum for element at sigma0, Thi
		("spechi", PDOUBLE),           
		
		# max spectral value across Tlo wavelength band
		("maxspeclo", DOUBLE),        
		
		# max spectral value across Thi wavelength band
		("maxspechi", DOUBLE),        
		
		#---------- Column densities and ratios for neutral OR ion per "kelem"

		# Warm column density (stored in neutral OR ion)
		("N_warm", DOUBLE),           

		# Hot column density  (stored in neutral OR ion)
		("N_hot", DOUBLE),            

		#---------- Column densities =====> STORED ONLY IN THE NEUTRAL ELEMENT'S STRUCTURE

		# Default warm neutral column density
		("N_warm_default", DOUBLE), 

		# Warm column density for neutrals + ions  
		("N_warm_total", DOUBLE),   

		# num ions / num neutrals ratio warm using Saha  
		("ions2neut_warm", DOUBLE), 

		# num ions / num neutrals ratio hot using Jones  
		("ions2neut_hot", DOUBLE),  

		# num of warm contribution atoms = neutrals + ions  
		("number_atoms", DOUBLE),   

		# Frame-by-frame running sum of number of atoms  
		("f2fsum_natoms", DOUBLE)  


									# N_warm(ion)     = N_warm(neutral) * ions2neut_warm
									# N_hot (neutral) = N_warm(neutral) * hot2warm * ( 1 + ions2neut_warm ) / ( 1 + ions2neut_hot )
									# N_hot (ion)     = N_hot (neutral) * ions2neut_hot
	]


################# Fitting of elements and associated data components ##################
class ElementsData(STRUCTURE):
	"""
	Structure containing all elements data relevant to fitting.
	Defines all associated variable types.
	"""
	
	_fields_ = [

		# Number of all meteor element species (neutrals and ions)
		("nelements", INT),          

		# Element index to adjust scalable coef
		("kelem_adj", INT),          

		# Element index of the abundance reference
		("kelem_ref", INT),          

		# Broadening in nm
		("sigma0", DOUBLE),               

		# Low excitation Temperature in K
		("Tlo", DOUBLE),                  

		# High excitation Temperature in K
		("Thi", DOUBLE),                  

		# Plasma radius at Tlo
		("plasma_radius_meters", DOUBLE), 

		# Range to meteor
		("range_meteor_meters", DOUBLE),  

		# Plasma volume at Tlo in meters^3
		("VolumeTlo", DOUBLE),            

		# Plasma volume at Thi in meters^3
		("VolumeThi", DOUBLE),            

		# Plasma volume depth at Tlo in meters
		("VolumeDepth", DOUBLE),          

		# Electron density (Jones estimation method)
		("ne_jones", DOUBLE),             

		# Electron density (Iterative method)
		("ne_iter", DOUBLE),              

		# iterations to reach ne convergence 0.1%
		("numiter", INT),            

		# Hot to warm plasma contribution ratio
		("hot2warm", DOUBLE),             

		# Elevation angle factor to scale airmass
		("Xfactor", DOUBLE),              

		# Number of desired wavelengths
		("nwave", INT),              

		# Desired wavelength in nm
		("wave", PDOUBLE),                

		# Combined element post-fit spectra for locked elements
		("lockfit", PDOUBLE),             

		# Combined element post-fit spectral flux corrected for saha
		("corrfit", PDOUBLE),             

		# Combined element post-fit spectral flux display subset
		("dispfit", PDOUBLE),             

		# Number of neutral only elements
		("nneutrals", INT),          

		# Index to the neutral elements
		("neutral_index", ct.POINTER(INT)),      

		# Number of fitted elements
		("nfit", INT),     

		# Indices of fitted elements
		("fit_element_index", ct.POINTER(INT)),  

		# Per element parameters
		("els", ct.POINTER(ElementLinesSpectrum))  

	]


######## Meteor spec; fit spectrum/column density for frame, multi-frame combo ########
class MeasFitSpectra(STRUCTURE): 
	"""
	Measured and fit are at the user defined wavelength_nm, in struct speccalibration.
	Defines all associated variable types.
	"""
	 
	_fields_ = [

		# Height of the meteor in km for specified frame
		("height", DOUBLE),         

		# Measured, integrated, mean-removed spectrum x num_wavelengths
		("meas_spectrum", PDOUBLE),  

		# Fit spectrum scaled by responsivity and extinction x num_wavelengths
		("fit_spectrum", PDOUBLE),  

		# Final coldensity x num_elements (neutrals or ions) 
		("columndensity", PDOUBLE)  
	]



################################## Star spectra #######################################



class StarInfo(STRUCTURE):
	"""
	Structure containing star information.
	Defines all associated variable types.
	"""

	_fields_ = [

		# Hipparcos star identifier
		("hip", LONG),             

		# Reserved for user specific identifier
		("userID", INT),          

		# Assigned spectral type from Gunn-Stryker catalog
		("specType", 8*CHAR),     

		# RA in degrees (J2000)
		("ra_deg", DOUBLE),          

		# DEC in degrees (J2000)
		("dec_deg", DOUBLE),         

		# V magnitude 
		("vmag", DOUBLE),            

		# B-V magnitude
		("bv", DOUBLE),              

		# V-I magnitude
		("vi", DOUBLE),              

		# B-V (0) magnitude
		("bv0", DOUBLE),             

		# Vector of spectral values (units?)
		("specval", PDOUBLE),         

		# Max spectral value
		("specmax", DOUBLE),       

		# Min spectral value > 0
		("specmin", DOUBLE),       

	]


class StarSpectraInfo(STRUCTURE):
	"""
	Structure containing spectral star information.  
	Defines all associated variable types.
	"""
	
	_fields_ = [

		# Number of wavelengths in the spectral data file
		("nwave", INT),           

		# Number of stars in the data file
		("nstars", INT),          

		# Vector of input file spectral wavelengths (nm)
		("wave", PDOUBLE),            

		# Structure to the star information
		("star", ct.POINTER(StarInfo)) 

	]



#######################################################################################
################################## SUPPORT FUNCTIONS ##################################
#######################################################################################


def double2ArrayToPointer(arr):
	""" Converts a 2D numpy to ctypes 2D array. 
	
	Arguments:
		arr: [ndarray] 2D numpy float64 array
	Return:
		arr_ptr: [ctypes double pointer]
	"""

	# Init needed data types
	ARR_DIMX = DOUBLE*arr.shape[1]
	ARR_DIMY = PDOUBLE*arr.shape[0]

	# Init pointer
	arr_ptr = ARR_DIMY()

	# Fill the 2D ctypes array with values
	for i, row in enumerate(arr):
		arr_ptr[i] = ARR_DIMX()

		for j, val in enumerate(row):
			arr_ptr[i][j] = val


	return arr_ptr


def double1pointerToArray(ptr, n):
	""" Converts ctypes 1D array into a 1D numpy array. 
	
	Arguments:
		ptr: [ctypes double pointer]
		n: [int] number of cameras
	Return:
		arr: [ndarrays] converted numpy array
		
	"""

	# Init a new empty data array
	arr = np.zeros(shape=n)

	# Go through every camera
	for i in range(n):
		arr[i] = ptr[i]

	return arr


def double2pointerToArray(ptr, n, m_sizes):
	""" Converts ctypes 2D array into a 2D numpy array. 
	
	Arguments:
		ptr: [ctypes double pointer]
		n: [int] number of cameras
		m_sizes: [list] number of measurements for each camera
	Return:
		arr_list: [list of ndarrays] list of numpy arrays, each list entry containing data for individual
			cameras
		
	"""

	arr_list = []

	# Go through every camera
	for i in range(n):

		# Init a new empty data array
		arr = np.zeros(shape=(m_sizes[i]))

		# Go through ctypes array and extract data for this camera
		for j in range(m_sizes[i]):
			arr[j] = ptr[i][j]

		# Add the data for this camera to the final list
		arr_list.append(arr)

	return arr_list


#######################################################################################
#################################### PRIMARY CLASS ####################################
#######################################################################################


class GuralSpectral(object):
	"""
	Class containing input parameters, spectral functions, adjustments to the model  
	spectra, calculation of elemental abundances, extinction  factors, 
	and warm/hot column densities.

	Loads spectral library, defines function inputs as C type variables to be 
	passed to Gural's C spectral library.
	"""


	def __init__(self, high_temp, low_temp, elem, extinction_factor, grat_info, \
		col_dense):
		"""
		Defines input parameters and loads library with function definitions 
		and conversion to C-type variables.
		"""

		########################## Initialize input parameters ########################

		self.high_temp = high_temp
		self.low_temp = low_temp
		self.elem = elem
		self.extinction_factor = extinction_factor
		self.grating_info_object = GratingInfo
		self.grat_info = [self.grating_info_object]  # grat_info
		self.col_dense = col_dense

		# Calls load library function with all spectral functions
		self.loadLibrary()

		# Define class names, must be done after loadLibrary
		self.spconfig = SpecConfiguration()
		self.spcalib = SpecCalibration()
		self.elemdata = ElementsData()
		self.spectra = MeasFitSpectra()
		self.starspectra = StarSpectraInfo()

		# Add variables needed by plasmaVolumes - MJM
		self.vinfinity_kmsec = 40.0     # km/sec
		self.approach_angle_radians = 55.0 * 3.141592654 / 180.0
		self.earth_radius_km = 6378.16  # WGS84
		self.site_height_km = 0.2       # Above WGS84 Earth surface
		self.meteor_height_km = 85.0    # Above WGS84 Earth surface
		self.altitude_deg = 45.0        # elevation angle above horizon
		Rsin = self.earth_radius_km * math.sin(self.altitude_deg * 3.141592654 / 180.0)
		self.meteor_range_km = math.sqrt( Rsin * Rsin + 2.0 * self.earth_radius_km * self.meteor_height_km + self.meteor_height_km * self.meteor_height_km) - Rsin

	def loadLibrary(self):

		### Initialize
		# Path to the spectral library
		SPECTRAL_LIBRARY = os.path.join("spectral_library", "SpectralTest")

		# Load the spectral library
		self.spectral_lib = npct.load_library(SPECTRAL_LIBRARY, os.path.dirname(__file__))

		######################## Configuration input function #########################

		self.spectral_lib.ReadSpectralConfigfile.restype = VOID
		self.spectral_lib.ReadSpectralConfigfile.argtypes = [
				ct.POINTER(CHAR),
				ct.POINTER(SpecConfiguration)#self.spconfig)
		]

		##### Spectral calibration file IO functions (responsivity, exctinction) ######

		self.spectral_lib.AllocMemorySpectralCalWavelengths.restype = INT 
		self.spectral_lib.AllocMemorySpectralCalWavelengths.argtypes = [
				ct.POINTER(SpecCalibration),#self.spcalib),
				INT, 
				INT,
				DOUBLE,
				DOUBLE
		]

		self.spectral_lib.ReadSpectralCALfile.restype = INT
		self.spectral_lib.ReadSpectralCALfile.argtypes = [ 
				ct.POINTER(CHAR),
				ct.POINTER(SpecConfiguration),#self.spconfig),
				ct.POINTER(SpecCalibration)#self.spcalib)
		]

		self.spectral_lib.WriteSpectralCALfile.restype = INT
		self.spectral_lib.WriteSpectralCALfile.argtypes = [
				ct.POINTER(CHAR),
				DOUBLE,
				ct.POINTER(SpecCalibration)#self.spcalib)
		]

		self.spectral_lib.FreeMemorySpectralCALfile.restype = VOID 
		self.spectral_lib.FreeMemorySpectralCALfile.argtypes = [ 
				ct.POINTER(SpecCalibration)#self.spcalib)
		] 

		self.spectral_lib.GetMostRecentSPCAL.restype = INT
		self.spectral_lib.GetMostRecentSPCAL.argtypes = [
				ct.POINTER(CHAR),
				ct.POINTER(CHAR),
				DOUBLE,
				ct.POINTER(CHAR)
		]

		##################### General purpose spectral IO functions ###################

		self.spectral_lib.WriteSpectrum.restype = VOID
		self.spectral_lib.WriteSpectrum.argtypes = [ 
				ct.POINTER(CHAR),    # In the C definition this is 'const car *pathname'
				INT,
				PDOUBLE,
				PDOUBLE,
				PDOUBLE
		]

		########################### Star Spectra Functions ############################

		self.spectral_lib.ReadStarSpectra.restype = VOID 
		self.spectral_lib.ReadStarSpectra.argtypes = [ 
				ct.POINTER(CHAR),
				ct.POINTER(StarSpectraInfo) #self.starspectra)
		]

		self.spectral_lib.GetStarIndexFromSPType.restype = INT
		self.spectral_lib.GetStarIndexFromSPType.argtypes = [ 
				ct.POINTER(StarSpectraInfo),#self.starspectra),
				ct.POINTER(CHAR)
		]

		self.spectral_lib.GetStarIndexFromUserID.restype = INT
		self.spectral_lib.GetStarIndexFromUserID.argtypes = [
				ct.POINTER(StarSpectraInfo), #self.starspectra),
				INT
		]

		self.spectral_lib.GetStarIndexFromHIP.restype = LONG
		self.spectral_lib.GetStarIndexFromHIP.argtypes = [ 
				ct.POINTER(StarSpectraInfo), #self.starspectra),
				LONG 
		]

		self.spectral_lib.InterpolateStarSpectrum.restype = VOID 
		self.spectral_lib.InterpolateStarSpectrum.argtypes = [ 
				ct.POINTER(StarSpectraInfo), #self.starspectra),
				LONG,
				INT,
				PDOUBLE,
				PDOUBLE
		]

		self.spectral_lib.FreeMemoryStarSpectra.restype = VOID
		self.spectral_lib.FreeMemoryStarSpectra.argtypes = [ 
				ct.POINTER(StarSpectraInfo), #self.starspectra)
		]

		#################### Responsivity and extinction functions ####################

		self.spectral_lib.ResponsivityExtinctionEstimate.restype = DOUBLE
		self.spectral_lib.ResponsivityExtinctionEstimate.argtypes = [ 
				PDOUBLE,
				PDOUBLE,
				DOUBLE,
				DOUBLE,
				ct.POINTER(SpecCalibration)#self.spcalib)
		] 

		self.spectral_lib.ResponsivityAccumulate.restype = VOID 
		self.spectral_lib.ResponsivityAccumulate.argtypes = [ 
				INT,
				DOUBLE,
				ct.POINTER(SpecCalibration)#self.spcalib)
		]

		self.spectral_lib.ExtinctionAccumulate.restype = VOID
		self.spectral_lib.ExtinctionAccumulate.argtypes = [ 
				INT,
				DOUBLE,
				ct.POINTER(SpecCalibration)#self.spcalib)
		]

		self.spectral_lib.Weighting.restype = DOUBLE
		self.spectral_lib.Weighting.argtypes = [ 
				INT,
				DOUBLE
		]

		#################### Elements and molecules input functions ###################

		self.spectral_lib.ReadElementsData.restype = VOID 
		self.spectral_lib.ReadElementsData.argtypes = [ 
				INT, 
				PDOUBLE,
				ct.POINTER(CHAR),
				ct.POINTER(ElementsData)#self.elemdata),
		]

		self.spectral_lib.FreeElementsMemory.restype = VOID 
		self.spectral_lib.FreeElementsMemory.argtypes = [ 
				ct.POINTER(ElementsData)#self.elemdata)
		]

		self.spectral_lib.AdjustableParametersDefaults.restype = VOID 
		self.spectral_lib.AdjustableParametersDefaults.argtypes = [ 
				ct.POINTER(SpecConfiguration),#self.spconfig),
				ct.POINTER(ElementsData)#self.elemdata)
		]

		################### Support functions for extinction fiting ###################

		self.spectral_lib.ExtinctionApprox.restype = DOUBLE
		self.spectral_lib.ExtinctionApprox.argtypes = [  
				DOUBLE,
				DOUBLE,
				DOUBLE
		]

		self.spectral_lib.Airmass.restype = DOUBLE
		self.spectral_lib.Airmass.argtypes = [ 
				DOUBLE,
				DOUBLE
		]
		
		self.spectral_lib.RhoIntegral.restype = DOUBLE
		self.spectral_lib.RhoIntegral.argtypes = [  
				DOUBLE,
				DOUBLE
		]

		self.spectral_lib.FitExtinctionCoefs.restype = VOID 
		self.spectral_lib.FitExtinctionCoefs.argtypes = [ 
				INT,
				PDOUBLE,
				PDOUBLE,
				ct.POINTER(SpecCalibration)#self.spcalib)
		]

		self.spectral_lib.TauOpticalDepth.restype = DOUBLE
		self.spectral_lib.TauOpticalDepth.argtypes = [ 
				INT,
				DOUBLE
		]

		self.spectral_lib.ExtinctionModel.restype = VOID 
		self.spectral_lib.ExtinctionModel.argtypes = [ 
				INT,
				PDOUBLE,
				PDOUBLE,
				DOUBLE,
				DOUBLE,
				ct.POINTER(SpecCalibration)#self.spcalib)
		]

		self.spectral_lib.AbsorptionLine.restype = DOUBLE
		self.spectral_lib.AbsorptionLine.argtypes = [ 
				DOUBLE,
				DOUBLE,
				DOUBLE,
				DOUBLE
		]

		################# Measured and fit spectral memory allocation #################

		self.spectral_lib.AllocMemorySpectrum.restype = INT
		self.spectral_lib.AllocMemorySpectrum.argtypes = [ 
				INT,
				INT,
				ct.POINTER(MeasFitSpectra)#self.spectra)
		]

		self.spectral_lib.FreeMemorySpectrum.restype = VOID 
		self.spectral_lib.FreeMemorySpectrum.argtypes = [ 
				ct.POINTER(MeasFitSpectra)#self.spectra)
		]

		#################### Electron and column density functions ####################

		self.spectral_lib.PlasmaRadius.restype = DOUBLE
		self.spectral_lib.PlasmaRadius.argtypes = [ 
				DOUBLE
		]

		self.spectral_lib.PlasmaVolumes.restype = VOID 
		self.spectral_lib.PlasmaVolumes.argtypes = [ 
				DOUBLE, 
				DOUBLE,
				DOUBLE,
				DOUBLE,
				ct.POINTER(ElementsData)#self.elemdata)
		]

		self.spectral_lib.GetElementIndex.restype = INT
		self.spectral_lib.GetElementIndex.argtypes = [ 
				INT,
				INT,
				ct.POINTER(ElementsData)#self.elemdata)
		]

		self.spectral_lib.GetPrimaryNeutralElement.restype = INT
		self.spectral_lib.GetPrimaryNeutralElement.argtypes = [ 
				ct.POINTER(ElementsData)#self.elemdata)
		]

		self.spectral_lib.JonesElectronDensity.restype = DOUBLE
		self.spectral_lib.JonesElectronDensity.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				INT
		]

		self.spectral_lib.IterativeElectronDensity.restype = DOUBLE
		self.spectral_lib.IterativeElectronDensity.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				INT,
				DOUBLE
		]

		self.spectral_lib.ColumnDensities_NumberAtoms.restype = VOID  
		self.spectral_lib.ColumnDensities_NumberAtoms.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				DOUBLE
		]

		####################### Abundance calculation functions #######################

		self.spectral_lib.ResetOneElementAbundance.restype = VOID 
		self.spectral_lib.ResetOneElementAbundance.argtypes = [ 
				INT,
				DOUBLE,
				ct.POINTER(ElementsData)#self.elemdata)
		]

		self.spectral_lib.ResetAllElementAbundances.restype = VOID 
		self.spectral_lib.ResetAllElementAbundances.argtypes = [ 
				DOUBLE,
				ct.POINTER(ElementsData)#self.elemdata)
		]

		self.spectral_lib.ComputeRelativeAbundance.restype = VOID 
		self.spectral_lib.ComputeRelativeAbundance.argtypes = [ 
				ct.POINTER(ElementsData)#self.elemdata)
		]

		################## Spectra calculation and support functions ##################

		self.spectral_lib.IronOnlySpectrumT.restype = VOID 
		self.spectral_lib.IronOnlySpectrumT.argtypes = [ 
				DOUBLE,
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE,
				PDOUBLE,
				DOUBLE,
				PDOUBLE
		]

		self.spectral_lib.ActiveElementsSpectraTlo.restype = VOID 
		self.spectral_lib.ActiveElementsSpectraTlo.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE,
				PDOUBLE,
				DOUBLE,
				PDOUBLE
		]

		self.spectral_lib.ActiveElementsSpectraThi.restype = VOID 
		self.spectral_lib.ActiveElementsSpectraThi.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE,
				PDOUBLE,
				DOUBLE,
				PDOUBLE
		]

		self.spectral_lib.FitSpectralCoefficients.restype = VOID 
		self.spectral_lib.FitSpectralCoefficients.argtypes = [ 
				DOUBLE,
				PDOUBLE,
				PDOUBLE,
				ct.POINTER(ElementsData),#self.elemdata),
				ct.POINTER(SpecConfiguration)#self.spconfig)
		]

		self.spectral_lib.SpectrumGivenOnlyLockedCoefs.restype = VOID 
		self.spectral_lib.SpectrumGivenOnlyLockedCoefs.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE
		]

		self.spectral_lib.SpectrumGivenOnlyFittingCoefs.restype = VOID 
		self.spectral_lib.SpectrumGivenOnlyFittingCoefs.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE 
		]

		self.spectral_lib.SpectrumGivenAllCoefs.restype = VOID 
		self.spectral_lib.SpectrumGivenAllCoefs.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE
		]

		self.spectral_lib.SpectrumGivenAllCoefs_Subset.restype = VOID 
		self.spectral_lib.SpectrumGivenAllCoefs_Subset.argtypes = [ 
				ct.POINTER(ElementsData),#self.elemdata),
				PDOUBLE,
				INT,
				INT
		]

		self.spectral_lib.FitFeSpectrum.restype = VOID 
		self.spectral_lib.FitFeSpectrum.argtypes = [ 
				DOUBLE, 
				PDOUBLE, 
				PDOUBLE, 
				ct.POINTER(ElementsData),#self.elemdata),
				ct.POINTER(SpecConfiguration),#self.spconfig),
				PDOUBLE, 
				PDOUBLE
		]

		self.spectral_lib.svd_NxN.restype = VOID 
		self.spectral_lib.svd_NxN.argtypes = [ 
				PDOUBLE, 
				PDOUBLE, 
				PDOUBLE, 
				PDOUBLE, 
				INT, 
				INT
		]

		##################### Julian date and calendar functions ######################

		self.spectral_lib.JulianDateAndTime.restype = DOUBLE
		self.spectral_lib.JulianDateAndTime.argtypes = [ 
				LONG,
				LONG,
				LONG,
				LONG,
				LONG,
				LONG,
				LONG
		]

		self.spectral_lib.CalendarDateAndTime.restype = VOID 
		self.spectral_lib.CalendarDateAndTime.argtypes = [ 
				DOUBLE,
				ct.POINTER(LONG),
				ct.POINTER(LONG),
				ct.POINTER(LONG),
				ct.POINTER(LONG),
				ct.POINTER(LONG),
				ct.POINTER(LONG),
				ct.POINTER(LONG)
		]

		### Initialize camera and grating parameters ###

		self.camos_camera_index = 0

		self.grat_info[self.camos_camera_index].camnum = 101

	
	############################# Library functions ###################################

	def writeSpectrum(self, index):
		# self.elemdata.els[46].spechi = 10000
		self.spectral_lib.WriteSpectrum(b"test.txt", self.spcalib.nwavelengths, self.spcalib.wavelength_nm, self.elemdata.els[index].speclo, self.elemdata.els[index].spechi)

	def writeFullSpectrum(self):
		self.spectral_lib.WriteSpectrum(b"test.txt", self.spcalib.nwavelengths, self.spcalib.wavelength_nm, self.spectra.fit_spectrum, self.spectra.fit_spectrum)

	def writeFullSpectrum2(self, out_file):
		with open(out_file, 'w') as f:
			### Write file header
			# file_id / name
			# 
			#
			#
			for w in range(self.spcalib.nwavelengths):
				f.write('%f %f %f\n' % (w, self.spcalib.wavelength_nm[w], self.spectra.fit_spectrum[w]))

	def loadMeteorTrajectory(self, traj_path):
		""" 
		Initalizes trajectory values to be taken from pickle file.
		Values will be reassigned from detailed trajectory file.
		Function must be re-run for every new frame.

		Arguments:
			traj_path: path to the trajectory pickle file for given event
		Return:
			self.vininity_kmsec: incoming velocity in m/s
			self.approach_angle_radians: approach angle in radians
			self.site_height_km: elevation of observing station above Earth's surface
			self.meteor_height_km: height of meteor above Earth's surface
			self.altitude_deg: deg above horizon (elevation angle) of meteor
			self._meteor_range_km: meteor range in km
			self.grat_info: grating info for the camera
		"""

		# Incoming velocity, m/s
		self.vinfinity_kmsec = 40 

		# Approach angle, in radians
		self.approach_angle_radians = np.radians(55) 

		# Site height in km, above the Earth's surface (WGS84)
		self.site_height_km = 0.2 

		# Meteor height in km, above Earth's surface (WGs84)
		self.meteor_height_km = 85

		# Elevation angle above horizon
		self.altitude_deg = 45 

		# Computational value to get meteor range in km
		self.Rsin = EARTH_RADIUS_KM * np.sin(np.radians(self.altitude_deg))

		# Meteor range, calculated with Rsin, Earth's radius, meteor height, etc.
		self.meteor_range_km = np.sqrt(self.Rsin * self.Rsin + 2.0 * EARTH_RADIUS_KM\
			 * self.meteor_height_km + self.meteor_height_km * self.meteor_height_km)\
				  - self.Rsin

		# Grating information, calculated value 
		self.grat_info[self.camos_camera_index].grating_area_scale = np.cos(np.radians\
			(self.spconfig.grating_offnormal_deg)) 
	
		return self.vinfinity_kmsec, self.approach_angle_radians, self.site_height_km,\
			self.meteor_height_km, self.altitude_deg, self.meteor_range_km, \
				self.grat_info


	def getElementIndex(self, ionization_state, atomic_mass_num):
		"""
		Retrieves an element from the elemdata class.
		
		Arguments:
			atomic_mass_num: atomic mass number of desired element
		Return:
			element_index: information for desired element
		"""

		# Extract element based on atomic number
		element_index = self.spectral_lib.GetElementIndex(atomic_mass_num, ionization_state, self.elemdata)
		# print(element_index)

		# Returns element
		return element_index

	
	def plasmaVolumes(self):
		"""
		If any of the arguments in plasmaVolumes changes at a later time, plasmaVolumes must
		be called again to infill the elemedata structure values for height, range, and volumes. 
		For example, if a different framee is being processed, this function must be called again.  

		Return:
			Plasma volumes 
		"""

		# Compute plasma volumes
		# self.vinfinity_kmsec = vinf
		# print('Vinf: %f' % self.vinfinity_kmsec)
		self.spectral_lib.PlasmaVolumes(self.meteor_height_km, self.meteor_range_km, self.approach_angle_radians, \
			self.elemdata.hot2warm, self.elemdata)


	
	def extinctionModel(self):
		"""
		Compute the model for extinction and the airmass given the event's metadata. 
		Update for each new frame analyzed.altitude_deg

		Return:
			Extinction model
		"""

		# Compute extinction model
		self.spectral_lib.ExtinctionModel(self.spcalib.nwavelengths, self.spcalib.wavelength_nm, \
			self.spcalib.modl_extn_spec, self.elemdata.Xfactor * self.altitude_deg, \
				self.earth_radius_km + self.site_height_km, self.spcalib)

	
	def resetAllElementalAbundances(self):
		"""
		Zeros all the elemental  abundances and number of atoms.
		Set all  element fitting flags to not selected for fitting = FITLESS.
		Compute Jones 1997 fraction of ionized atoms
			Beta = n+ / (n+ + no) = function( Vinf )
		"""

		# Reset all elemental abundances
		self.spectral_lib.ResetAllElementAbundances(self.vinfinity_kmsec, self.elemdata)


	def computeWarmPlasmaSpectrum(self):
		"""
		Computes spectrum for warm plasma. 
		Must be re-called every time something is changed.

		Return:
			Model spectra for warm plasma
		"""

		# Calls Gural's ActiveElementsSPectraTlo function to get spectra for warm plasma
		self.spectral_lib.ActiveElementsSpectraTlo(self.elemdata,
								self.spcalib.ord1_resp_spec, 
								self.spcalib.ord2_resp_spec, 
								self.spcalib.gratinfo[self.camos_camera_index].grating_area_scale,
								self.spcalib.modl_extn_spec)	                      


	def computeHotPlasmaSpectrum(self):
		"""
		Computes spectrum for hot plasma.
		Must be re-called every time something is changed. 
		
		Return:
			Model spectra for hot plasma
		"""

		# Calls Gural's ActElementsSpectraThi function to get spectra for hot plasma
		self.spectral_lib.ActiveElementsSpectraThi(self.elemdata,
								self.spcalib.ord1_resp_spec, 
								self.spcalib.ord2_resp_spec, 
								self.spcalib.gratinfo[self.camos_camera_index].grating_area_scale,
								self.spcalib.modl_extn_spec)

	
	def changeHot2WarmRatio(self, hot2warm):
		"""
		Changes ratio between hot and warm plasma. 
		
		Arguments: 
			hot2warm: ratio of hot to warm plasma
		Return:
			updated model spectra for hot plasma, with new hot to warm ratio
		"""

		# Set ratio of hot to warm plasma
		self.elemdata.hot2warm = hot2warm 

		# Adjust plasma volume,  must  be called every time  hot2warm is changed
		self.plasmaVolumes()
		
		# Update the hot plasma spectrum
		self.computeHotPlasmaSpectrum()


	def changeBroadening(self, sigma0):
		"""
		Obtain a new model if broadening changes.

		Arguments:
			sigma0: changed broadening factor
		Return:
			updated  model spectra for hot and warm plasma
		"""

		# Set sigma0
		self.sigma0 = sigma0

		# Update warm plasma spectrum
		self.computeWarmPlasmaSpectrum()

		# Update hot plasma spectrum
		self.computeHotPlasmaSpectrum()
	
	
	def elemFitting(self, elem):
		"""
		Defines which element is being fitted

		Arguments:
			elem: element to be fitted
		"""

		self.elemdata.els[elem].user_fitflag = FITTING
	

	def removeElemFromModel(self, elem):
		"""
		Remove an element from the model spectra
		
		Arguments:
			elem: element to be turned off
		Return: 
			updated model spectra for warm plasma
			updated model spectra for hot plasma
		"""

		# Remove the selected element
		self.elemdata.els[elem].user_fitflag = FITLESS

		# Update warm and hot temperatures
		self.elemdata.Tlo = self.spconfig.nominal_lo_exc_temp
		self.elemdata.Thi = self.spconfig.nominal_hi_exc_temp

		# Update hot and warm spectra
		self.computeWarmPlasmaSpectrum
		self.computeHotPlasmaSpectrum


	def lockElemFit(self, elem):
		"""
		Lock an element's abundance (column density?) in the model spectra
		
		Arguments:
			elem: element to be locked
		Return:
			updated model spectra for warm plasma
			updated model spectra for hot plasma
		"""

		# Lock the selected element
		self.elemdata.els[elem].user_fitlag = FITLOCK

		# Update warm and hot temperatures
		self.elemdata.Tlo = self.spconfig.nominal_lo_exc_temp
		self.elemdata.Thi = self.spconfig.nominal_hi_exc_temp

		# Update hot and warm spectra
		self.computeWarmPlasmaSpectrum
		self.computeHotPlasmaSpectrum

	
	def setReferenceElem(self, elem = None):
		"""
		Sets the reference element for relative abundances. Can either
		provide reference elements in order of Mg, Na, Fe, or the user can select
		the element they'd like to use.
		
		Arguments:
			elem: (optional) element to be used as reference, should be neutral. If none
			is selected, the function will default to Mg, then Na, then Fe
		Return:
			kelem_ref: reference element
		"""
		if elem is None:

			# Selects first active fitting element in the order of Mg, Na, Fe
			self.elemdata.kelem_ref = self.spectral_lib.GetPrimaryNeutralElement(self.elemdata)
			# self.elemdata.kelem_ref = self.kelem_ref
			# print(self.elemdata.kelem_ref)
			self.kelem_ref = self.elemdata.kelem_ref

		else: 

			self.elemdata.kelem_ref = elem  
			# self.elemdata.kelem_ref = self.kelem_ref
			self.kelem_ref = self.elemdata.kelem_ref
		
		return self.kelem_ref
	
	
	def computeRelativeAbundances(self, kelem_ref):
		"""
		Compute the abundances relative to a reference element.
		"""

		# Compute relative abundance
		self.spectral_lib.ComputeRelativeAbundance(self.elemdata)

		for kneutral in range(0, self.elemdata.nneutrals):
			kelem = self.elemdata.neutral_index[kneutral]
			if (self.elemdata.els[kelem].user_fitflag == FITLESS):
				continue #... do not report non-active elements
	
	
	def computeFullSpec(self, elem1, elem2, elem3=None):
		"""
		Compute the full spectrum which requires column densities for warm, hot, neutrals, ions
		Artificially set the neutral warm column densities of the 2 elements
		Get a first guess at the electron density from Jones
		Iteratively refine the electron density via Jenniskens method
		Get the remaining column densities for warm-ions, hot-neutrals, hot-ions
		Compute spectrum for all active elements.
		
		Arguments:  
			elem1: first element for which a fit is desired
			elem2: second element for which a fit is desired
		Return:
			A model spectrum calculated with coefficients determined by the function
		"""

		# Normally get col density from a fit (see later)
		self.elemdata.els[elem1].N_warm = 3.0e+9 # Check this...
		# print("** %s" % self.elemdata.els[elem1].N_warm)

		# Normally get col density from a fit (see later)
		self.elemdata.els[elem2].N_warm = 1.2e+09 # Check this... MJM
		# print("** %s" % self.elemdata.els[elem2].N_warm)

		if elem3 != None:
			# self.elemdata.els[elem3].N_warm = 1.0e+09
			self.elemdata.els[elem3].N_warm = 3.0e+10
			# print("Elem3 n_warm: %s" % self.elemdata.els[elem3].N_warm)

		# Sets ne_jones
		ne_guess = self.spectral_lib.JonesElectronDensity(self.elemdata, self.elemdata.kelem_ref)
		print("+++++ ne_Jones %s" % ne_guess)         

		# Sets ne_iter
		self.ne = float(self.spectral_lib.IterativeElectronDensity(self.elemdata, self.elemdata.kelem_ref, ne_guess))
		print("+++++ ne_iter %s" % self.ne)  
		
		# Caclulate column densities
		self.spectral_lib.ColumnDensities_NumberAtoms(self.elemdata, self.ne)
		print("Column density #1 = %s" % self.elemdata.els[elem1].N_warm)
		print("Column density #2 = %s" % self.elemdata.els[elem2].N_warm)

		# Calculate the  spectrum, given all coefficients previously calculated
		print('ComputeFullSpec')
		self.spectral_lib.SpectrumGivenAllCoefs(self.elemdata, self.spectra.fit_spectrum)
		print('ComputeFullSpec')

	
	def fitMeasSpec(self, elem1, elem2, elem3=None):
		"""
		Fit the measured spectrum to the active elements.
		Note that the active elements were previously set to fitting and
		the active element's model spectra have been calculated.
		First computes the coefficients for the warm column densities.
		Next computes the remaining column densities.
		
		Return:
			Model spectrum calculated with coefficients determined by the function
		"""
		print("--- fitMeasSpec ---")
		# Normally get col density from a fit (see later)
		self.elemdata.els[elem1].N_warm = 3.0e+9 # Check this...
		print("Num warm element 1: %s" % self.elemdata.els[elem1].N_warm)

		# Normally get col density from a fit (see later)
		self.elemdata.els[elem2].N_warm = 1.2e+09 # Check this... MJM
		print("Num warm element 2: %s" % self.elemdata.els[elem2].N_warm)

		if elem3 != None:
			# self.elemdata.els[elem3].N_warm = 1.0e+09
			self.elemdata.els[elem3].N_warm = 3.0e+10

		# Sets ne_jones
		ne_guess = self.spectral_lib.JonesElectronDensity(self.elemdata, self.elemdata.kelem_ref)   
		print('Calculated ne_guess: %s' % ne_guess)      

		# Sets ne_iter
		self.ne = float(self.spectral_lib.IterativeElectronDensity(self.elemdata, \
			self.elemdata.kelem_ref, ne_guess))
		print('Calculated ne_iter: %e' % self.ne)

		# Fit spectral coefficients
		print('--- Fitting spectral coefficients...')
		self.spectral_lib.FitSpectralCoefficients(self.spcalib.nwavelengths, \
			self.spcalib.wavelength_nm, self.spectra.meas_spectrum, self.elemdata, self.spconfig)

		# Calculate column densities
		print('--- Calculating column densities...')
		self.spectral_lib.ColumnDensities_NumberAtoms(self.elemdata, self.ne)

		# Calculate model spectrum
		print('--- Calculating model spectrum...')
		self.spectral_lib.SpectrumGivenAllCoefs(self.elemdata, self.spectra.fit_spectrum)
		print('FitMeasuredSpec')
	
		print("Num warm element 1: %s" % self.elemdata.els[elem1].N_warm)
		print("Num warm element 2: %s" % self.elemdata.els[elem2].N_warm)
	
	def scaleWarmColumnDensity(self, elem, scale_factor):
		"""
		Scales up warm column density for a particular element by a given factor
		Then updates the spectrum
		
		Arguments:
			elem: element for which the user wants to scale up the warm column density
			scale_factor: the factor by which the element should be scaled up
		Return:
			Model spectrum calculated with coefficients determined by the function
		"""

		# Adujust but the scale factor
		self.elemdata.els[elem].N_warm *= scale_factor

		# Calculate column densities
		self.spectral_lib.ColumnDensities_NumberAtoms(self.elemdata, self.ne)

		# Calculate model spectra
		self.spectral_lib.SpectrumGivenAllCoefs(self.elemdata, self.spectra.fit_spectrum)

	
	def fitElem1LockElem2(self, elem1, elem2):
		"""
		Spectral fit with elem1 FITTED and elem2 FITLOCK
		Column density for elem1 is fitted, column density for elem2 is locked
		
		Arguments:
			elem1: element that is to be FITTED
			elem2: element that is to be LOCKED at a constant value
		Return:
			A model spectra calculated with the coefficients determined by the function

		"""

		# Set elem1 to fitting and elem2 to fitlock
		self.elemdata.els[elem1].user_fitflag = FITTING; #... fits this element's Nwarm col dens
		self.elemdata.els[elem2].user_fitflag = FITLOCK; #... no change to Nwarm col dens

		# Fit relevant element
		self.spectral_lib.FitSpectralCoefficients(self.spcalib.nwavelengths, \
			self.spcalib.wavelength_nm, self.spectra.meas_spectrum, self.elemdata, self.spconfig)

		# Get column densities
		self.spectral_lib.ColumnDensities_NumberAtoms(self.elemdata, self.ne)

		# Calculate the model spectrum
		self.spectral_lib.SpectrumGivenAllCoefs(self.elemdata, self.spectra.fit_spectrum)

	def getStarIndexFromHIP(self, hip):
		# Define hipparcos number
		self.hip = hip

		# Define star index
		# self.star_index = self.spectral_lib.GetStarIndexFromHIP(self.starspectra, self.hip)
		# Interpolate star spectrum
		self.spectral_lib.InterpolateStarSpectrum(self.starspectra, self.star_index, self.spcalib.nwavelengths, \
			self.spcalib.wavelength_nm, self.spcalib.ref_star_spec)

		# print('Yay')


	def interpolateCatalogedSpectrum(self, hipno):
		"""
		Get index for the desired star from the star catalog, show the function
		call to interpolate its catalogued spectrum to the user defined wavelengths. 

		Arguments:
			hip: Hipparcos number
		Return:
			Interpolated star spectrum
		"""

		# print('Testing 4 5 6')

		# Define hipparcos number
		self.hip = hipno
		# print(self.hip)

		# Define star index
		self.star_index = self.spectral_lib.GetStarIndexFromHIP(self.starspectra, self.hip)
		# print(self.star_index)

		# Interpolate star spectrum
		self.spectral_lib.InterpolateStarSpectrum(self.starspectra, self.star_index, self.spcalib.nwavelengths, \
			self.spcalib.wavelength_nm, self.spcalib.ref_star_spec)
		
		# Apply delay 
		self.spectral_lib.Delay_msec(1000)


	def clearRespExtinctVectors(self):
		"""
		Clear responsivity and extinction accumulator vectors.
		"""
		
		# Clear responsivity vector
		self.spectral_lib.ResponsivityAccumulate(RECLEAR, 1.0, self.spcalib)

		# Clear extinction vector
		self.spectral_lib.ExtinctionAccumulate(RECLEAR, 1.0, self.spcalib)


	def updateRespExtinct(self, counter):
		"""
		Loop over a certain number of stars to update the responsivity and extinction
		for the night. Find the star's index by its name. 

		Arguments:
			counter: number of stars to count
		Return:
			Updated responsivity and extinction
		"""

		self.starcounter = counter

		while self.starcounter < 3:

			star_index = self.starcounter  #... Normally get this index from the star name collected


			#-------- Interpolate this cataloged star's spectrum and place in the spcalib structural
			#         element ref_star_spec at the user specified wavelengths via linear interpolation.

			self.spectral_lib.InterpolateStarSpectrum(self.starspectra, star_index, self.spcalib.nwavelengths, self.spcalib.wavelength_nm, self.spcalib.ref_star_spec)

			#-------- Generate an "artificial" star measurement by convolving the star with the 
			#         current 1st order responsivity, latest extinction model and added random noise,
			#         at the user specified wavelengths.
			#         NOTE: Normally you would extract the star spectrum from a star collection.

			for kwave in range(0, self.spcalib.nwavelengths):
				noise_component = self.noise_multiplier * self.starspectra.star[star_index].specmax * (2.0 * np.random.rand() - 1.0)
				self.spectra.meas_spectrum[kwave] = self.spcalib.ref_star_spec[kwave] * self.spcalib.ord1_resp_spec[kwave] * self.spcalib.modl_extn_spec[kwave] + noise_component

			#-------- Compute the individual star responsivity and extinction estimate and let user decide
			#         to add it to the cummulative running results or skip to next star.

			Xairmass = self.spectral_lib.ResponsivityExtinctionEstimate(self.spectra.meas_spectrum, self.spcalib.ref_star_spec, self.altitude_deg, self.site_height_km, self.spcalib)


			#-------- Compute the weighting function

			weight = self.spectral_lib.Weighting( UNITYWEIGHTS, 0.0)  #... uses unity weighting


			#.... Display the estimated responsivity and extinction over the ord1_resp and modl_extn

				#**** Plot the display of responsivity and extinction ****


			#.... IF airmass below the configuration file limit and ...
			#     IF the user decides to include the estimated responsivity to the daily running response

			if(Xairmass <= self.spconfig.airmass_limit):
				self.spectral_lib.ResponsivityAccumulate( REDAILY, weight, self.spcalib )


			#.... IF the user decides to include the estimated extinction to the daily running extinction tally

			self.spectral_lib.ExtinctionAccumulate(REDAILY, weight, self.spcalib)


			#-------- Go to next star spectrum collected

			self.starcounter += 1

			#... end of star collection loop for the night


	def displayRespExtinctEstimate(self, include = False, save = True):
		"""
		Display the estimated responsivity and extinction (ord1_resp, modl_resp).
		User can decide to save the updated results or pass on writing/updating SPCAL file
		for this data if it was of poor quality. 
		SPCAL file has date/time in the filename to make it unique, so it will not
		overwrite previous SPCAL.

		Arguments:
			include: default = False. Set to True if the user decides to include
				daily responsivity to the multi-day cummulative response.
			save: defautl = True. Set to False if the user does not want to save
				the results of the updated SPCAL for this data. 
		"""

		# IF the user decides to include the daily responsivity to the multi-day cummulative response
		if include == True:
			self.spectral_lib.ResponsivityAccumulate(REMULTI, 1.0, self.spcalib)
		
		# IF the user decides to include the daily extinction to the multi-day cummulative extinction
		else:
			self.spectral_lib.ExtinctionAccumulate( REMULTI, 1.0, self.spcalib)

		# Fit the extinction model to the 5 known contributors to update the extinction coefficients
		self.spectral_lib.FitExtinctionCoefs(self.spcalib.nwavelengths, self.spcalib.wavelength_nm, \
			self.spcalib.cumm_extn_spec, self.spcalib)

		# Calculate Julian date and time
		jdt = self.spectral_lib.JulianDateAndTime(2013, 4, 30, 8, 0, 0, 0)

		# Write updated responsivity and extinction to a file
		if save == True:
			self.spectral_lib.WriteSpectralCALfile(" spectral_library/SPCAL/".encode('ascii'), self.jdt, self.spcalib)

		# Calculate model fit extinction
		self.spectral_lib.ExtinctionModel(self.spcalib.nwavelengths, self.spcalib.wavelength_nm, \
			self.spcalib.modl_extn_spec, self.elemdata.Xfactor * self.altitude_deg, \
				self.earth_radius_km + self.site_height_km, self.spcalib)


	def clearAllMemory(self):
		"""
		Free all allocated memory.
		"""
		
		# Free memory associated with spectra
		self.spectral_lib.FreeMemorySpectrum(self.spectra)

		# Free memory associated with element data
		self.spectral_lib.FreeElementsMemory(self.elemdata)

		# Free memory associated with spcalib
		self.spectral_lib.FreeMemorySpectralCALfile(self.spcalib)

		# Free memory associated with starspectra
		self.spectral_lib.FreeMemoryStarSpectra(self.starspectra)


	"""
	def run(self):
		Function that runs the simulation

		... 


		return self 
	"""