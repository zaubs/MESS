"""

GUI for Camo-S Analysis

Naming conventions:
Class: ClassName
Functions: functionName
Variable names: variable, variable_name
Widgets: Name_widgetype, i.e. Upload_button, Direct_linedit, Spectral_label

"""
###################################################################################################
################################## IMPORTS AND SUPPORT FUNCTIONS ##################################
###################################################################################################

#################### STANDARD IMPORTS ####################

import numpy as np
import pyqtgraph as pg
import sys 
sys.path.append('../RMS/RMS/Routines')
import matplotlib.pyplot as plt
import scipy.ndimage
import os
import sys
import imageio 
imread = imageio.imread
import math
import bz2
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import (r2_score, mean_absolute_error)
from scipy.signal import savgol_filter
from scipy import interpolate

#################### PYQT5 LIBRARY/PACKAGE IMPORTS ####################

from PyQt5 import QtWidgets, uic, QtGui, QtCore
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from calibration import Ui_CalibrationDialog
from stellarcalibration import Ui_StellarCalibrationDialog
from columndensity import Ui_ColumnDensityDialog

#################### FROM WMPL ####################

from wmpl.Formats.Vid import readVid
from wmpl.Utils.TrajConversions import unixTime2Date
from wmpl.Formats.Plates import loadScale, plateScaleMap

# Cython init
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from BinImageCy import binImage as binImageCy

spectral_library = __import__("CAMO-Spectral_Library")

#################### SUPPORT FUNCTIONS AND CLASSES ####################
def twoDGaussian(params, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y, saturation = params
    if isinstance(saturation, np.ndarray):
        saturation = saturation[0,0]

    xo = float(xo)
    yo = float(yo)

    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp(-(a*((x - xo)**2) + 2*b*(x - xo)*(y - yo) + c*((y-yo)**2)))

    g[g > saturation] = saturation

    return g.ravel()

def fitPSF(imarray, avepixel_mean, x2, y2):
    segment_radius = 25
    roundness_threshold = 0.5
    max_feature_ratio = 0.8

    x_fitted = []
    y_fitted = []
    amplitude_fitted = []
    intensity_fitted = []
    sigma_y_fitted = []
    sigma_x_fitted = []

    # Set the initial guess
    initial_guess = (10.0, segment_radius, segment_radius, 1.0, 1.0, 0.0, avepixel_mean)

    for star in zip(list(y2), list(x2)):
        hot_pixel = False
        small_star = False
        zero_intensity = False
        fitting_failed = False
        max_ratio = False

        y, x = star

        y_min = y - segment_radius
        y_max = y + segment_radius
        x_min = x - segment_radius
        x_max = x + segment_radius

        if y_min < 0:
            y_min = np.array([0])
        if y_max > np.shape(imarray)[0]:
            y_max = np.array(np.shape(imarray)[0])
        if x_min < 0:
            x_min = np.array([0])
        if x_max > np.shape(imarray)[1]:
            x_max = np.array(np.shape(imarray[1]))

        x_min = int(x_min)
        x_max = int(x_max)
        y_min = int(y_min)
        y_max = int(y_max)

        star_seg = imarray[y_min:y_max,x_min:x_max]

        y_ind, x_ind = np.indices(star_seg.shape)

        saturation = (2**(8*star_seg.itemsize) - 1)*np.ones_like(y_ind)

        try:
            popt, pcov = scipy.optimize.curve_fit(twoDGaussian, (y_ind, x_ind, saturation), star_seg.ravel(), \
                p0=initial_guess, maxfev=200)
        except RuntimeError:
            fitting_failed = True

        if fitting_failed == False:
            amplitude, yo, xo, sigma_y, sigma_x, theta, offset = popt
        else:
            amplitude, yo, xo, sigma_y, sigma_x, theta, offset = (0,0,0,1.0,1.0,0,0)

        if min(sigma_y/sigma_x, sigma_x/sigma_y) < roundness_threshold:
            hot_pixel = True

        if (4*sigma_x*sigma_y/segment_radius**2 > max_feature_ratio):
            max_ratio = True

        crop_y_min = int(yo - 3*sigma_y) + 1
        if crop_y_min < 0: crop_y_min = 70

        crop_y_max = int(yo + 3*sigma_y) + 1
        if crop_y_max >= star_seg.shape[0]: crop_y_max = star_seg.shape[0] - 1

        crop_x_min = int(xo - 3*sigma_x) + 1
        if crop_x_min < 0: crop_x_min = 0

        crop_x_max = int(xo + 3*sigma_x) + 1
        if crop_x_max >= star_seg.shape[1]: crop_x_max = star_seg.shape[1] - 1


        if (y_max - y_min) < 3:
            crop_y_min = int(yo - 2)
            crop_y_max = int(yo + 2)
        if (x_max - x_min) < 3:
            crop_x_min = int(xo - 2)
            crop_x_max = int(xo + 2)

        star_seg_crop = star_seg[crop_y_min:crop_y_max,crop_x_min:crop_x_max]

        if (star_seg_crop.shape[0] == 0) or (star_seg_crop.shape[1] == 0):
            small_star = True

        bg_corrected = offset
        intensity = np.sum(star_seg_crop - bg_corrected)

        if intensity <= 0:
            zero_intensity = True

        if (hot_pixel == True) or (small_star == True) or (zero_intensity == True) or (max_ratio == True):
            x_fitted.append(x)
            y_fitted.append(y)
            amplitude_fitted.append(-999)
            intensity_fitted.append(0)
            sigma_y_fitted.append(-999)
            sigma_x_fitted.append(-999)
        else:
            x_fitted.append(x_min + xo)
            y_fitted.append(y_min + yo)
            amplitude_fitted.append(amplitude)
            intensity_fitted.append(intensity)
            sigma_y_fitted.append(sigma_y)
            sigma_x_fitted.append(sigma_x)

    return x_fitted, y_fitted, amplitude_fitted, intensity_fitted, sigma_y_fitted, sigma_x_fitted

# Allows user to adjust image levels 
def adjustLevels(img_array, minv, gamma, maxv, nbits=None, scaleto8bits=False):
    """
     Adjusts levels on image with given parameters.

    Arguments:
        img_array: [ndarray] Input image array.
        minv: [int] Minimum level.
        gamma: [float] gamma value
        Mmaxv: [int] maximum level.

    Keyword arguments:
        nbits: [int] Image bit depth.
        scaleto8bits: [bool] If True, the maximum value will be scaled to 255 and the image will be converted
            to 8 bits.

    Return:
        [ndarray] Image with adjusted levels.
    """
    
    if nbits is None:
       
        # Get the bit depth from the image type
        nbits = 8*img_array.itemsize

    input_type = img_array.dtype

    # Calculate maximum image level
    max_lvl = 2**nbits - 1.0

    # Limit the maximum level
    if maxv > max_lvl:
        maxv = max_lvl

    # Check that the image adjustment values are in fact given
    if (minv is None) or (gamma is None) or (maxv is None):
        return img_array

    minv = minv/max_lvl
    maxv = maxv/max_lvl
    interval = maxv - minv
    invgamma = 1.0/gamma

    # Make sure the interval is at least 10 levels of difference
    if interval*max_lvl < 10:

        minv *= 0.9
        maxv *= 1.1

        interval = maxv - minv
       
    # Make sure the minimum and maximum levels are in the correct range
    if minv < 0:
        minv = 0

    if maxv*max_lvl > max_lvl:
        maxv = 1.0
   
    img_array = img_array.astype(np.float64)

    # Reduce array to 0-1 values
    img_array = np.divide(img_array, max_lvl)

    # Calculate new levels
    img_array = np.divide((img_array - minv), interval)

    # Cut values lower than 0
    img_array[img_array < 0] = 0

    img_array = np.power(img_array, invgamma)

    img_array = np.multiply(img_array, max_lvl)

    # Convert back to 0-maxval values
    img_array = np.clip(img_array, 0, max_lvl)

    # Scale the image to 8 bits so the maximum value is set to 255
    if scaleto8bits:
        img_array *= 255.0/np.max(img_array)
        img_array = img_array.astype(np.uint8)

    else:

        # Convert the image back to input type
        img_array = img_array.astype(input_type)

    return img_array

# Loads the image into the script
def loadImage(img_path, flatten=-1):
    """ Load the given image. Handle loading it using different libraries. 
    
    Arguments:
        img_path: [str] Path to the image.
    Keyword arguments:
        flatten: [int] Convert color image to grayscale if -1. -1 by default.
    """

    img = imageio.imread(img_path, as_gray=bool(flatten))

    return img

# Bins a given image, provides a numpy array 
def binImage(img, bin_factor, method='avg'):
    """ Bin the given image. The binning has to be a factor of 2, e.g. 2, 4, 8, etc.
    This is just a wrapper function for a cythonized function that does the binning.
    
    Arguments:
        img: [ndarray] Numpy array representing an image.
        bin_factor: [int] The binning factor. Has to be a factor of 2 (e.g. 2, 4, 8).
    Keyword arguments:
        method: [str] Binning method.  'avg' by default.
            - 'sum' will sum all values in the binning window and assign it to the new pixel.
            - 'avg' will take the average.
    Return:
        out_img: [ndarray] Binned image.
    """

    input_type = img.dtype

    # Make sure the input image is of the correct type
    if img.dtype != np.uint16:
        img = img.astype(np.uint16)
    
    # Perform the binning
    img = binImageCy(img, bin_factor, method=method)

    # Convert the image back to the input type
    img = img.astype(input_type)

    return img

# Structure containing flat field
class FlatStruct(object):
    def __init__(self, flat_img, dark=None):
        """ Structure containing the flat field.
        Arguments:
            flat_img: [ndarray] Flat field.
        """

        # Convert the flat to float64
        self.flat_img = flat_img.astype(np.float64)

        # Store the original flat
        self.flat_img_raw = np.copy(self.flat_img)

        # Apply the dark, if given
        self.applyDark(dark)

        # Compute the flat median
        self.computeAverage()

        # Fix values close to 0
        self.fixValues()


    def applyDark(self, dark):
        """ Apply a dark to the flat. """

        # Apply a dark frame to the flat, if given
        if dark is not None:
            self.flat_img = applyDark(self.flat_img_raw, dark)
            self.dark_applied = True

        else:
            self.flat_img = np.copy(self.flat_img_raw)
            self.dark_applied = False

        # Compute flat median
        self.computeAverage()

        # Fix values close to 0
        self.fixValues()


    def computeAverage(self):
        """ Compute the reference level. """


        # # Bin the flat by a factor of 4 using the average method
        # flat_binned = binImage(self.flat_img, 4, method='avg')

        # # Take the maximum average level of pixels that are in a square of 1/4*height from the centre
        # radius = flat_binned.shape[0]//4
        # img_h_half = flat_binned.shape[0]//2
        # img_w_half = flat_binned.shape[1]//2
        # self.flat_avg = np.max(flat_binned[img_h_half-radius:img_h_half+radius, \
        #     img_w_half-radius:img_w_half+radius])

        self.flat_avg = np.median(self.flat_img)

        # Make sure the self.flat_avg value is relatively high
        if self.flat_avg < 1:
            self.flat_avg = 1


    def fixValues(self):
        """ Handle values close to 0 on flats. """

        # Make sure there are no values close to 0, as images are divided by flats
        self.flat_img[(self.flat_img < self.flat_avg/10) | (self.flat_img < 10)] = self.flat_avg


    def binFlat(self, binning_factor, binning_method):
        """ Bin the flat. """

        # Bin the processed flat
        self.flat_img = binImage(self.flat_img, binning_factor, binning_method)

        # Bin the raw flat image
        self.flat_img_raw = binImage(self.flat_img_raw, binning_factor, binning_method)

# Load flat
# def loadFlat(dir_path, file_name, dtype=None, byteswap=True, dark=None):
#     """ Load the flat field image. 
#     Arguments:
#         dir_path: [str] Directory where the flat image is.
#         file_name: [str] Name of the flat field file.
#     Keyword arguments:
#         dtype: [bool] A given file type fill be force if given (e.g. np.uint16).
#         byteswap: [bool] Byteswap the flat image. False by default.
#     Return:
#         flat_struct: [Flat struct] Structure containing the flat field info.
#     """

#     # Load the flat image
#     flat_img = loadImage(os.path.join(dir_path, file_name), -1)

#     # Change the file type if given
#     if dtype is not None:
#         flat_img = flat_img.astype(dtype)

#     # If the flat isn't a 8 bit integer, convert it to uint16
#     elif flat_img.dtype != np.uint8:
#         flat_img = flat_img.astype(np.uint16)


#     if byteswap:
#         flat_img = flat_img.byteswap()
        

#     # Init a new Flat structure
#     flat_struct = FlatStruct(flat_img, dark=dark)

#     return flat_struct

# def uploadFlat(self, dtype=None, byteswap=True, dark=None):


            

# Apply flat
def applyFlat(img, flat_struct):
    """ Apply a flat field to the image.
    Arguments:
        img: [ndarray] Image to flat field.
        flat_struct: [Flat struct] Structure containing the flat field.
        
    Return:
        [ndarray] Flat corrected image.
    """

    # Check that the input image and the flat have the same dimensions, otherwise do not apply it
    if img.shape != flat_struct.flat_img.shape:
        return img

    input_type = img.dtype

    # Apply the flat
    img = flat_struct.flat_avg*img.astype(np.float64)/flat_struct.flat_img

    # Limit the image values to image type range
    dtype_info = np.iinfo(input_type)
    img = np.clip(img, dtype_info.min, dtype_info.max)

    # Make sure the output array is the same as the input type
    img = img.astype(input_type)

    return img

# Get handle positions in ROI
def getHandlePositions(self):
    """Return the positions of all handles in local coordinates."""
    pos = [self.mapFromScene(self.lines[0].getHandles()[0].scenePos())]
    for l in self.lines:
        pos.append(self.mapFromScene(l.getHandles()[1].scenePos()))
    return pos

###################################################################################################
########################################## PRIMARY CLASS ##########################################
###################################################################################################

class Ui(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):

        # call inherited classes __init__ method
        super(Ui, self).__init__(*args, **kwargs)    
        
        # Load the .ui file
        uic.loadUi('MESS.ui', self)                
        self.title = "Meteor Elemental Spectra Software"
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)

        # display the GUI 
        self.show()                                  

        ###########################################################################################
        ################################ /// GUI MODIFICATIONS /// ################################
        ###########################################################################################

        #################### SPECTRAL FILE IMAGE VIEW ####################

        self.scene = pg.GraphicsScene()
        # Create the widget
        self.spectral_imagewidget = pg.GraphicsView(self.scene)
        # Set Graphics View to ViewBox                                         
        self.spectral_imageframe = pg.ViewBox()
        # Set Graphics view as central widget                                               
        self.spectral_imagewidget.setCentralWidget(self.spectral_imageframe) 
        # Lock aspect ratio for the frame                   
        self.spectral_imageframe.setAspectLocked()        
        # Disable menu                      
        self.spectral_imageframe.setMenuEnabled(False)   
        # Invert the image along the y axis                       
        self.spectral_imageframe.invertY()
        # Add image item to ViewBox
        self.spectral_image = pg.ImageItem()           
        # add item to image frame (ViewBox)                                       
        self.spectral_imageframe.addItem(self.spectral_image)                                   
        # Location of widget in layout
        self.Spectral_layout.addWidget(self.spectral_imagewidget, 0, 0)

        #################### SPECTRAL ROI IMAGE VIEW ####################

        self.spectralROI_imagewidget = pg.GraphicsView()
        self.spectralROI_imageframe = pg.ViewBox()
        self.spectralROI_imagewidget.setCentralWidget(self.spectralROI_imageframe)
        self.spectralROI_imageframe.setMenuEnabled(False)
        self.spectralROI_image = pg.ImageItem()
        self.spectralROI_imageframe.addItem(self.spectralROI_image)
        self.SpectralROI_layout.addWidget(self.spectralROI_imagewidget,0,0)

        #################### DIRECT IMAGE HISTOGRAM ####################

        # Create histogram widget
        self.spectral_hist = pg.HistogramLUTWidget()   
        # Connect histogram to image in ViewBox                                       
        self.spectral_hist.setImageItem(self.spectral_image)                                    
        # location of widget in layout
        self.Spectral_layout.addWidget(self.spectral_hist, 0, 20)                               

        #################### INITIALIZE MOUSE ####################

        self.spectral_roi = None
        # Change to crosshair cursor in the direct image view widget
        self.spectral_imagewidget.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        # Enable mouse click event in the image view widget
        # self.spectral_image.mousePressEvent = self.getSpectralPosition
        self.spectralROI_image.mousePressEvent = self.getSpectralPosition


        #################### SPECTRAL ROI MARKERS ####################

        # Init affine marker on spectral image
        self.spectralROI_markers = pg.ScatterPlotItem()
        self.spectralROI_markers.setData(pxMode=False, symbol='+', size=10, pen='r', brush='r')
        self.spectralROI_imageframe.addItem(self.spectralROI_markers) 

        #################### DIRECT MARKERS ####################

        # Init affine marker on spectral image
        self.spectral_markers = pg.ScatterPlotItem()
        self.spectral_markers.setData(pxMode=False, symbol='+', size=15, pen='b', width=5, brush=None)
        self.spectral_imageframe.addItem(self.spectral_markers) 

        #################### INITIALIZE PLOT MOUSE #####################

        self.Plot.scene().sigMouseClicked.connect(self.mouse_clicked)
                        
        # ################# INITIALIZE SPECTRAL ROI/MOUSE #################

        # # Initalize ROI to None
        # self.spectral_roi = None
        # # Change to crosshair cursor in the spectral image widget
        # self.spectral_imagewidget.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))

        #################### SPECTRAL MARKERS ####################

        # Init affine marker on spectral image
        # self.affine_markers = pg.ScatterPlotItem()
        # self.affine_markers.setData(pxMode=False, symbol='o', size=10, pen='r', brush='r')
        # self.spectralROI_imageframe.addItem(self.affine_markers)

        # Init projected affine marker  on  spectral  image
        self.affine_markers = pg.ScatterPlotItem()
        self.affine_markers.setData(pxMode=False, symbol='+', size=10, pen='b', brush='b')
        self.spectralROI_imageframe.addItem(self.affine_markers)  
        
        ################# LOAD SPECTRAL FLAT ####################

        # Init file path and name
        # flat_path = self.FlatPath_linedit.text()
        # flat_name = self.FlatName_linedit.text()

        # self.FlatName_linedit.returnPressed.connect(lambda: self.updateFlatName())

        # Call the function to load the flat
        # if flat_name != '':
        #     self.flat_structure  = loadFlat(flat_path, flat_name)

        # initialize the GuralSpectral object
        self.spectral = spectral_library.GuralSpectral(10000, 4500, None, None, None, None)

        ###########################################################################################
        ################################# /// BUTTON TRIGGERS /// #################################
        ###########################################################################################

        ################# Spectral FILE CONTROL BUTTONS #################

        # Upload spectral file 
        self.UploadSpectral_button.clicked.connect(self.uploadSpectralVid)                      
        # Next spectral frame        
        self.NextSpectral_button.clicked.connect(self.nextSpectralFrame)  
        # Jump ahead 5 frames
        self.ForwardFiveSpectral_button.clicked.connect(self.forwardFiveSpectralFrames)                      
        # Last spectral frame        
        self.LastSpectral_button.clicked.connect(self.lastSpectralFrame)   
        # Jump back 5 frames
        self.BackFiveSpectral_button.clicked.connect(self.backFiveSpectralFrames)

        # Rotate spectral vid file
        # self.RotateVid_button.clicked.connect(self.rotateVid)
        # self.SpectralRotation_rollbox.valueChanged.connect(self.rotateVid)
        self.SpectralRotation_rollbox.setToolTip('Vid file frames will be rotated CW by this amount before display.')


        # Run affine transform        
        # self.Affine_button.clicked.connect(self.affineTransform)
        # Update affine transform
        # self.UpdateAffine_button.clicked.connect(self.updateTransform)
        # Pick a feature for the transform
        # self.PickFeature_button.clicked.connect(self.pickFeature)        


        # ################# SPECTRAL FILE CONTROL BUTTONS #################
    
        # Load files from event.txt file
        self.LoadEventFile_button.clicked.connect(self.loadEventFile)
        
        # Make ROI box appear        
        self.SelectSpectralRegion_button.clicked.connect(self.spectralROI)                        

        # Clear ROI box        
        self.ClearSpectralRegion_button.clicked.connect(self.clearSpectralROI)             
        # Clear the affine transform point
        # self.ClearAffine_button.clicked.connect(self.clearAffine)
        # Apply the image flat
        self.RemoveFlat_button.clicked.connect(self.removeSpectralFlat)

        # Auto pick ROI
        self.AutoPick_button.clicked.connect(self.autoPickROI)
        # Auto pick meteor
        self.AutoPickDirect_button.clicked.connect(self.autoPickDirect)

        self.AutoTrackDirect_button.clicked.connect(self.autoTrackDirect)

        ################# SpectralFlat
        self.UploadSpectralFlat_button.clicked.connect(self.uploadSpectralFlat)
        self.AutoSpectralFlat_button.clicked.connect(self.autoSpectralFlat)

        ################# PLOTTING BUTTONS #################

        # Plot the measured spectrum    
        self.MeasuredSpec_button.clicked.connect(self.plotMeasuredSpec)
        # Calibrate the spectrum
        self.CalibrateSpectrum_button.clicked.connect(lambda: self.calibrationClicked())  
        # Set reference spectrum
        self.SetReference_button.clicked.connect(self.setReference)            
        # Clear the plot    
        self.Clear_button.clicked.connect(self.clearSpec)
        # Save the plot
        self.SavePlot_button.clicked.connect(self.savePlot)
        #MJM

        #################### Elemental Abundance Buttons #################

        self.FitMeasuredSpectrum_button.clicked.connect(self.fitMeasuredSpectrum)
        self.ResetElements_button.clicked.connect(self.resetAllElementalAbundances)

        ######################## Commands Buttons #########################

        # Rollboxes
        self.Extinction_rollbox.valueChanged.connect(lambda: self.updateExtinctionValue())
        self.Roll_rollbox.valueChanged.connect(lambda: self.updateRollValue())
        self.Lmm_rollbox.valueChanged.connect(lambda: self.updateLmmValue())
        self.HighTemp_rollbox.valueChanged.connect(lambda: self.updateHighTempValue())
        self.LowTemp_rollbox.valueChanged.connect(lambda: self.updateLowTempValue())
        self.Sigma_rollbox.valueChanged.connect(lambda: self.updateSigmaValue())
        self.Hot2WarmRatio_rollbox.valueChanged.connect(lambda: self.updateHot2WarmRatio())

        self.MeteorSpeed_rollbox.valueChanged.connect(lambda: self.updateMeteorSpeed())
        self.MeteorHeight_rollbox.valueChanged.connect(lambda: self.updateMeteorHeight())
        self.ColumnDensity_rollbox.valueChanged.connect(lambda: self.updateColumnDensity())
        self.PlasmaRadius_rollbox.valueChanged.connect(lambda: self.updatePlasmaRadius())

        # Element buttons
        self.Na_button.clicked.connect(self.elementButtonClicked)
        self.Mg_button.clicked.connect(self.elementButtonClicked)
        self.Ca_button.clicked.connect(self.elementButtonClicked)
        self.Fe_button.clicked.connect(self.elementButtonClicked)
        self.K_button.clicked.connect(self.elementButtonClicked)
        self.O_button.clicked.connect(self.elementButtonClicked)
        self.N_button.clicked.connect(self.elementButtonClicked)
        self.Si_button.clicked.connect(self.elementButtonClicked)
        self.H_button.clicked.connect(self.elementButtonClicked)
        self.Cr_button.clicked.connect(self.elementButtonClicked)
        self.Mn_button.clicked.connect(self.elementButtonClicked)
        self.Al_button.clicked.connect(self.elementButtonClicked)
        self.Ti_button.clicked.connect(self.elementButtonClicked)
        self.Ni_button.clicked.connect(self.elementButtonClicked)
        self.C_button.clicked.connect(self.elementButtonClicked)
        self.Co_button.clicked.connect(self.elementButtonClicked)
        self.Cu_button.clicked.connect(self.elementButtonClicked)
        self.V_button.clicked.connect(self.elementButtonClicked)
        self.Sc_button.clicked.connect(self.elementButtonClicked)

        self.CaII_button.clicked.connect(self.elementButtonClicked)
        self.MgII_button.clicked.connect(self.elementButtonClicked)
        self.NaII_button.clicked.connect(self.elementButtonClicked)
        self.SiII_button.clicked.connect(self.elementButtonClicked)
        self.FeII_button.clicked.connect(self.elementButtonClicked)
        self.NII_button.clicked.connect(self.elementButtonClicked)
        self.CII_button.clicked.connect(self.elementButtonClicked)

        self.N2_button.clicked.connect(self.elementButtonClicked)
        self.O2_button.clicked.connect(self.elementButtonClicked)
        self.CH_button.clicked.connect(self.elementButtonClicked)
        self.CO2_button.clicked.connect(self.elementButtonClicked)

        self.CaO_button.clicked.connect(self.elementButtonClicked)
        self.FeO_button.clicked.connect(self.elementButtonClicked)
        self.MgO_button.clicked.connect(self.elementButtonClicked)
        self.AlO_button.clicked.connect(self.elementButtonClicked)

        # toggle-able command buttons
        self.HotTempOn_button.setCheckable(True)
        self.WarmTempOn_button.setCheckable(True)
        self.Ions_button.setCheckable(True)
        self.Neutral_button.setCheckable(True)

        self.HotTempOn_button.clicked.connect(lambda: self.hotTempToggle())
        self.HotTempOn_button.toggled.connect(self.refreshPlot)
        self.WarmTempOn_button.clicked.connect(lambda: self.warmTempToggle())
        self.WarmTempOn_button.toggled.connect(self.refreshPlot)
        self.Ions_button.clicked.connect(lambda: self.ionsToggle())
        self.Neutral_button.clicked.connect(lambda: self.neutralToggle())
        self.Responsivity_check.toggled.connect(lambda: self.responsivityToggle())
        self.Extinction_check.toggled.connect(lambda: self.extinctionToggle())

        self.RefreshPlot_button.clicked.connect(self.refreshPlot)
        self.CalculateFullSpectrum_button.clicked.connect(self.calculateFullSpectrum)

        noise_multiplier = 0.0 # 0 = no noise, 0.1 = defaults

        spectral_library.readSpectralConfig(self.spectral)
        spectral_library.allocMemory(self.spectral)

        # Assign camera numbers to the grating structure
        for i in range(spectral_library.MAXGRATINGS):
            self.spectral.spcalib.gratinfo[i].camnum = self.spectral.spconfig.camnum[i]

        camos_camera_index = 0;
        self.spectral.spcalib.gratinfo[camos_camera_index].grating_area_scale = math.cos(self.spectral.spconfig.grating_offnormal_deg * 3.141592654 / 180.0)
        self.spectral.spcalib.gratinfo[camos_camera_index].camnum = 101

        spectral_library.readSpectralCALFile(self.spectral)
        spectral_library.loadElementsData(self.spectral)

        # Update controls based on config
        # self.SpectralConfigVersion_label.setText(str(self.spectral.spconfig.version))
        # self.SpectralConfigOrder_label.setText(str(self.spectral.spconfig.order4spcalib))
        # self.SpectralConfigRowDelta_label.setText(str(self.spectral.spconfig.rowdelta))
        # self.SpectralConfigMinCalWavelengthNm_label.setText(str(self.spectral.spconfig.min_cal_wavelength_nm))
        # self.SpectralConfigMaxCalWavelengthNm_label.setText(str(self.spectral.spconfig.max_cal_wavelength_nm))
        # self.SpectralConfigDelWavelengthNm_label.setText(str(self.spectral.spconfig.del_wavelength_nm))
        # self.SpectralConfigMinFitWavelengthNm_label.setText(str(self.spectral.spconfig.min_fit_wavelength_nm))
        # self.SpectralConfigMaxFitWavelengthNm_label.setText(str(self.spectral.spconfig.max_fit_wavelength_nm))
        # self.SpectralConfigMinspacWavelengthNm_label.setText(str(self.spectral.spconfig.minspan_wavelength_nm))
        # self.SpectralConfigSmoothBandwidthNm_label.setText(str(self.spectral.spconfig.smooth_bandwidth_nm))
        # self.SpectralConfigFaintestStarMag_label.setText(str(self.spectral.spconfig.faintest_star_vmag))
        # self.SpectralConfigAirmassLimit_label.setText(str(self.spectral.spconfig.airmass_limit))
        # self.SpectralConfigFadingCoef_label.setText(str(self.spectral.spconfig.fading_coef))
        # self.SpectralConfigCoinTimeTolerance_label.setText(str(self.spectral.spconfig.coin_time_tolerance))
        # self.SpectralConfigMinLoExcTemp_label.setText(str(self.spectral.spconfig.min_lo_exc_temp))
        # self.SpectralConfigMaxLoExcTemp_label.setText(str(self.spectral.spconfig.max_lo_exc_temp))
        # self.SpectralConfigStepLoExcTemp_label.setText(str(self.spectral.spconfig.step_lo_exc_temp))
        # self.SpectralConfigNominalLoExcTemp_label.setText(str(self.spectral.spconfig.nominal_lo_exc_temp))
        # self.SpectralConfigGratingOffnormal_label.setText(str(self.spectral.spconfig.grating_offnormal_deg))
        # self.SpectralConfigDefaultRollAngle_label.setText(str(self.spectral.spconfig.default_roll_deg))
        # self.SpectralConfigDefaultPitchAngle_label.setText(str(self.spectral.spconfig.default_pitch_deg))
        # self.SpectralConfigDefaultYawDeg_label.setText(str(self.spectral.spconfig.default_yaw_deg))
        # self.SpectralConfigDefaultNe_label.setText(str(self.spectral.spconfig.default_ne))
        # self.SpectralConfigDefaultHot2warm_label.setText(str(self.spectral.spconfig.default_hot2warm))
        # self.SpectralConfigNcamsGrating_label.setText(str(self.spectral.spconfig.ncams_grating))
        # self.SpectralConfigNominalHiExcTemp_label.setText(str(self.spectral.spconfig.nominal_hi_exc_temp))
        # self.SpectralConfigNominalSigma0_label.setText(str(self.spectral.spconfig.nominal_sigma0))

        self.LowTemp_rollbox.setValue(int(self.spectral.spconfig.nominal_lo_exc_temp))
        self.HighTemp_rollbox.setValue(int(self.spectral.spconfig.nominal_hi_exc_temp))
        self.Sigma_rollbox.setValue(self.spectral.spconfig.nominal_sigma0)
        self.Hot2WarmRatio_rollbox.setValue(self.spectral.spconfig.default_hot2warm)

        # Define elements for fitting
        self.elementButtons = []
        self.elementDeets = [] # element name, fit state, atomic number, index
        self.kelem = {}
        self.fitState = {}
        self.elementButtons.append(self.Na_button)
        self.elementDeets.append(['Na',0,11,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 11)])
        self.elementButtons.append(self.Mg_button)
        self.elementDeets.append(['Mg',0,12,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 12)])
        self.elementButtons.append(self.Ca_button)
        self.elementDeets.append(['Ca',0,20,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 20)])
        self.elementButtons.append(self.Fe_button)
        self.elementDeets.append(['Fe',0,26,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 26)])
        self.elementButtons.append(self.K_button)
        self.elementDeets.append(['K',0,19,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 19)])
        self.elementButtons.append(self.O_button)
        self.elementDeets.append(['O',0,8,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 8)])
        self.elementButtons.append(self.N_button)
        self.elementDeets.append(['N',0,7,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 7)])
        self.elementButtons.append(self.N2_button)
        self.elementDeets.append(['N2',0,7,spectral_library.GuralSpectral.getElementIndex(self.spectral, 50, 7)])
        self.elementButtons.append(self.Si_button)
        self.elementDeets.append(['Si',0,14,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 14)])
        self.elementButtons.append(self.H_button)
        self.elementDeets.append(['H',0,1,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 1)])
        self.elementButtons.append(self.Mn_button)
        self.elementDeets.append(['Mn-I',0,25,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 25)])
        self.elementButtons.append(self.Cr_button)
        self.elementDeets.append(['Cr-I',0,24,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 24)])
        self.elementButtons.append(self.Al_button)
        self.elementDeets.append(['Al-I',0,13,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 13)])
        self.elementButtons.append(self.Ti_button)
        self.elementDeets.append(['Ti-I',0,22,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 22)])
        # self.elementButtons.append(self.H_button)
        # self.elementDeets.append(['H-I',0,1,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 1)])
        self.elementButtons.append(self.C_button)
        self.elementDeets.append(['C-I',0,6,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 6)])
        self.elementButtons.append(self.Sc_button)
        self.elementDeets.append(['Sc-I',0,21,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 21)])
        self.elementButtons.append(self.V_button)
        self.elementDeets.append(['V-I',0,23,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 23)])
        self.elementButtons.append(self.Co_button)
        self.elementDeets.append(['Co-I',0,27,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 27)])
        self.elementButtons.append(self.Ni_button)
        self.elementDeets.append(['Ni-I',0,28,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 28)])
        self.elementButtons.append(self.Cu_button)
        self.elementDeets.append(['Cu-I',0,29,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 29)])


        self.elementButtons.append(self.CaII_button)
        self.elementDeets.append(['CaII',0,20,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 20)])
        self.elementButtons.append(self.MgII_button)
        self.elementDeets.append(['MgII',0,12,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 12)])
        self.elementButtons.append(self.NaII_button)
        self.elementDeets.append(['NaII',0,11,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 11)])
        self.elementButtons.append(self.SiII_button)
        self.elementDeets.append(['Si-II',0,14,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 14)])

        # Extra elements
        self.elementDeets.append(['Al-II',0,13,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 13)])
        # self.elementDeets.append(['Al-I',0,13,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 13)])
        self.elementDeets.append(['Ar-II',0,18,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 18)])
        self.elementDeets.append(['Ar-I',0,18,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 18)])
        self.elementDeets.append(['Ba-II',0,56,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 56)])
        self.elementDeets.append(['Ba-I',0,56,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 56)])
        self.elementDeets.append(['Be-II',0,4,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 4)])
        self.elementDeets.append(['Be-I',0,4,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 4)])
        # self.elementDeets.append(['Ca-II',0,20,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 20)])
        self.elementDeets.append(['C-II',0,6,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 6)])
        self.elementDeets.append(['Cl-II',0,17,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 17)])
        self.elementDeets.append(['Cl-I',0,17,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 17)])
        self.elementDeets.append(['Co-II',0,27,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 27)])
        self.elementDeets.append(['Cr-II',0,24,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 24)])
        # self.elementDeets.append(['Cr-I',0,24,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 24)])
        self.elementDeets.append(['Cs-II',0,55,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 55)])
        self.elementDeets.append(['Cs-I',0,55,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 55)])
        self.elementDeets.append(['Cu-II',0,29,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 29)])
        self.elementDeets.append(['Cu-I',0,29,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 29)])
        self.elementDeets.append(['Fe-II',0,26,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 26)])
        self.elementDeets.append(['F-II',0,9,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 9)])
        self.elementDeets.append(['F-I',0,9,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 9)])
        self.elementDeets.append(['Ge-II',0,32,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 32)])
        self.elementDeets.append(['Ge-I',0,32,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 32)])
        self.elementDeets.append(['H-II',0,1,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 1)])
        # self.elementDeets.append(['H-I',0,1,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 1)])
        self.elementDeets.append(['K-II',0,19,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 19)])
        self.elementDeets.append(['Li-II',0,3,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 3)])
        self.elementDeets.append(['Li-I',0,3,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 3)])
        # self.elementDeets.append(['Mg-II',0,12,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 12)])
        self.elementDeets.append(['Mn-II',0,25,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 25)])
        # self.elementDeets.append(['Mn-I',0,25,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 25)])
        self.elementDeets.append(['Mn-II',0,12,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 12)])
        self.elementDeets.append(['N2-II',0,7,spectral_library.GuralSpectral.getElementIndex(self.spectral, 51, 7)])
        self.elementDeets.append(['N2-I',0,7,spectral_library.GuralSpectral.getElementIndex(self.spectral, 50, 7)])
        # self.elementDeets.append(['Na-II',0,11,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 11)])
        self.elementDeets.append(['Ne-II',0,10,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 10)])
        self.elementDeets.append(['Ne-I',0,10,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 10)])
        self.elementDeets.append(['Ni-II',0,28,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 28)])
        self.elementDeets.append(['Ni-I',0,28,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 28)])
        self.elementDeets.append(['O-II',0,8,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 8)])
        self.elementDeets.append(['Pb-II',0,82,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 82)])
        self.elementDeets.append(['Pb-I',0,82,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 82)])
        self.elementDeets.append(['P-II',0,15,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 15)])
        self.elementDeets.append(['P-I',0,15,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 15)])
        self.elementDeets.append(['Rb-II',0,37,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 37)])
        self.elementDeets.append(['Rb-I',0,37,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 37)])
        self.elementDeets.append(['Sc-II',0,21,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 21)])
        self.elementDeets.append(['Sc-I',0,21,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 21)])
        self.elementDeets.append(['Sc-II',0,21,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 21)])
        
        self.elementDeets.append(['S-II',0,16,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 16)])
        self.elementDeets.append(['S-I',0,16,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 16)])
        self.elementDeets.append(['Sr-II',0,38,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 38)])
        self.elementDeets.append(['Sr-I',0,38,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 38)])
        self.elementDeets.append(['Ti-II',0,22,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 22)])
        self.elementDeets.append(['Ti-I',0,22,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 22)])
        self.elementDeets.append(['V-II',0,23,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 23)])
        self.elementDeets.append(['V-I',0,23,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 23)])
        self.elementDeets.append(['Y-II',0,39,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 39)])
        self.elementDeets.append(['Y-I',0,39,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 39)])
        self.elementDeets.append(['Zn-II',0,30,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 30)])
        self.elementDeets.append(['Zn-I',0,30,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 30)])
        self.elementDeets.append(['Zr-II',0,40,spectral_library.GuralSpectral.getElementIndex(self.spectral, 1, 40)])
        self.elementDeets.append(['Zr-I',0,40,spectral_library.GuralSpectral.getElementIndex(self.spectral, 0, 40)])

        self.elemIndex = 0

        # ============================================================================
        #    At this stage, you should extract an integrated spectrum from the 
        #       imagery (for a frame or aggregated frames) that is to be fit and
        #       place the spectrum in the vector "spectrum.integ_spec" at the 
        #       corresponding vector of wavelengths of spcal.wavelength_nm.
        #       NOTE: For the purposes of this example we will infill integ_spec
        #             later with a model spectrum plus noise to be fit.
        # 
        #    You will also need the corresponding metadata of the event such as
        #       heights above earth's surface, range to the meteor, the approach 
        #       angle (look direction angle off the radiant), entry velocity, 
        #       altitude angle (elevation above the horizon).   
        # ============================================================================

        # self.vinfinity_kmsec = 40.0     # km/sec
        # self.approach_angle_radians = 55.0 * 3.141592654 / 180.0
        # self.earth_radius_km = 6378.16  # WGS84
        # self.site_height_km = 0.2       # Above WGS84 Earth surface
        # self.meteor_height_km = 85.0    # Above WGS84 Earth surface
        # self.altitude_deg = 45.0        # elevation angle above horizon
        # Rsin = self.earth_radius_km * math.sin(self.altitude_deg * 3.141592654 / 180.0)
        # self.meteor_range_km = math.sqrt( Rsin * Rsin + 2.0 * self.earth_radius_km * self.meteor_height_km + self.meteor_height_km * self.meteor_height_km) - Rsin

        #========== Set the grating camera with its scaling factor based off the config file input. 
        #           CAMO-S has one grating camera so only index [0] is set. 
        #           On other systems this scaling factor varies with each event due the different angle
        #           of incidence of the light path onto the grating.
        camos_camera_index = 0
        self.spectral.spcalib.gratinfo[camos_camera_index].grating_area_scale = math.cos(self.spectral.spconfig.grating_offnormal_deg * 3.141592654 / 180.0)
        print('Grating area scale: %s' % self.spectral.spcalib.gratinfo[camos_camera_index].grating_area_scale)

        #========== Set user adjustable values in the elemdata structure to their starting defaults
        #              such as sigma, temperatures, electron density, hot-to-warm, airmass factor
        spectral_library.adjustableParametersDefaults(self.spectral)
        self.spectral.vinfinity_kmsec = 35.0
        spectral_library.GuralSpectral.resetAllElementalAbundances(self.spectral)

        self.PlottedSpectrumNumber = 0

        # Set rollboxes to default values
        self.MeteorHeight_rollbox.setValue(self.spectral.meteor_height_km)
        self.MeteorSpeed_rollbox.setValue(self.spectral.vinfinity_kmsec)
        self.ColumnDensity_rollbox.setValue(self.spectral.elemdata.ne_jones/1E13)
        self.ColumnDensityIter_rollbox.setValue(self.spectral.elemdata.ne_iter/1E13)
        self.PlasmaRadius_rollbox.setValue(self.spectral.elemdata.plasma_radius_meters)


        self.Scale_label.setText(str(self.SpectralScale_rollbox.value()))
        self.nm0_label.setText(self.nm0_edit.text())

        self.extra_elements = ['Al-II', 'Ar-II', 'Ar-I', 'Ba-II', 'Ba-I', 'Be-II', 'Be-I', \
        'Cl-II', 'Cl-I', 'Co-II', 'Cr-II', 'Cr-I', 'Cs-II', 'Cs-I', 'Cu-II', \
        'F-II', 'F-I', 'Ge-II', 'Ge-I', 'H-II', 'K-II', 'Li-II', 'Li-I', 'Mg-II', 'Mn-II',\
        'N2-II', 'N2-I', 'Na-II', 'Ne-II', 'Ne-I', 'Ni-II', 'O-II', 'Pb-II', 'Pb-I', 'P-II', \
        'P-I', 'Rb-II', 'Rb-I', 'Sc-II', 'S-I', 'Sr-II', 'Sr-I', 'Ti-II', 'V-II', \
        'Y-II', 'Y-I', 'Zn-II', 'Zn-I', 'Zr-II', 'Zr-I']
        self.ExtraElements_combo.addItems(self.extra_elements)
        self.ExtraElements_combo.currentIndexChanged.connect(self.plotExtraElement)
        self.extraElementsCount = 0
        self.DisplayedExtraElements_combo.currentIndexChanged.connect(self.removeExtraElement)
        # self.PlotExtraElements_button.clicked.connect(self.plotExtraElement)

        ########## Setup Stellar Calibration Tab ############
        spectral_library.readStarSpectra(self.spectral)
        starlist = ['Alderamin', 'Caph', 'Dubhe', 'Kappa Dra', 'Kochab', 'Megrez', 'Merak', 'Polaris', 'Shedir', 'Thuban', 'Yildun']
        self.ChooseStar_combo.addItems(starlist)

        self.UploadStar_button.clicked.connect(self.uploadStarVid)
        self.UploadStarFlat_button.clicked.connect(self.uploadStarFlat)
        self.PlotStar_button.clicked.connect(self.plotStar)
        self.ShowStellarSpectrum_button.clicked.connect(self.plotStellarSpectrum)
        self.PickStarReference_button.clicked.connect(lambda: self.stellarCalibrationClicked())
        self.ComputeCalibration_button.clicked.connect(self.computeCalibration)

        self.elementsState = np.zeros(3) # unlocked, fitting, locked

        try:
            self.responsivityDefault = np.genfromtxt('DefaultResponsivity-EMCCD.csv', delimiter=',', skip_header=2)
            # self.responsivityDefault = np.genfromtxt('DefaultResponsivity-EMCCD2.csv', delimiter=',', skip_header=2)
            # self.responsivityDefault = np.genfromtxt('FS9910.csv', delimiter=',', skip_header=2)
            print('Loaded default responsivity curve.')
        except:
            print('Unable to load default responsivity curve!')


        self.dir_x = 0
        self.dir_y = 0

    ###############################################################################################
    ###################################### /// FUNCTIONS /// ######################################
    ###############################################################################################

    def computeCalibration(self):

        self.calibration = np.divide(self.stellarIntensity, self.spectral_profile)
        self.calibrationX = self.spectral_x

        # Set axis titles 
        self.Plot.setLabel('left', 'Intensity')
        self.Plot.setLabel('bottom', 'Wavelength (nm)')

        # Create the plot
        pen = pg.mkPen(width = 2)
        self.Plot.plot(self.calibrationX, self.calibration, pen = pen)

    def plotStellarSpectrum(self):
        # nwave, nstars, wave, star
        # hip, spectype, vmag, specval,

        compStars = {
            'Alderamin': 105199,
            'Caph': 746,
            'Dubhe': 54061,
            'Kappa Dra': 61281,
            'Kochab': 72607,
            'Megrez': 59774,
            'Merak': 53910,
            'Polaris': 11767,
            'Shedir': 3179,
            'Thuban': 68756,
            'Yildun': 85822}

        starHipNum = int(compStars[self.ChooseStar_combo.currentText()])

        # for i in range(self.spectral.starspectra.nstars):
        #     print(self.spectral.starspectra.wave[0])
        #     if self.spectral.starspectra.star[i].hip == starHipNum:
        #         print('Match: %s' %  self.spectral.starspectra.star[i].hip)
        #         self.spectral.star_index = i

        spectral_library.GuralSpectral.interpolateCatalogedSpectrum(self.spectral, starHipNum)

        star_x = []
        star_y = []
        for i in range(self.spectral.starspectra.nwave):
            star_x = np.append(star_x, self.spectral.starspectra.wave[i])
            star_y = np.append(star_y, self.spectral.starspectra.star[self.spectral.star_index].specval[i])

        self.calibMax = np.max(star_y)
        self.calScale = self.calibMax / self.starMax
        self.stellarIntensity = star_y

        print('Scaler: %s' % self.calScale)

        self.clearSpec()
        self.plotStar()

        print('x min: %s   x max: %s' % (np.min(star_x), np.max(star_x)))
        self.star_x_min = np.min(star_x)
        self.star_x_max = np.max(star_x)

        # index_min = np.where(star_x == np.round(self.scaled_spectral_x_min))
        # index_max = np.where(star_x == np.round(self.scaled_spectral_x_max))
        index_min = (np.abs(star_x - np.round(self.scaled_spectral_x_min))).argmin()
        index_max = (np.abs(star_x - np.round(self.scaled_spectral_x_max))).argmin()

        # print(star_x)

        print('Index of min observed: %s' % np.where(star_x == np.round(self.scaled_spectral_x_min)))
        print('Index of max observed: %s' % np.where(star_x == np.round(self.scaled_spectral_x_max)))

        self.star_x_clipped = star_x[index_min:index_max]
        # print(len(self.star_x_clipped))
        # print(len(self.spectral_x))

        # Set axis titles 
        self.Plot.setLabel('left', 'Intensity')
        self.Plot.setLabel('bottom', 'Wavelength (nm)')
        print('Calib Max: %s' % np.max(star_y))
        # Create the plot
        pen = pg.mkPen(width = 2)
        self.Plot.plot(star_x, star_y, pen = pen)
        # self.Plot.setXRange(np.min(scaled_spectral_profile),np.max(scaled_spectral_profile))
        # self.CalibrateSpectrum_button.setEnabled(True)
        # self.SetReference_button.setEnabled(True)

    def stellarCalibrationClicked(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_StellarCalibrationDialog()
        self.ui.setupUi(self.window)
        self.ui.nm0StarSet_button.clicked.connect(self.nm0StarSet)
        self.window.show()

    def columnDensityClicked(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_ColumnDensityDialog()
        self.ui.setupUi(self.window)
        self.ui.columnDensityAccept_button.clicked.connect(self.columnDensitySet)
        self.window.show()

    def columnDensitySet(self):
        self.spectral.elemdata.els[self.elemIndex].N_warm = float(self.ui.columnDensity_edit.text()) * (10 ** self.ui.exponent_rollbox.value())
        print(float(self.ui.columnDensity_edit.text()), self.ui.exponent_rollbox.value())
        print('N warm: %s' % self.spectral.elemdata.els[self.elemIndex].N_warm)

    def nm0StarSet(self):
        self.nm0Star = float(self.ui.nm0Star_edit.text())
        self.x0Star = float(self.ui.x0Star_label.text())
        self.clearSpec()
        self.plotStar()
        # print(self.x0Star)

    def calibrationClicked(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_CalibrationDialog()
        self.ui.setupUi(self.window)
        self.ui.CalculateScale_button.clicked.connect(self.calculateScale)
        self.ui.UpdateScale_button.clicked.connect(self.updateScale)
        self.ui.nm0Set_button.clicked.connect(self.nm0Set)
        self.window.show()

    def nm0Set(self):
        # nm0 =int(self.ui.nm0_edit.text())
        self.nm0_label.setText(self.ui.nm0_edit.text())
        self.nm0_edit.setText(self.ui.nm0_edit.text())
        self.clearSpec()
        self.plotMeasuredSpec()

    def uploadStarVid(self, byteswap=True, dark=None):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec():
            spectral_file_name = dlg.selectedFiles()

            if spectral_file_name[0].endswith('.vid'):
                
                star_path = os.path.split(spectral_file_name[0])[0]
                star_name = os.path.split(spectral_file_name[0])[1]

                self.spectral_vid = readVid(star_path, star_name)
                self.spectral_currentframe = int(len(self.spectral_vid.frames)/2)
                self.spectral_vidlength = len(self.spectral_vid.frames)

                self.updateStarFrames()
                self.autoPickROI()
                # self.uploadStarFlat()

    def uploadStarFlat(self, byteswap=True, dark=None):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec():
            flat_file_name = dlg.selectedFiles()
            
            if flat_file_name[0].endswith('.png'):
                # print(flat_file_name[0])
               
                # Load the flat image
                flat_img = loadImage(flat_file_name[0], -1)

                # # Change the file type if given
                # if dtype is not None:
                #     flat_img = flat_img.astype(dtype)
                #  # If the flat isn't a 8 bit integer, convert it to uint16
                # elif flat_img.dtype != np.uint8:
                #     flat_img = flat_img.astype(np.uint16)

                if flat_img.dtype != np.uint8:
                    flat_img = flat_img.astype(np.uint16)

                # if byteswap:
                #     flat_img = flat_img.byteswap()
                #     print('swapped bytes')

                flat_img = flat_img.byteswap()

                # Init a new Flat structure
                self.flat_structure = FlatStruct(flat_img, dark=dark)

                self.removeStarFlat()
            else:
                pass

    def removeStarFlat(self):
        """
        Allows user to remove flat from the current frame only. 
        Flat can be restored by moving to next/last frame, and then 
        returning to desired frame. 

        Background and spectrum will still be calculted with flat applied,
        the flat can only be removed for user convenience.
        """

        # Set frame 
        # self.spectral_frame = self.spectral_vid.frames[self.spectral_currentframe]
        # self.spectral_frame_img = self.spectral_frame.img_data

        self.spectral_frame_img = applyFlat(self.spectral_frame_img, self.flat_structure)

        # # Display time
        # self.st = unixTime2Date(self.spectral_frame.ts, self.spectral_frame.tu, dt_obj=False)
        # self.st = str(self.st)
        # self.SpectralTime_label.setText(self.st)
        # self.update()

        # # Display frame number
        # self.SpectralFrame_label.setNum(self.spectral_currentframe)
        # self.update()

        # Set image levels
        minv = np.percentile(self.spectral_frame_img, 0.1)
        maxv = np.percentile(self.spectral_frame_img, 99.95)
        gamma = 1

        # Create an image with properly adjust levels
        spectral_frame_img = adjustLevels(self.spectral_frame_img, minv, gamma, maxv, scaleto8bits=True)

        # Set spectral image
        self.spectral_image.setImage(spectral_frame_img.T)

    def updateStarFrames(self):

        self.avg_frame_img = np.zeros((self.spectral_vid.frames[0].img_data.shape[0],self.spectral_vid.frames[0].img_data.shape[1]))
        for i in range(len(self.spectral_vid.frames)):
            self.avg_frame_img += self.spectral_vid.frames[i].img_data

        self.avg_frame_img = self.avg_frame_img / len(self.spectral_vid.frames)
        self.avg_frame_img = self.avg_frame_img.astype(int)

        self.spectral_frame_img = self.avg_frame_img

        # Display time
        # self.st = unixTime2Date(self.spectral_frame.ts, self.spectral_frame.tu, dt_obj=False)

        # self.SpectralTime_label.setText(' at ' + str(self.st[3]) + ':' + str(self.st[4]) + \
        #  ':' + str(self.st[5]) + '.' + str(self.st[6]) + 'UT on ' + str(self.st[0]) + '/' + \
        #  str(self.st[1]) + '/' + str(self.st[2]))
        # self.update()

        # # Display frame number
        # # self.SpectralFrame_label.setNum(self.spectral_currentframe)
        # self.SpectralFrame_label.setText('Frame # ' + str(self.spectral_currentframe))
        # self.update()

        # Call the function to apply the flat
        # if self.flat_structure is not None: 
        #     self.spectral_frame_img = applyFlat(self.spectral_frame_img, self.flat_structure)

        # Set spectral image
        self.spectral_image.setImage(self.spectral_frame_img.T)

    def plotStar(self):

        self.spectral_array = self.spectral_roi.getArrayRegion(self.spectral_frame_img.T, self.spectral_image)


        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        # plt.subplots(figsize=(676*px, 100*px))
        # plt.imshow(spectral_array_unzoomed.T)
        # plt.show()

        # Set spectral profile
        self.spectral_profile = np.sum(self.spectral_array, axis=1)
        self.spectral_profile = self.spectral_profile[5:-5]

        # Init array for the scaled profile
        global scaled_spectral_profile
        scaled_spectral_profile = np.zeros(len(self.spectral_profile))

        # Scaling parameters
        s = 1/2.15 # px/nm
        # s = self.SpectralScale_rollbox.value() * self.spectral_vid.frames[0].img_data.shape[1]/676.0 # px/nm
        # nm0 = float(self.nm0_edit.text()) # nm 
        print('s = %s' % s)

        try:
            if self.ui.x0Star_label.text() == 'not set':
                pass
        except:
            self.x0Star = 0
            self.nm0Star = 0
            self.calScale = 1

        self.dnm = (self.nm0Star-self.x0Star)


        # Calculate wavelength values as they correspond to each pixel
        for i in range(len(scaled_spectral_profile)):
            # nmt = (((i - self.x) / s) + nm0)
            nmt = i/s + 316.0 + self.dnm
            scaled_spectral_profile = np.append(scaled_spectral_profile, nmt)

        # Scale observed spectrum
        self.starMax = np.max(self.spectral_profile)
        self.spectral_profile = self.spectral_profile * self.calScale
        self.spectral_x = scaled_spectral_profile

        # Take the part of the array with the desired values
        length = len(scaled_spectral_profile)
        middle_index = length//2
        
        # Reset array
        scaled_spectral_profile = scaled_spectral_profile[middle_index:]

        print('x min: %s  x max: %s' % (np.min(scaled_spectral_profile), np.max(scaled_spectral_profile)))
        self.scaled_spectral_x_min = np.min(scaled_spectral_profile)
        self.scaled_spectral_x_max = np.max(scaled_spectral_profile)

        # Set axis titles 
        self.Plot.setLabel('left', 'Intensity')
        self.Plot.setLabel('bottom', 'Wavelength (nm)')

        # Create the plot
        pen = pg.mkPen('r', width=2)
        self.Plot.plot(scaled_spectral_profile, self.spectral_profile, pen = pen)
        self.Plot.setXRange(np.min(scaled_spectral_profile),np.max(scaled_spectral_profile))
        # self.CalibrateSpectrum_button.setEnabled(True)
        # self.SetReference_button.setEnabled(True)

    def calculateScale(self):

        old_scale = self.ui.NewScale_rollbox.value()
        w1 = int(self.ui.Wave1_edit.text())
        w2 = int(self.ui.Wave2_edit.text())
        x1 = float(self.ui.CalibX1_label.text())
        x2 = float(self.ui.CalibX2_label.text())
        new_scale = old_scale / np.abs((w2-w1)/(x2-x1))
        self.ui.NewScale_rollbox.setValue(new_scale)
        nm0 = float(self.ui.nm0_edit.text())
        self.ui.UpdateScale_button.setEnabled(True)

    def updateMeteorHeight(self):
        print('Updating meteor height...')

    def updateMeteorSpeed(self):
        print('Updating meteor speed...')

    def updateColumnDensity(self):
        print('Updating column density...')

    def updatePlasmaRadius(self):
        print('Updating plasma radius...')

    def updateScale(self):
        self.SpectralScale_rollbox.setValue(self.ui.NewScale_rollbox.value())
        self.nm0_label.setText(self.ui.nm0_edit.text())
        self.nm0_edit.setText(self.ui.nm0_edit.text())
        self.Scale_label.setText(str(self.ui.NewScale_rollbox.value()))
        self.plotMeasuredSpec()

    def pickFeature(self):
        print('Pick feature set.')
        self.spectralROI_image.mousePressEvent = self.getSpectralPosition

    def mouse_clicked(self, evt):

        vb = self.Plot.plotItem.vb
        scene_coords = evt.scenePos()
        if self.Plot.sceneBoundingRect().contains(scene_coords):
            mouse_point = vb.mapSceneToView(scene_coords)
            print(f'clicked plot X: {mouse_point.x()}, Y: {mouse_point.y()}, event: {evt}')
            self.statusBar.showMessage(f'clicked plot X: {mouse_point.x()}, Y: {mouse_point.y()}')

            try:
                if self.ui.CalibX1_label.text() == 'not set':
                    global num_clicks
                    num_clicks = 0
                if num_clicks == 0:
                    self.ui.CalibX1_label.setText(str(mouse_point.x()))
                    self.statusBar.showMessage('Setting X1 to %f' % mouse_point.x())
                    num_clicks = 1
                else:
                    self.ui.CalibX2_label.setText(str(mouse_point.x()))
                    self.statusBar.showMessage('Setting X2 to %f' % mouse_point.x())
                    num_clicks = 0
            except:
                pass

            try:
                self.ui.x0Star_label.setText(str(mouse_point.x()))
                self.statusBar.showMessage('Setting x0 to %f' % mouse_point.x())
            except:
                pass

    def plotExtraElement(self):
        if self.ExtraElements_combo.currentIndex() == 0: return
        print('Plotting extra element')

        for i in range(len(self.elementDeets)):
            if self.elementDeets[i][0] == self.ExtraElements_combo.currentText():
                self.elemName = self.elementDeets[i][0]
                self.elemNumber = self.elementDeets[i][2]
                self.elemIndex = self.elementDeets[i][3]
                self.calculateElementSpectrum()

        # scaled_extra_element_array = np.zeros(len(self.element_array))

        # There are 64 extra elements. Let's make the colors...
        try:
            pen_color = (self.elemIndex,64)
        except:
            pen_color = 'w'

        plotName = str(self.elemName)

        try: globals()[plotName]
        except KeyError:  globals()[plotName] = None

        if globals()[plotName] is None:
            self.Plot.addLegend()
            self.Plot.enableAutoRange(axis='y', enable=True)
            self.DisplayedExtraElements_combo.addItem(plotName)
            globals()[plotName] = self.Plot.plot(self.element_array[:,0], self.element_array[:,2], pen=pg.mkPen(pen_color, width=2, style=QtCore.Qt.DotLine), name=self.elemName)
        else:
            # globals()[plotName].setData(self.element_array[:,0], self.element_array[:,2])
            self.Plot.removeItem(globals()[self.elemName])
            globals()[self.elemName] = None

        self.ExtraElements_combo.setCurrentIndex(0)
        self.extraElementsCount += 1

    def removeExtraElement(self):
        # if self.DisplayedExtraElements_combo.currentIndex == 0:
        if self.DisplayedExtraElements_combo.count() == self.extraElementsCount:
            self.Plot.removeItem(globals()[self.DisplayedExtraElements_combo.currentText()])
            self.DisplayedExtraElements_combo.removeItem(self.DisplayedExtraElements_combo.currentIndex())
            self.extraElementsCount -= 1
        else:
            return

    def plotElement(self, event):

        # scaled_element_array = np.zeros(len(self.element_array))

        # # Scaling parameters
        # s = 2.85 # px/nm
        # nm0 = 410 # nm

        if self.elemName == 'Na':
            pen_color = (255,255,0)
        elif self.elemName == 'K':
            pen_color = (238,130,238)
        elif self.elemName == 'Mg':
            pen_color = (0,255,0)
        elif self.elemName == 'O':
            pen_color = (255,69,0)
        elif self.elemName == 'N2':
            pen_color = (1,18)
        elif self.elemName == 'Si':
            pen_color = (2,18)
        elif self.elemName == 'Fe':
            pen_color = (30,144,255)
        elif self.elemName == 'Ca':
            pen_color = (3,18)
        elif self.elemName == 'N':
            pen_color = (139,0,0)
        elif self.elemName == 'Al':
            pen_color = (4,18)
        elif self.elemName == 'Cr':
            pen_color = (5,18)
        elif self.elemName == 'Mn':
            pen_color = (6,18)
        elif self.elemName == 'Ti':
            pen_color = (7,18)
        elif self.elemName == 'Ni':
            pen_color = (8,18)
        elif self.elemName == 'C':
            pen_color = (9,18)
        elif self.elemName == 'Co':
            pen_color = (10,18)
        elif self.elemName == 'Cu':
            pen_color = (11,18)
        elif self.elemName == 'V':
            pen_color = (12,18)
        elif self.elemName == 'Sc':
            pen_color = (13,18)
        elif self.elemName == 'H':
            pen_color = (255,0,0)

        elif self.elemName == 'CaII':
            pen_color = (160,32,240)
        elif self.elemName == 'FeII':
            pen_color = (14,18)
        elif self.elemName == 'MgII':
            pen_color = (15,18)
        elif self.elemName == 'SiII':
            pen_color = (210,180,140)
        elif self.elemName == 'NII':
            pen_color = (16,18)
        elif self.elemName == 'CII':
            pen_color = (18,18)

        # elif self.elemName == 'O2':
        #     pen_color = (17,23)
        # elif self.elemName == 'CH':
        #     pen_color = (17,23)
        # elif self.elemName == 'CO2':
        #     pen_color = (17,23)

        # elif self.elemName == 'CaO':
        #     pen_color = (18,23)
        # elif self.elemName == 'FeO':
        #     pen_color = (19,23)
        # elif self.elemName == 'MgO':
        #     pen_color = (20,23)
        # elif self.elemName == 'AlO':
        #     pen_color = (21,23)
        else:
            pen_color = 'k' 

        plotName = str(self.elemName)

        try:
            globals()[plotName]
        except KeyError: 
            globals()[plotName] = None

        if globals()[plotName] is None:
            self.Plot.addLegend()
            globals()[plotName] = self.Plot.plot(self.element_array[:,0], self.element_array[:,3], pen=pg.mkPen(pen_color, width=4), name=self.elemName)
        else:
            # self.Plot.enableAutoRange(axis='y', enable=True)
            # self.Plot.setAutoVisible(y=True)
            globals()[plotName].setData(self.element_array[:,0], self.element_array[:,3])  

    def refreshPlot(self):
        self.spectral.elemdata.hot2warm = self.Hot2WarmRatio_rollbox.value()
        self.spectral.elemdata.sigma0 = self.Sigma_rollbox.value()
        self.spectral.elemdata.Tlo = self.LowTemp_rollbox.value()
        self.spectral.elemdata.Thi = self.HighTemp_rollbox.value()
        self.spectral.elemdata.ne_jones = self.ColumnDensity_rollbox.value()
        self.spectral.elemdata.Xfactor = self.Extinction_rollbox.value()
        self.spectral.elemdata.plasma_radius_meters = self.PlasmaRadius_rollbox.value()
        self.spectral.meteor_height_km = self.MeteorHeight_rollbox.value()
        self.spectral.vinfinity_kmsec = self.MeteorSpeed_rollbox.value()


        spectral_library.GuralSpectral.changeHot2WarmRatio(self.spectral, self.spectral.elemdata.hot2warm)
        spectral_library.GuralSpectral.plasmaVolumes(self.spectral)
        spectral_library.GuralSpectral.computeWarmPlasmaSpectrum(self.spectral)
        spectral_library.GuralSpectral.computeHotPlasmaSpectrum(self.spectral)

        if self.Extinction_check.isChecked():
            spectral_library.GuralSpectral.extinctionModel(self.spectral)


        if self.HotTempOn_button.isChecked() and (self.WarmTempOn_button.isChecked() == False):
            # print('Hot is on. Warm is not.')
            
            self.element_array = np.zeros((self.spectral.spcalib.nwavelengths,4))
            for i in range(self.spectral.spcalib.nwavelengths):
                self.element_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
                self.element_array[i][1] = 0
                self.element_array[i][2] = self.spectral.elemdata.els[self.elemIndex].spechi[i]
                self.element_array[i][3] = self.element_array[i][2]

        elif (self.HotTempOn_button.isChecked() == False) and self.WarmTempOn_button.isChecked():
            # print('Warm is on. Hot is not.')
            # spectral_library.GuralSpectral.computeWarmPlasmaSpectrum(self.spectral)
            self.element_array = np.zeros((self.spectral.spcalib.nwavelengths,4))
            for i in range(self.spectral.spcalib.nwavelengths):
                self.element_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
                self.element_array[i][1] = self.spectral.elemdata.els[self.elemIndex].speclo[i]
                self.element_array[i][2] = 0
                self.element_array[i][3] = self.element_array[i][1]

        else:
            # print('Warm and Hot are both on.')
            # spectral_library.GuralSpectral.computeWarmPlasmaSpectrum(self.spectral)
            # spectral_library.GuralSpectral.computeHotPlasmaSpectrum(self.spectral)
            self.element_array = np.zeros((self.spectral.spcalib.nwavelengths,4))
            for i in range(self.spectral.spcalib.nwavelengths):
                self.element_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
                self.element_array[i][1] = self.spectral.elemdata.els[self.elemIndex].speclo[i]
                self.element_array[i][2] = self.spectral.elemdata.els[self.elemIndex].spechi[i]
                self.element_array[i][3] = self.element_array[i][1] + self.element_array[i][2]

        self.element_array[:,1] = self.element_array[:,1] * 10**self.Scale_rollbox.value()
        self.element_array[:,2] = self.element_array[:,2] * 10**self.Scale_rollbox.value()
        self.element_array[:,3] = self.element_array[:,3] * 10**self.Scale_rollbox.value()
        print('Max values: %s %s %s' % (np.max(self.element_array[:,1]), np.max(self.element_array[:,2]), np.max(self.element_array[:,3])))

        self.plotElement(self)

    def fitMeasuredSpectrum(self):
        """ Call fitMeasSpec from the Gural spectral library which does the following...
            1. Sets Jones' electron density by calling JonesElectronDensity from the C library
            2. Uses Jones' n_e as an initial guess for an iterative method (IterativeElectronDensity function)
            3. Spectral coefficients are then fit by calling FitSpectralCoefficients
            4. Column densities are then calculated with ColumnDensities_NumberAtoms
            5. The model spectrum is then calculated with SpectrumGivenAllCoefs
        """
        self.kelem_ref = spectral_library.GuralSpectral.setReferenceElem(self.spectral)

        fittingElems = []

        for i in range(self.spectral.elemdata.nelements):
            if (self.elementDeets[i][1] == 1):
                fittingElems.append(self.elementDeets[i][3])

        if len(fittingElems) == 2:
            spectral_library.GuralSpectral.fitMeasSpec(self.spectral, fittingElems[0], fittingElems[1])
        else:
            spectral_library.GuralSpectral.fitMeasSpec(self.spectral, fittingElems[0], fittingElems[1], fittingElems[2])

    def calculateFullSpectrum(self):
        self.kelem_ref = spectral_library.GuralSpectral.setReferenceElem(self.spectral)
        print('Reference Element: %s' % self.kelem_ref)


        # print(self.spectral.elemdata.kelem_ref)
        # print(self.spectral.elemdata.els[0].N_warm)
        fittingElems = []

        for i in range(self.spectral.elemdata.nelements):
            # print(self.elementDeets[i][0],self.elementDeets[i][3], self.elementDeets[i][1])
            if (self.elementDeets[i][1] == 1):
                fittingElems.append(self.elementDeets[i][3])

        print('%s elements being fit...' % len(fittingElems))
        print(fittingElems)

        if len(fittingElems) == 2:
            spectral_library.GuralSpectral.computeFullSpec(self.spectral, fittingElems[0], fittingElems[1])
        elif len(fittingElems) == 3:
            spectral_library.GuralSpectral.computeFullSpec(self.spectral, fittingElems[0], fittingElems[1], fittingElems[2])

        self.fullspec_array = np.zeros((self.spectral.spcalib.nwavelengths,2))

        for i in range(self.spectral.spcalib.nwavelengths):
            self.fullspec_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
            self.fullspec_array[i][1] = self.spectral.spectra.fit_spectrum[i]

        try:
            self.fullSpecScaler = self.plotMax / np.max(self.fullspec_array[:,1])
        except:
            self.fullSpecScaler = 1.0
            print('Could not scale spectrum.')

        self.plotFullSpectrum()

        spectral_library.GuralSpectral.writeFullSpectrum2(self.spectral, './FullSpectrumTest.txt')

        # print(self.spectral.elemdata.els[17].N_warm)
        print("self.spectral.ne %s" % self.spectral.ne)

        # for i in range(self.spectral.elemdata.nelements):
        #     print(self.spectral.elemdata.els[i].element_string, self.elementDeets[i][3], self.elementDeets[i][1])

    def plotFullSpectrum(self):

        plotName = 'Full Spectrum'

        try:
            globals()[plotName]
        except KeyError:
            globals()[plotName] = None

        if globals()[plotName] is None:
            self.Plot.addLegend()
            globals()[plotName] = self.Plot.plot(self.fullspec_array[:,0], self.fullspec_array[:,1]*self.fullSpecScaler, pen=pg.mkPen('k', width=4), name='Full')
        else:
            globals()[plotName].setData(self.fullspec_array[:,0], self.fullspec_array[:,1])
        
    def calculateElementSpectrum(self):
        self.spectral.spconfig.default_hot2warm = self.Hot2WarmRatio_rollbox.value()
        self.spectral.elemdata.hot2warm = self.Hot2WarmRatio_rollbox.value()
        self.spectral.elemdata.sigma0 = self.Sigma_rollbox.value()
        self.spectral.elemdata.Tlo = self.LowTemp_rollbox.value()
        self.spectral.elemdata.Thi = self.HighTemp_rollbox.value()
        self.spectral.elemdata.ne_jones = self.ColumnDensity_rollbox.value()
        self.spectral.elemdata.Xfactor = self.Extinction_rollbox.value()
        self.spectral.elemdata.plasma_radius_meters = self.PlasmaRadius_rollbox.value()
        self.spectral.meteor_height_km = self.MeteorHeight_rollbox.value()
        self.spectral.vinfinity_kmsec = self.MeteorSpeed_rollbox.value()

        #========== If any of the input arguments in the call below to PlasmaVolumes changes
        #              at a later time, you must call PlasmaVolumes again to infill the  
        #              elemdata structure values for height, range, and volumes. 
        #              For example, a different frame (height and range) is to be 
        #              processed, or the user adjusts the hot-to-warm ratio for 
        #              fitting purposes.
        spectral_library.GuralSpectral.plasmaVolumes(self.spectral)

        #========== Compute the model for extinction and the airmass given the event's metadata.
        #           You may want to update for each altitude change or keep it fixed for all frames.
        if self.Extinction_check.isChecked():
            spectral_library.GuralSpectral.extinctionModel(self.spectral)

        #========== Zero all the elemental abundances and #atoms
        #           Set all element fitting flags to not selected for fitting = FITLESS 
        #           Compute Jones 1997 fraction of ionized atoms Beta = n+ / ( n+ + no ) = function( Vinf )
        # spectral_library.GuralSpectral.resetAllElementalAbundances(self.spectral)

        # #========== Obtain model spectra for warm and hot temperature plasmas
        # # self.spectral.elemdata.els[kelem_Fe].user_fitflag = FITTING
        

        # spectral_library.GuralSpectral.elemFitting(self.spectral, self.elemIndex)

        spectral_library.GuralSpectral.computeWarmPlasmaSpectrum(self.spectral)
        spectral_library.GuralSpectral.computeHotPlasmaSpectrum(self.spectral)

        if self.HotTempOn_button.isChecked() and (self.WarmTempOn_button.isChecked() == False):
            print('Hot is on. Warm is not.')
            # spectral_library.GuralSpectral.computeHotPlasmaSpectrum(self.spectral)
            self.element_array = np.zeros((self.spectral.spcalib.nwavelengths,4))
            for i in range(self.spectral.spcalib.nwavelengths):
                self.element_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
                self.element_array[i][1] = 0
                self.element_array[i][2] = self.spectral.elemdata.els[self.elemIndex].spechi[i]
                self.element_array[i][3] = self.element_array[i][2]

        elif (self.HotTempOn_button.isChecked() == False) and self.WarmTempOn_button.isChecked():
            print('Warm is on. Hot is not.')
            # spectral_library.GuralSpectral.computeWarmPlasmaSpectrum(self.spectral)
            self.element_array = np.zeros((self.spectral.spcalib.nwavelengths,4))
            for i in range(self.spectral.spcalib.nwavelengths):
                self.element_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
                self.element_array[i][1] = self.spectral.elemdata.els[self.elemIndex].speclo[i]
                self.element_array[i][2] = 0
                self.element_array[i][3] = self.element_array[i][1]

        else:
            print('Warm and Hot are both on.')
            # spectral_library.GuralSpectral.computeWarmPlasmaSpectrum(self.spectral)
            # spectral_library.GuralSpectral.computeHotPlasmaSpectrum(self.spectral)
            self.element_array = np.zeros((self.spectral.spcalib.nwavelengths,4))
            for i in range(self.spectral.spcalib.nwavelengths):
                self.element_array[i][0] = self.spectral.spcalib.wavelength_nm[i]
                self.element_array[i][1] = self.spectral.elemdata.els[self.elemIndex].speclo[i]
                self.element_array[i][2] = self.spectral.elemdata.els[self.elemIndex].spechi[i]
                self.element_array[i][3] = self.element_array[i][1] + self.element_array[i][2]

        self.element_array[:,1] = self.element_array[:,1] * 10**self.Scale_rollbox.value()
        self.element_array[:,2] = self.element_array[:,2] * 10**self.Scale_rollbox.value()
        self.element_array[:,3] = self.element_array[:,3] * 10**self.Scale_rollbox.value()

    def resetAllElementalAbundances(self):
        spectral_library.GuralSpectral.resetAllElementalAbundances(self.spectral)

    def messageBox(self, message):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Warning)
        msgBox.setText(message)
        msgBox.setStandardButtons(QMessageBox.Ok)
        # msgBox.buttonClicked.connect(msgButtonClick)
        msgBox.exec()

    def elementButtonClicked(self):

        # Detect whether shift key is held down
        mods = QtGui.QApplication.keyboardModifiers()
        isShiftPressed = mods & QtCore.Qt.ShiftModifier

        # if bool(isShiftPressed) == True:
        #     self.sender().setStyleSheet('background-color:#FFFFFF;font:bold')
        #     print('Shift Pressed')

        self.buttonClicked = self.sender().objectName()
      
        self.elemName = str(self.buttonClicked).split('_')[0]
        print(self.elemName)
        self.buttonIndex = self.elementButtons.index(self.sender())
        self.elemNumber = self.elementDeets[self.buttonIndex][2]
        self.elemIndex = self.elementDeets[self.buttonIndex][3]
        print(self.elemNumber, self.elemIndex, self.buttonIndex)

        self.elementDeets[self.buttonIndex][1] += 1

        # print(self.spectral.elemdata.els[0:-1].element_string)
        # for i in range(self.spectral.elemdata.nelements):
        #     print(self.spectral.elemdata.els[i].element_string)

        if bool(isShiftPressed) == True and self.elementDeets[self.buttonIndex][1] != 1:
            if self.elementDeets[self.buttonIndex][1] == 2:
                self.elementsState[1] -= 1
            elif self.elementDeets[self.buttonIndex][1] == 3:
                self.elementsState[2] -= 1

            self.elementDeets[self.buttonIndex][1] = 0
            # self.spectral.elemdata.els[self.elemIndex].user_fitflag = 0
            self.spectral.elemdata.els[self.elemIndex].user_fitflag = 0
            
            spectral_library.GuralSpectral.removeElemFromModel(self.spectral, self.elemNumber)          
            self.statusBar.showMessage('Unlocked and %s removed from plot.' % self.elemName,2000)
            self.Plot.removeItem(globals()[self.elemName])
            globals()[self.elemName] = None
        
        if self.elementsState[1] > 2 and self.elementDeets[self.buttonIndex][1] == 2:
            self.sender().setStyleSheet('background-color:#FFFF00;color:#000000;')
            # spectral_library.GuralSpectral.lockElemFit(self.spectral, self.elemNumber)
            # self.spectral.elemdata.els[self.elemIndex].user_fitflag = 2
            self.spectral.elemdata.els[self.elemIndex].user_fitflag = 2
            print('Locked element and added to Fit')
            self.statusBar.showMessage('Locked %s fit' % self.elemName,2000)
            self.elementsState[2] += 1
            self.elementsState[1] -= 1
            # print(self.elementsState)
            return

        if self.elementsState[1] > 2:
            print('You can only have 3 elements in fitting mode at a time! Lock (left-click) or remove (shift-left-click) an element.')
            self.elementDeets[self.buttonIndex][1] -= 1
            self.messageBox('You can only have 3 elements in fitting mode at a time! Lock (left-click) or remove (shift-left-click) an element to proceed.')
            return

        
        if self.elementDeets[self.buttonIndex][1] == 1:
            self.sender().setStyleSheet('background-color:#00FF00;color:#000000;')
            spectral_library.GuralSpectral.elemFitting(self.spectral, self.elemIndex)
            print("elemNumber: %s" % self.elemNumber)
            print("elemIndex: %s" % self.elemIndex)
            # self.spectral.elemdata.els[self.elemIndex].user_fitflag = 1
            self.statusBar.showMessage('Set to fitting mode. %s' % self.elemName,2000)
            
            self.calculateElementSpectrum()
            self.plotElement(self)
            # print("N warm: %s" % self.spectral.elemdata.els[self.elemIndex].N_warm)
            self.columnDensityClicked()
            # print("N warm: %s" % self.spectral.elemdata.els[self.elemIndex].N_warm)

            # print('N warm: %s' % self.spectral.elemdata.els[elemIndex].N_warm)

            self.elementsState[1] += 1

            # Need to add calculated spectrum to the sum of the displayed spectra
            # and then (re)plot the summed spectrum

        elif self.elementDeets[self.buttonIndex][1] == 2:
            self.sender().setStyleSheet('background-color:#FFFF00;color:#000000;')
            # spectral_library.GuralSpectral.lockElemFit(self.spectral, self.elemNumber)
            self.spectral.elemdata.els[self.elemIndex].user_fitflag = 2
            print('Locked element and added to Fit')
            self.statusBar.showMessage('Locked %s fit' % self.elemName,2000)
            self.elementsState[2] += 1
            self.elementsState[1] -= 1

        elif self.elementDeets[self.buttonIndex][1] == 3:
            self.sender().setStyleSheet('background-color:#00FF00;color:#000000;')
            self.elementDeets[self.buttonIndex][1] = 1
            self.spectral.elemdata.els[self.elemIndex].user_fitflag = 1
            self.elementsState[1] += 1
            self.elementsState[2] -= 1

        else:
            self.sender().setStyleSheet('background-color:#FFFFFF;color:#000000;')

        # self.fitlessElems = []
        # self.fittingElems = []
        # self.lockedElems = []

        i = 0
        for i in range(len(self.elementDeets)):
            self.fitState[self.elementDeets[i][0]] = self.elementDeets[i][1]

        # print(self.elementsState)

    ################# DIRECT FILE CONTROL FUNCTIONS #################

    # Upload direct .vid file
    # def uploadDirectVid(self):
    #     """
    #     Uploads the direct video to the GUI.
    #     Uses readVid function from wmpl.
    #     Will display 0th frame first, unless script is changed manually.
    #     """
    
    #     # File path and name
    #     direct_path = self.DirectFilePath_linedit.text()
    #     print('Direct path is %s' % direct_path)
    #     direct_name = self.DirectFileName_linedit.text()
        
    #     # Read in the .vid file
    #     self.direct_vid = readVid(direct_path, direct_name)
    #     # Set frame number to be displayed
    #     self.direct_currentframe = 0
    #     # Define length of video in frames
    #     self.direct_vidlength = len(self.direct_vid.frames)

    #     # Update direct frame image view
    #     self.updateDirectFrames()

    def loadEventFile(self):
        self.spectral.vinfinity_kmsec = 40.0
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec():
            event_file_name = dlg.selectedFiles()
            print(event_file_name)

            if event_file_name[0].endswith('.txt'):
                file_path = os.path.split(event_file_name[0])[0]
                file_name = os.path.split(event_file_name[0])[1]

                if event_file_name == 'event.txt':            
                    event_name = file_path.split('/')[-1]
                else:
                    event_name = file_name.split('.')[0]


                with open(os.path.join(file_path,file_name)) as f:
                    lines = f.readlines()

                    print(lines[2].split())
                    print(lines[3].split())
                    print(lines[4].split())
                    
                    
                    if lines[3].split()[13] == "'KTJ'" or lines[3].split()[15] == "'KTJ'" or lines[4].split()[15] == "'KTJ'":
                        height = []
                        event_date = event_name.split('_')[0]
                        try:
                            self.spectral.vinfinity_kmsec = float(lines[4].split()[13])
                        except:
                            self.spectral.vinfinity_kmsec = 35.0
                        print('Vinf1: %f' % self.spectral.vinfinity_kmsec)
                        for line in lines:
                            if line.split()[0] == 'fit':
                                theta = float(line.split()[11])
                                mag = float(line.split()[23]) # Magnitude
                                try:
                                    height = float(line.split()[31])
                                except:
                                    height = 85.0
                                # print(theta, mag, height)

                        ### Now, use the event_name to open the vid file. Will need to decompress bz2, first
                        with bz2.open('/srv/meteor/klingon/events/' + event_date + '/ev_' + event_name + '_02I.vid.bz2') as f:
                            spectral_vid = f.read()

                        spectral_file_name = '/tmp/ev_%s_02I.vid' % event_name
                        with open(spectral_file_name, 'wb') as f:
                            f.write(spectral_vid)
                            print('Spectral file written to /tmp ...')

                        with bz2.open('/srv/meteor/klingon/events/' + event_date + '/ev_' + event_name + '_02J.vid.bz2') as f:
                            direct_vid = f.read()
                            # print(type(test_vid))

                        direct_file_name = '/tmp/ev_%s_02J.vid' % event_name
                        with open(direct_file_name, 'wb') as f:
                            f.write(direct_vid)
                            print('Direct file written to /tmp ...')


                        event_split = event_name.split('_')
                        
                        self.uploadDirectVid(direct_file_name)
                        # self.uploadSpectralVid(spectral_file_name)

                        os.remove(direct_file_name)
                        print('Direct file removed from /tmp ...')
                        os.remove(spectral_file_name)
                        print('Spectral file removed from /tmp ...')

                        print('Running auto-spectrum-find...')
                        self.findSpectrumFrame()

                    else:
                        print('This event has no spectral files.')

            else:
                pass

    def findSpectrumFrame(self):
        immean = []
        for i in range(self.spectral_vidlength):
            # Set frame 
            self.spectral_frame = self.spectral_vid.frames[i]
            self.spectral_frame_img = self.spectral_frame.img_data
            immean.append(np.max(np.mean(self.spectral_frame_img, axis=1)))

        maxindex = np.argmax(immean)
        # self.spectral_frame = self.spectral_vid.frames[9]
        # self.spectral_frame_img = self.spectral_frame.img_data
        # y=np.mean(self.spectral_frame_img, axis=1)

        # print(np.argmax(immean))
        # ysmooth = np.convolve(y, np.ones(10)/10, mode='valid')

        # plt.plot(immean)
        # plt.show()

        # Code below goes to best frame
        # self.direct_currentframe += 5
        # self.direct_currentframe = self.direct_currentframe%self.direct_vidlength
        
        # # Update frame shown in region
        # self.updateDirectFrames()

        return maxindex

    def uploadSpectralVid(self, file=None):
        if file == False:
            dlg = QFileDialog()
            dlg.setFileMode(QFileDialog.AnyFile)

            if dlg.exec():
                spectral_file_name = dlg.selectedFiles()

                if spectral_file_name[0].endswith('.vid'):
                    spectral_path = os.path.split(spectral_file_name[0])[0]
                    spectral_name = os.path.split(spectral_file_name[0])[1]

                    # self.DirectFilePath_linedit.setText(direct_path)
                    # self.DirectFileName_linedit.setText(direct_name)

                    self.SpectralFileName_label.setText(os.path.split(spectral_file_name[0])[1])

                    self.spectral_vid = readVid(spectral_path, spectral_name)
                    self.spectral_currentframe = int(len(self.spectral_vid.frames)/2)
                    self.spectral_vidlength = len(self.spectral_vid.frames)

                    self.updateSpectralFrames()

                    ### Enable some buttons ###
                    # self.FlattenSpectral_button.setEnabled(True)
                    self.AutoPick_button.setEnabled(True)
                    # self.AutoTrackSpectral_button.setEnabled(True)
                    self.SelectSpectralRegion_button.setEnabled(True)
                    self.ClearSpectralRegion_button.setEnabled(True)
                    # self.CheckSpectralRegion_button.setEnabled(True)
                    # self.CheckSpectralBackground_button.setEnabled(True)
                    # self.UploadSpectralBias_button.setEnabled(True)
                    self.UploadSpectralFlat_button.setEnabled(True)
                    self.AutoSpectralFlat_button.setEnabled(True)
                    # self.Affine_button.setEnabled(True)
                    # self.PickFeature_button.setEnabled(True)
                    self.AutoPick_button.setEnabled(True)
                    # self.AutoTrackSpectral_button.setEnabled(True)
                    # self.ManualPickDirect_button.setEnabled(True)
                    # self.ClearPicksSpectral_button.setEnabled(True)
                    self.DeltaX_edit.setEnabled(True)
                    self.DeltaY_edit.setEnabled(True)
                    self.nm0_edit.setEnabled(True)
                    # self.UpdateAffine_button.setEnabled(True)
                                        
                else:
                    pass

        else:
            spectral_file_name = file

            if spectral_file_name.endswith('.vid'):
                spectral_path = os.path.split(spectral_file_name)[0]
                spectral_name = os.path.split(spectral_file_name)[1]

                # self.DirectFilePath_linedit.setText(direct_path)
                # self.DirectFileName_linedit.setText(direct_name)

                self.SpectralFileName_label.setText('Spectral camera file: ' + os.path.split(spectral_file_name)[1])

                self.spectral_vid = readVid(spectral_path, spectral_name)
                self.spectral_currentframe = int(len(self.spectral_vid.frames)/2)
                self.spectral_vidlength = len(self.spectral_vid.frames)

                self.updateSpectralFrames()

                ### Enable some buttons ###
                self.FlattenSpectral_button.setEnabled(True)
                self.AutoPick_button.setEnabled(True)
                # self.AutoTrackSpectral_button.setEnabled(True)
                self.SelectSpectralRegion_button.setEnabled(True)
                self.ClearSpectralRegion_button.setEnabled(True)
                self.CheckSpectralRegion_button.setEnabled(True)
                self.CheckSpectralBackground_button.setEnabled(True)
                self.UploadSpectralBias_button.setEnabled(True)
                self.UploadSpectralFlat_button.setEnabled(True)
                self.AutoSpectralFlat_button.setEnabled(True)
                # self.Affine_button.setEnabled(True)
                # self.PickFeature_button.setEnabled(True)
                self.AutoPick_button.setEnabled(True)
                # self.AutoTrackSpectral_button.setEnabled(True)
                # self.ManualPickDirect_button.setEnabled(True)
                # self.ClearPicksSpectral_button.setEnabled(True)


    # # Update direct frame
    # def updateDirectFrames(self):
    #     """
    #     Updates frame shown in the direct file image view. 
    #     Updates time stamp and frame number. 
    #     Adjusts image levels.
    #     """

    #     # Set frame to be displayed
    #     self.direct_frame = self.direct_vid.frames[self.direct_currentframe]
    #     self.direct_frame_img = self.direct_frame.img_data

    #     # Display time
    #     self.dt = unixTime2Date(self.direct_frame.ts, self.direct_frame.tu, dt_obj=False)

    #     # self.dt = str(self.dt)
    #     dt_h = str(self.dt[3])
    #     dt_m = str(self.dt[4])
    #     dt_s = str(self.dt[5])
    #     dt_us = str(round(self.dt[6]*1000))

    #     if len(dt_us) == 5:
    #         dt_us = '0' + dt_us
    #     elif len(dt_us) == 4:
    #         dt_us = '00' + dt_us

    #     self.DirectTime_label.setText(' at ' + dt_h + ':' + dt_m + \
    #      ':' + dt_s + '.' + dt_us + 'UT on ' + str(self.dt[0]) + '/' + \
    #      str(self.dt[1]) + '/' + str(self.dt[2]))
    #     self.update()

    #     # Display frame number
    #     # self.DirectFrame_label.setNum(self.direct_currentframe)

    #     self.DirectFrame_label.setText('Frame # ' + str(self.direct_currentframe))
    #     # self.DirectFrame_label.setText('Viewing frame #' + self.direct_currentframe)
    #     self.update()

    #     # Set image levels
    #     minv = np.percentile(self.direct_frame_img, 0.1)
    #     maxv = np.percentile(self.direct_frame_img, 99.95)
    #     gamma = 1

    #     # Create an image with properly adjusted levels
    #     self.direct_frame_img = adjustLevels(self.direct_frame_img, minv, gamma, maxv, scaleto8bits=True)
        
    #     # Set the image to be displayed
    #     self.direct_image.setImage(self.direct_frame_img.T)


    # Move to next Direct frame
    # def nextSpectralFrame(self):
    #     """
    #     Increases the direct frame number by 1 frame to show the next frame.
    #     """

    #     # Increase frame number by one
    #     self.direct_currentframe += 1
    #     self.direct_currentframe = self.direct_currentframe%self.direct_vidlength
        
    #     # Update frame shown in region
    #     self.updateDirectFrames()

    # def forwardFiveSpectralFrames(self):
    #     """
    #     Increases the direct frame number by 5 frames.
    #     """

    #     # Increase frame number by five
    #     self.direct_currentframe += 5
    #     self.direct_currentframe = self.direct_currentframe%self.direct_vidlength
        
    #     # Update frame shown in region
    #     self.updateDirectFrames()

    # # Move to last Direct frame
    # def lastSpectralFrame(self):
    #     """
    #     Decrease the direct frame number by 1 frame to show the previous frame.
    #     """

    #     # Decrease frame number by one
    #     self.direct_currentframe -= 1
    #     self.direct_currentframe = self.direct_currentframe%self.direct_vidlength
        
    #     # Update frame shown in region
    #     self.updateDirectFrames()

    # def backFiveSpectralFrames(self):
    #     """
    #     Decrease the direct frame number by 5 frames.
    #     """

    #     # Decrease frame number by five
    #     self.direct_currentframe -= 5
    #     self.direct_currentframe = self.direct_currentframe%self.direct_vidlength
        
    #     # Update frame shown in region
    #     self.updateDirectFrames()

    def rotateVid(self):
        print('Rotating vid file...')
        print('Rotation: %s' % self.SpectralRotation_rollbox.value())


    # Get position of the mouse
    def getSpectralPosition(self, event):
        """
        Will display the coordinates of a mouse click in the direct image view.
        Information will be updated each time the mouse is clicked within the direct image.
        """
        # # Detect whether shift key is held down
        # mods = QtGui.QApplication.keyboardModifiers()
        # isShiftPressed = mods & QtCore.Qt.ShiftModifier

        # if bool(isShiftPressed) == True:
        # #     self.sender().setStyleSheet('background-color:#FFFFFF;font:bold')
        #     print('Shift Pressed')

        self.updateSpectralROI()

        # if event.button() == QtCore.Qt.RightButton:
        # print('Right-click, setting index wavelength.')
        # Set x and y coordinates of click
        self.dir_xpos = event.pos().x()
        self.dir_ypos = event.pos().y()

        # Update label with new click coordinates
        self.SpectralXYCoordsDisplay_label.setText('Mouse coords: ( %d : %d )' % (self.dir_xpos, self.dir_ypos))

        self.hu = self.dir_xpos
        self.hv = self.dir_ypos

        self.DeltaX_edit.setText(str(int(self.hu - self.dir_x)))
        self.DeltaY_edit.setText(str(int(self.hv - self.dir_y)))
        print(self.hu-self.dir_x)

        self.affine_markers.setData(x = [self.hu], y = [self.hv])

        # self.PickFeature_button.setChecked(False)
        return

    # Run the affine transform, plot point on spectral image
    def affineTransform(self):
        """
        Runs affine transform from a given point on the direct image
        to a corresponding point on the spectral image, with information 
        from getDirectPosition. 

        Sets data for the affine_markers, and displays the point on the 
        spectral image.
        """
        dlg = QFileDialog(filter="Affine files (*.aff)")
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec():
            scale_file_name = dlg.selectedFiles()
            # print(scale_file_name)
        
        # Scale plate file path
        scale_dir_path = "."
        # scale_file_name = "direct_spectral_20210526_02J.aff"

        # Load the scale plate
        # self.scale = loadScale(scale_dir_path, scale_file_name[0])

        # Enable the update button
        # self.UpdateAffine_button.setEnabled(True)

        # Convert image (X, Y) to encoder (Hu, Hv)
        # self.hu, self.hv = plateScaleMap(self.scale, self.dir_x, self.dir_y)

        self.hu = self.dir_x + 300
        self.hv = self.dir_y

        # Set data for markers and plot on 
        self.affine_markers.setData(x = [self.hu], y = [self.hv])

    def updateTransform(self):
        self.nm0_label.setText(self.nm0_edit.text())
        deltaX = int(self.DeltaX_edit.text())
        deltaY = int(self.DeltaY_edit.text())
        # self.hu, self.hv = plateScaleMap(self.scale, self.dir_x + deltaX, self.dir_y + deltaY)
        self.hu = self.dir_x + deltaX
        self.hv = self.dir_y + deltaY
        self.affine_markers.setData(x = [self.hu], y = [self.hv])

    ################# SPECTRAL FILE CONTROL FUNCTIONS #################

    # Upload spectral .vid file
    # def uploadSpectralVid(self):
    #     """
    #     Uploads the direct video to the GUI.
    #     Uses readVid function from wmpl.
    #     Will display 0th frame first, unless script is changed manually.
    #     """

    #     # File path and name
    #     # spectral_path = self.SpectralFilePath_linedit.text()
    #     # spectral_name = self.SpectralFileName_linedit.text()
        
    #     # Read in the .vid file
    #     self.spectral_vid = readVid(spectral_path, spectral_name)

    #     # Initialize the current frame
    #     self.spectral_currentframe = 0
    #     self.spectral_vidlength = len(self.spectral_vid.frames)

    #     # Update display
    #     self.updateSpectralFrames()

    def uploadSpectralFlat(self, dtype=None, byteswap=True, dark=None):

        dlg = QFileDialog(filter='PNG files (*.png)')
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec():
            flat_file_name = dlg.selectedFiles()
            
            if flat_file_name[0].endswith('.png'):
                # print(flat_file_name[0])
               
                # Load the flat image
                flat_img = loadImage(flat_file_name[0], -1)

                # # Change the file type if given
                # if dtype is not None:
                #     flat_img = flat_img.astype(dtype)
                #  # If the flat isn't a 8 bit integer, convert it to uint16
                # elif flat_img.dtype != np.uint8:
                #     flat_img = flat_img.astype(np.uint16)

                if flat_img.dtype != np.uint8:
                    flat_img = flat_img.astype(np.uint16)

                if byteswap:
                    flat_img = flat_img.byteswap()
            
                # plt.figure()
                # plt.imshow(flat_img)
                # plt.show()

                # Init a new Flat structure
                self.flat_structure = FlatStruct(flat_img, dark=dark)
                self.spectral_flat_exists = True
            else:
                pass

    def autoSpectralFlat(self):

        flats = [self.spectral_vid.frames[i].img_data for i in range(len(self.spectral_vid.frames))]
        flat_stack = np.stack(flats)

        # auto_flat = np.median(flat_stack, axis=0)
        auto_flat = np.ones(np.shape(self.spectral_vid.frames[0].img_data))

        self.flat_structure = FlatStruct(auto_flat, dark=None)

        self.spectral_frame_img = applyFlat(self.spectral_frame_img, self.flat_structure)

        # Set image levels
        self.minv = np.percentile(self.spectral_frame_img, 0.05)
        self.maxv = np.percentile(self.spectral_frame_img, 99.95)
        gamma = 1

        # Create an image with properly adjust levels
        spectral_frame_img = adjustLevels(self.spectral_frame_img, self.minv, gamma, self.maxv, scaleto8bits=True)

        # Set spectral image
        self.spectral_image.setImage(spectral_frame_img.T)
        self.spectral_flat_exists = True

    # def uploadSpectralVid(self, file=None):
    #     if file == False:
    #         dlg = QFileDialog()
    #         dlg.setFileMode(QFileDialog.AnyFile)

    #         if dlg.exec():
    #             spectral_file_name = dlg.selectedFiles()

    #             if spectral_file_name[0].endswith('.vid'):
    #                 # print(spectral_file_name)
    #                 spectral_path = os.path.split(spectral_file_name[0])[0]
    #                 spectral_name = os.path.split(spectral_file_name[0])[1]

    #                 self.SpectralFileName_label.setText('Spectral camera file: ' + os.path.split(spectral_file_name[0])[1])

    #                 self.spectral_vid = readVid(spectral_path, spectral_name)

    #                 self.spectral_vidlength = len(self.spectral_vid.frames)
    #                 self.spectral_currentframe = findSpectrumFrame()
    #                 # print(self.spectral_currentframe)
    #                 # self.spectral_currentframe = int(len(self.spectral_vid.frames)/2)
                    

    #                 self.updateSpectralFrames()

    #                 ### Enable some buttons ###
    #                 self.FlattenSpectral_button.setEnabled(True)
    #                 self.AutoPick_button.setEnabled(True)
    #                 self.AutoTrackSpectral_button.setEnabled(True)
    #                 self.SelectSpectralRegion_button.setEnabled(True)
    #                 self.ClearSpectralRegion_button.setEnabled(True)
    #                 self.CheckSpectralRegion_button.setEnabled(True)
    #                 self.CheckSpectralBackground_button.setEnabled(True)
    #                 self.UploadSpectralBias_button.setEnabled(True)
    #                 self.UploadSpectralFlat_button.setEnabled(True)
    #                 self.AutoSpectralFlat_button.setEnabled(True)
    #                 self.Affine_button.setEnabled(True)
    #             else:
    #                 pass
    #     else:
    #         spectral_file_name = file

    #         if spectral_file_name.endswith('.vid'):
    #             spectral_path = os.path.split(spectral_file_name)[0]
    #             spectral_name = os.path.split(spectral_file_name)[1]

    #             # self.DirectFilePath_linedit.setText(direct_path)
    #             # self.DirectFileName_linedit.setText(direct_name)

    #             self.SpectralFileName_label.setText('Spectral camera file: ' + os.path.split(spectral_file_name)[1])

    #             self.spectral_vid = readVid(spectral_path, spectral_name)
    #             self.spectral_vidlength = len(self.spectral_vid.frames)
    #             self.spectral_currentframe = self.findSpectrumFrame()
    #             # self.spectral_currentframe = int(len(self.spectral_vid.frames)/2)
    #             # self.spectral_vidlength = len(self.spectral_vid.frames)

    #             self.updateSpectralFrames()

    #             self.FlattenSpectral_button.setEnabled(True)
    #             self.AutoPick_button.setEnabled(True)
    #             self.AutoTrackSpectral_button.setEnabled(True)
    #             self.SelectSpectralRegion_button.setEnabled(True)
    #             self.ClearSpectralRegion_button.setEnabled(True)
    #             self.CheckSpectralRegion_button.setEnabled(True)
    #             self.CheckSpectralBackground_button.setEnabled(True)
    #             self.UploadSpectralBias_button.setEnabled(True)
    #             self.UploadSpectralFlat_button.setEnabled(True)
    #             self.AutoSpectralFlat_button.setEnabled(True)

    # Update spectral frame
    def updateSpectralFrames(self):
        """
        Updates frame shown in the spectral file image view. 
        Updates time stamp and frame number. 
        Applies image flat.
        """

        # Set frame 
        self.spectral_frame = self.spectral_vid.frames[self.spectral_currentframe]

        self.spectral_frame_img_og = self.spectral_frame.img_data
        
        # self.spectral_frame_img = self.avg_frame_img

        self.spectral_frame_img = scipy.ndimage.rotate(self.spectral_frame.img_data, self.SpectralRotation_rollbox.value())
        # print('Spectral frame shape: %s x %s' % (np.shape(self.spectral_frame_img)[0], np.shape(self.spectral_frame_img)[1]))

        try:
            if self.spectral_flat_exists:
                self.spectral_frame_img = applyFlat(self.spectral_frame_img, self.flat_structure)
        except:
            self.spectral_flat_exists = False
            pass

        # Set image levels
        if self.spectral_flat_exists == False:
            self.minv = np.percentile(self.spectral_frame_img_og, 0.05)
            self.maxv = np.percentile(self.spectral_frame_img_og, 99.95)
        
        gamma = 1

        # Create an image with properly adjust levels
        spectral_frame_img = adjustLevels(self.spectral_frame_img, self.minv, gamma, self.maxv, scaleto8bits=True)

        # Display time
        self.st = unixTime2Date(self.spectral_frame.ts, self.spectral_frame.tu, dt_obj=False)
        # print('DT-spec')
        # print(self.st)
        # self.st = str(self.st)
        # self.SpectralTime_label.setText(self.st)

        st_h = str(self.st[3])
        st_m = str(self.st[4])
        st_s = str(self.st[5])
        st_us = str(round(self.st[6]*1000))

        if len(st_us) == 5:
            st_us = '0' + st_us
        elif len(st_us) == 4:
            st_us = '00' + st_us

        self.SpectralTime_label.setText(' at ' + st_h + ':' + st_m + \
         ':' + st_s + '.' + st_us + 'UT on ' + str(self.st[0]) + '/' + \
         str(self.st[1]) + '/' + str(self.st[2]))
        self.update()

        # Display frame number
        self.SpectralFrame_label.setText('Frame # ' + str(self.spectral_currentframe))
        self.update()

        # Call the function to apply the flat
        # if self.flat_structure is not None: 
        #     self.spectral_frame_img = applyFlat(self.spectral_frame_img, self.flat_structure)

        # Set spectral image
        self.spectral_image.setImage(spectral_frame_img.T)

    def updateSpectralROI(self):
        # print('Updating Spectral ROI image')
        self.checkSpectralRegion()
        self.spectralROI_image.setImage(self.spectral_array)

    # Click to next Spectral frame
    def nextSpectralFrame(self):
        """
        Increases the spectral frame number by 1 to show the next frame
        """

        # Increase spectral frame number by 1
        self.spectral_currentframe += 1
        self.spectral_currentframe = self.spectral_currentframe%self.spectral_vidlength
        
        # Update frame shown in region
        self.updateSpectralFrames() 

    def forwardFiveSpectralFrames(self):
        """
        Increases the spectral frame number by 1 to show the next frame
        """

        # Increase spectral frame number by 1
        self.spectral_currentframe += 5
        self.spectral_currentframe = self.spectral_currentframe%self.spectral_vidlength
        
        # Update frame shown in region
        self.updateSpectralFrames() 

    # Click to last Spectral frame
    def lastSpectralFrame(self):
        """
        Decrease the spectral frame number by 1 to show the previous frame.
        """

        # Decrease frame number by 1
        self.spectral_currentframe -= 1
        self.spectral_currentframe = self.spectral_currentframe%self.spectral_vidlength
        
        # Update frame shown in region
        self.updateSpectralFrames()

    def backFiveSpectralFrames(self):
        """
        Decrease the spectral frame number by 1 to show the previous frame.
        """

        # Decrease frame number by 1
        self.spectral_currentframe -= 5
        self.spectral_currentframe = self.spectral_currentframe%self.spectral_vidlength
        
        # Update frame shown in region
        self.updateSpectralFrames()

    # Introduces ROI box
    def spectralROI(self):
        """
        Introduces square ROI box on spectral frame. 
        Includes handles to rotate, scale, or translate the ROI. 
        """

        # Introduce ROI, provided ROI is None 
        if self.spectral_roi is None:
            self.spectral_roi = pg.ROI((0,0), size = (500, 20), angle = 0, invertible = False, \
                maxBounds = None, snapSize = 1, scaleSnap = False, \
                    translateSnap = False, rotateSnap = False, \
                        parent = self.spectral_image, \
                            pen = None, movable = True, \
                            rotatable = True, resizable = True, removable = True)
            self.spectral_roi.addRotateHandle([0.5,0.5], [0.25, 0.25])
            self.spectral_roi.addScaleHandle([1,0.5], [0,0])
            self.spectral_roi.addTranslateHandle([0,0.5],  [0,0])
            self.angle = self.spectral_roi.angle()
            self.updateSpectralROI()


    def spectralAutoROI(self,width, height, roll, intercept):
        """
        Introduces square ROI box on spectral frame. 
        Includes handles to rotate, scale, or translate the ROI. 
        """
        if self.spectral_roi is not None:
            self.spectral_roi.deleteLater()
            self.spectral_roi = None


        # Introduce ROI, provided ROI is None 
        if self.spectral_roi is None:
            self.spectral_roi = pg.ROI((0,intercept-1/2*height), size = (width, height), angle = roll, invertible = False, \
                maxBounds = None, snapSize = 1, scaleSnap = False, \
                    translateSnap = False, rotateSnap = False, \
                        parent = self.spectral_image, \
                            pen = None, movable = True, \
                            rotatable = True, resizable = True, removable = True)
            self.spectral_roi.addRotateHandle([0.5,0.5], [0.25, 0.25])
            self.spectral_roi.addScaleHandle([1,0.5], [0,0])
            self.spectral_roi.addTranslateHandle([0,0.5],  [0,0])
            self.angle = self.spectral_roi.angle()
            self.updateSpectralROI()
    
    def autoTrackDirect(self):
        
        # Set frame to be displayed
        # print(self.direct_currentframe)
        # self.direct_frame = self.direct_vid.frames[self.direct_currentframe]
        # self.direct_frame_img = self.direct_frame.img_data

        dimage = self.direct_vid.frames[23].img_data

        plt.imshow(dimage)
        plt.show()


    def autoPickDirect(self):

        dimage = self.spectral_frame_img

        # Calculate the global mean and stddev
        global_mean = np.mean(dimage)
        global_stddev = np.std(dimage)
        
        # Change data type to 32 bit
        data = dimage.astype(np.float32)
        # data = dimage

        # Apply a mean filter
        fx = 3
        fy = 3
        data = scipy.ndimage.filters.convolve(data, weights=np.full((fx,fy), 1.0/4))

        # Locate local maxima
        neighborhood_size = 80
        intensity_threshold = 400
        data_max = scipy.ndimage.filters.maximum_filter(data, neighborhood_size)

        maxima = (data == data_max)
        data_min = scipy.ndimage.filters.minimum_filter(data, neighborhood_size)
        diff = ((data_max - data_min) > intensity_threshold)

        maxima[diff == 0] = 0

        # Find and label the maxima
        labeled, num_objects = scipy.ndimage.label(maxima)

        # Find centres of mass
        xy = np.array(scipy.ndimage.center_of_mass(data, labeled, range(1, num_objects+1)))

        # Unpack coordinates
        y, x = np.hsplit(xy,2)

        x2, y2, amplitude, intensity, sigma_y_fitted, sigma_x_fitted = fitPSF(dimage, global_mean, x, y)

        extractions = np.array(list(zip(x2, y2, amplitude, intensity, sigma_x_fitted, sigma_y_fitted)), dtype=object)

        # # Delete the direct ROI
        # self.direct_roi.deleteLater()

        # # Re-initialize the direct ROI
        # self.direct_roi = None

        self.dir_x = extractions[0,0]
        self.dir_y = extractions[0,1]

        self.spectral_markers.setData(x = [self.dir_x], y = [self.dir_y])
        self.spectral_circle.setData(x = [self.dir_x], y = [self.dir_y])

        # print(self.dir_x, self.dir_y)

        # for i in range(-25,25):
        #     for j in range(-25,25):
        #         self.direct_vid.frames[self.direct_currentframe].img_data.T[int(self.dir_x)+i,int(self.dir_y)+j] = 30

        # print(np.amin(self.direct_vid.frames[self.direct_currentframe].img_data))

        # self.updateDirectFrames()


    def autoPickROI(self):

       # image = imageio.imread('TestSpectrum1.png', as_gray=True)
        image = self.spectral_frame_img

        # fig = plt.figure(figsize=(10,6))
        # ax = fig.add_subplot(111)

        y,x = np.indices(image.shape)

        mae = []
        scores = []
        roll_ransac = []
        m_ransac = []
        b_ransac = []

        # print(np.amax(image))

        for i in range(1,np.amax(image),int(np.amax(image)*0.02)):
            # print(i)
            # axins = ax.inset_axes([-0.1,0.8,0.3,0.3])
            # ax.set_xlim(0,image.shape[1])
            # ax.set_ylim(0,image.shape[0])
            # plt.gca().invert_yaxis()

            valid_z = (y.ravel()>0) & (image.ravel()>(50+i))
            x_valid = x.ravel()[valid_z]
            y_valid = y.ravel()[valid_z]
            z_valid = image.ravel()[valid_z]

            if len(z_valid) > 200:
                # sc1 = ax.scatter(x_valid,y_valid, color='yellowgreen', marker='.')
                # sc2 = ax.scatter(x_valid,y_valid, color='gold', marker='.')

                ransac = RANSACRegressor(residual_threshold=5).fit(x_valid.reshape(-1,1), y_valid.reshape(-1,1), sample_weight=z_valid)
                # ransac.fit(x_valid.reshape(-1,1), y_valid.reshape(-1,1), sample_weight=z_valid**2)
                inlier_mask = ransac.inlier_mask_
                outlier_mask = np.logical_not(inlier_mask)

                line_X = np.arange(x_valid.min(), x_valid.max())[:,np.newaxis]
                line_y_ransac = ransac.predict(line_X)

                prediction = ransac.predict(x_valid.reshape(-1,1))

                mae.append(mean_absolute_error(y_valid,prediction))
                scores.append(ransac.score(x_valid.reshape(-1,1), y_valid.reshape(-1,1)))

                # sc1.set_offsets(np.c_[x_valid[inlier_mask], y_valid[inlier_mask]])
                # sc2.set_offsets(np.c_[x_valid[outlier_mask], y_valid[outlier_mask]])

                z = np.polyfit(x_valid,y_valid, w=z_valid, deg=1)
                p = np.poly1d(z)

                x_plot = np.linspace(x_valid.min(), x_valid.max(), 100)
                y_plot = p(x_plot)

                # ax.plot(x_plot, y_plot, '-r', lw=2, label='LR')
                # ax.plot(line_X, line_y_ransac, label='RANSAC')

                # axins.plot(mae, label='MAE')

                # ax.legend(loc='upper right')
                # axins.legend(loc='upper right')

                this_roll = -1*math.degrees(math.atan2((line_y_ransac.max()-line_y_ransac.min()),(line_X.max()-line_X.min())))                
                roll_ransac.append(this_roll)

                this_b = line_y_ransac.max()-((line_y_ransac.max()-line_y_ransac.min())/(line_X.max()-line_X.min()))*line_X.max()
                b_ransac.append(this_b)
                
                # ax.text(int(image.shape[1]/4),-10,f'This roll = {this_roll:.3} degrees')

                if len(mae) > 1:
                    best_roll = roll_ransac[np.argmin(mae)]
                    # ax.text(int(image.shape[1]/2),-10,f'Best roll = {best_roll:.3} degrees')

                    self.Roll_rollbox.setValue(best_roll)

                    # x_vals = np.linspace(0,image.shape[1],2)
                    # y_vals = b_ransac[np.argmin(mae)]+20+ np.tan(math.radians(best_roll))*x_vals
                    # ax.plot(x_vals, y_vals, '--')

                    if self.spectral_roi is not None:
                        self.spectral_roi.deleteLater()
                        self.spectral_roi = None

                    self.spectralAutoROI(image.shape[1],10,best_roll,b_ransac[np.argmin(mae)])

            else:
                break
            # # print(line_y_ransac.max()-((line_y_ransac.max()-line_y_ransac.min())/(line_X.max()-line_X.min()))*line_X.max())
            # plt.axis('equal')

            # fig.canvas.draw_idle()
            # plt.pause(0.01)
            # ax.cla()

            # plt.show()

        # print(self.dir_x, self.dir_y)
        # plt.waitforbuttonpress()

        # self.spectralAutoROI(image.shape[1],20,best_roll,b_ransac[np.argmin(mae)])

    # Check image background
    def checkSpectralBackground(self):
        """
        Uses frame sets 1 and 2 to build a spectral background.
        This background will be subtracted from the measured spectrum to reduce noise. 
        Background calculation will include image flat, if one exists. 
        """
        
        # First frame set start
        # self.spectral_background_startframe_beg = int(self.Frame1SpectralStart_linedit.text())
        # # First frame set end
        # self.spectral_background_startframe_end = int(self.Frame1SpectralEnd_linedit.text())
        # # Last frame set start
        # self.spectral_background_lastframe_beg = int(self.Frame2SpectralStart_linedit.text())
        # # Last frame set end
        # self.spectral_background_lastframe_end = int(self.Frame2SpectralEnd_linedit.text())

        self.spectral_background_startframe_beg = 10
        # First frame set end
        self.spectral_background_startframe_end = 20
        # Last frame set start
        self.spectral_background_lastframe_beg = 35
        # Last frame set end
        self.spectral_background_lastframe_end = 45


        # Define frame range
        frame_range = list(range(self.spectral_background_startframe_beg, self.spectral_background_startframe_end))
        frame_range += list(range(self.spectral_background_lastframe_beg, self.spectral_background_lastframe_end))
        
        # Build array 
        self.spectral_background = np.zeros(shape=(len(frame_range), \
            self.spectral_vid.frames[0].img_data.shape[0], self.spectral_vid.frames[0].img_data.shape[1]))

        # Fill array
        for k, i in enumerate(frame_range):
            frame = self.spectral_vid.frames[i].img_data
            self.spectral_background[k] = frame

        # Take median value of each entry, convert to float
        self.spectral_background = np.median(self.spectral_background, axis = 0)
        self.spectral_background = scipy.ndimage.median_filter(self.spectral_background, size = 10)

        # Apply to flat, if flat exists
        if self.flat_structure is not None: 
            self.spectral_background = applyFlat(self.spectral_background.astype(np.uint16), self.flat_structure) 
        
        # Define spectral background array, each entry a float
        self.spectral_background = self.spectral_background.astype(np.float64)

        # plt.figure()
        # plt.imshow(self.spectral_background)
        # plt.show()

    # Makes spectral background appear in a new window
    def showSpectralBackground(self):
        """
        Displays spectral background in a pop-up window for the user. 
        If the user is not satisfied with the background, they can enter
        new frame set numbers and check the background again. 

        This function call is not necessary before plotting. 
        """

        # Calculate spectral background
        self.checkSpectralBackground()

        # Display background in a pop-up window
        plt.imshow(self.spectral_background, cmap = 'gray')
        plt.show()

    # Check selected region
    def checkSpectralRegion(self):
        """
        Get the array region in the spectralROI to define the image data.
        """

        # Calculate image background
        self.checkSpectralBackground()

        # Re define image data
        spectral_frame = self.spectral_vid.frames[self.spectral_currentframe]
        spectral_frame_img = spectral_frame.img_data.astype(np.float64) - self.spectral_background
        spectral_frame_img[spectral_frame_img < 0] = 0
        spectral_frame_img = self.spectral_frame_img.astype(np.uint16)
        # print(np.mean(self.spectral_background))
        # spectral_frame_img = spectral_frame.img_data
        # spectral_frame_img = self.spectral_frame_img.astype(np.uint16)

        # Get array region from ROI
        self.spectral_array = self.spectral_roi.getArrayRegion(spectral_frame_img.T, self.spectral_image)

    # Makes ROI appear in a new window
    def showSpectralRegion(self):
        """"
        Displays spectral ROI in a pop-up window for the user. 
        If the user is not satisfied, they can re-adjust the size/angle 
        of the ROI and check again. 

        This function call is not necessary before plotting.
        """

        # Calculate spectral region
        self.checkSpectralRegion()        
        
        # Display region in a pop-up window
        # plt.imshow(self.spectral_array.T,cmap = 'gray')
        # plt.show()

    # Clear ROI box
    def clearSpectralROI(self):
        """
        Clears spectral ROI box from spectral image view. 
        Deletes all associated data. 
        Re-initializes spectral_roi to None, so it can be 
        called again later. 

        This function means only one ROI box can be present
        on the image view at a time.
        """

        # Delete the spectral ROI
        self.spectral_roi.deleteLater()

        # Re-initialize the spectral ROI
        self.spectral_roi = None

    # Clear the affine marker - ***
    def clearAffine(self):
        """
        Clears the affine marker from the spectral image view. 
        Re-initializes the marker so the affine transform function
        can be performed again, with a new marker. 
        """

        # Clear the marker
        self.spectralROI_imageframe.removeItem(self.affine_markers)

        # Re-initalize
        self.affine_markers = pg.ScatterPlotItem()
        self.affine_markers.setPen('r')
        self.affine_markers.setSymbol('o')
        self.spectralROI_imageframe.addItem(self.affine_markers)   
    
    # Load spectral flat loadFlat(dir_path, file_name, dtype=None, byteswap=True, dark=None


    # Remove image flat - ***
    def removeSpectralFlat(self):
        """
        Allows user to remove flat from the current frame only. 
        Flat can be restored by moving to next/last frame, and then 
        returning to desired frame. 

        Background and spectrum will still be calculted with flat applied,
        the flat can only be removed for user convenience.
        """
        # Set frame 
        self.spectral_frame = self.spectral_vid.frames[self.spectral_currentframe]
        self.spectral_frame_img = self.spectral_frame.img_data

        # Display time
        self.st = unixTime2Date(self.spectral_frame.ts, self.spectral_frame.tu, dt_obj=False)
        self.st = str(self.st)
        self.SpectralTime_label.setText(self.st)
        self.update()

        # Display frame number
        self.SpectralFrame_label.setNum(self.spectral_currentframe)
        self.update()

        # Set image levels
        minv = np.percentile(self.spectral_frame_img, 0.05)
        maxv = np.percentile(self.spectral_frame_img, 99.95)
        gamma = 1

        # Create an image with properly adjust levels
        spectral_frame_img = adjustLevels(self.spectral_frame_img, minv, gamma, maxv, scaleto8bits=True)

        # Set spectral image
        self.spectral_image.setImage(spectral_frame_img.T)

    ################# JOINT FILE CONTROL FUNCTIONS #################

    # # Next frames
    # def nextFrame(self):
    #     """
    #     Moves both spectral and direct frames forward by 1 frame.
    #     """
    #     # Direct frame
    #     self.nextDirectFrame()

    #     # Spectral frame
    #     self.nextSpectralFrame()

    # # Last frames
    # def lastFrame(self):
    #     """
    #     Moves both spectral and direct frames back by 1 frame.
    #     """
    #     # Direct frame
    #     self.lastDirectFrame()

    #     # Spectral Frame
    #     self.lastSpectralFrame()

    # # Next frames that are as close in time as possible
    # def nextTimeFrame(self):
    #     """
    #     Finds next frame set that is as close in time as possible.
    #     Evaluates tenths and hundredths of a second
    #     portion of timestamp, finds the difference
    #     between the direct and spectral timestamps, then recognizes the time 
    #     gap and finds the appropriate method by which to get the next closest
    #     frames. 

    #     Does not work if the time gap is too great - user must have the 
    #     timestamp within 100 milliseconds.
    #     """

    #     # Tenths place value of seconds for the direct timestamp
    #     self.dir_ms_a = int(self.dt[25])
    #     # Hundredths place value of seconds for the direct timestamp
    #     self.dir_ms_b = int(self.dt[26])
    #     # Tenths place value of the seconds for the spectral timestamp
    #     self.spec_ms_a = int(self.st[25])
    #     # Hundredths place value of seconds for the spectral timestamp
    #     self.spec_ms_b = int(self.st[26])

    #     # Subtract tenths value of direct timestamp from tenths value of spectral timestamp
    #     self.a = self.dir_ms_a - self.spec_ms_a
    #     # Repeat above for hundredths
    #     self.b = self.dir_ms_b - self.spec_ms_b

    #     # If the tenths values are the same
    #     if self.a == 0:
    #         if 0 < self.b < 5:
    #             self.nextFrame()

    #         if 5 <= self.b < 9:
    #             self.nextSpectralFrame()
            
    #         if 0 >= self.b > -5:
    #             self.nextSpectralFrame()
    #             self.nextDirectFrame()
            
    #         if -5 >= self.b > -10:
    #             self.nextDirectFrame()
        
    #     # If the direct timestamp is ahead of the spectral timestamp
    #     if self.a ==1:

    #         if self.b == 0:
    #             self.nextSpectralFrame()
    #             self.nextSpectralFrame()

    #         if 0 < self.b < 5:
    #             self.nextSpectralFrame()
    #             self.nextSpectralFrame()
            
    #         if self.b > 5:
    #             self.nextSpectralFrame()
    #             self.nextSpectralFrame()
    #             self.nextSpectralFrame()
            
    #         if self.b < 0:
    #             self.nextSpectralFrame()

    #     # IF the direct timestamp is behind the spectral timestamp
    #     if self.a == -1:
            
    #         if self.b == 0:
    #             self.nextDirectFrame()
    #             self.nextDirectFrame()

    #         if 0 < self.b < 5:
    #             self.nextDirectFrame()
    #             self.nextDirectFrame()
            
    #         if self.b > 5:
    #             self.nextDirectFrame()
    #             self.nextDirectFrame()
    #             self.nextDirectFrame()
            
    #         if self.b < 0:
    #             self.nextDirectFrame()
        
    #     # If the direct timestamp is ahead by more than a tenth of a second
    #     if self.a > 1:

    #         print ("Time difference too great to apply time lock")

    #     # If the spectral timestamp is ahead by more than a tenth of a second
    #     if self.a < -1:
            
    #         print("Time difference too great to apply time lock")
    
    # # Last frames that are as close in time as possible
    # def lastTimeFrame(self):
    #     """
    #     Finds last frame set that is as close in time as possible.
    #     Evaluates tenths and hundredths of seconds
    #     portion of timestamp, finds the difference
    #     between the direct and spectral timestamps, then recognizes the time 
    #     gap and finds the appropriate method by which to get the next closest
    #     frames. 

    #     Does not work if the time gap is too great - user must have the 
    #     timestamp within 1 tenth of a second.
    #     """

    #     # If the tenths values are the same
    #     if self.a == 0:

    #         if 0 < self.b < 5:
    #             self.lastFrame()

    #         if 5 <= self.b < 9:
    #             self.lastSpectralFrame()
            
    #         if 0 >= self.b > -5:
    #             self.lastSpectralFrame()
    #             self.lastDirectFrame()
            
    #         if -5 >= self.b > -10:
    #             self.lastDirectFrame()
        
    #     # If the direct timestamp is ahead of the spectral timestamp
    #     if self.a ==1:

    #         if self.b == 0:
    #             self.lastSpectralFrame()
    #             self.lastSpectralFrame()

    #         if 0 < self.b < 5:
    #             self.lastSpectralFrame()
    #             self.lastSpectralFrame()
            
    #         if self.b > 5:
    #             self.lastSpectralFrame()
    #             self.lastSpectralFrame()
    #             self.lastSpectralFrame()
            
    #         if self.b < 0:
    #             self.lastSpectralFrame()

    #     # IF the direct timestamp is behind the spectral timestamp
    #     if self.a == -1:
            
    #         if self.b == 0:
    #             self.lastDirectFrame()
    #             self.lastDirectFrame()

    #         if 0 < self.b < 5:
    #             self.lastDirectFrame()
    #             self.lastDirectFrame()
            
    #         if self.b > 5:
    #             self.lastDirectFrame()
    #             self.lastDirectFrame()
    #             self.lastDirectFrame()
            
    #         if self.b < 0:
    #             self.lastDirectFrame()
        
    #     # If the direct timestamp is ahead by more than a tenth of a second
    #     if self.a > 1:

    #         print ("Time difference too great to apply time lock")

    #     # If the spectral timestamp is ahead by more than a tenth of a second
    #     if self.a < -1:
            
    #         print("Time difference too great to apply time lock")
           
    ################# PLOTTING FUNCTIONS #################

     # Projects affine marker onto spectrum
    def projectAffine(self):
        """
        Takes the mapped affine marker and projects it onto the spectrum. 
        Uses the scene handles to parameterize the spectrum as a line. 
        """       

        # Obtaion coordinates of ROI handles relative to scene
        self.handles = self.spectral_roi.getSceneHandlePositions()

        # Extract scale and translation  handle information
        scale_handle = list(self.handles[1])
        translate_handle = list(self.handles[2])

        # Extract scale and translation handle positions as QpointF objects
        scale_handle_position  = scale_handle[1]
        translate_handle_position = translate_handle[1]

        # Get equation of line using the two points
        m = (scale_handle_position.y() - translate_handle_position.y() ) / (scale_handle_position.x() - translate_handle_position.x() )

        b1 = scale_handle_position.y() - m*scale_handle_position.x()


        b2 = ((self.hu / m) + self.hv)

        self.x = (b2 - b1) / (m + (1 / m))

        y = m*self.x + b1

    # plot the measured spectrum
    def plotMeasuredSpec(self):
        """
        Runs background and region  of interest calculations.
        Converts spectrum from pixels to nanometers. 
        Plots measured spectrum in graph area
        """

        # Check background and set region to be plotted
        self.checkSpectralBackground()
        self.checkSpectralRegion()
        self.projectAffine()
        # self.x = 400
        
        # Set pen  
        if self.PlottedSpectrumNumber == 0:
            pen = pg.mkPen(width = 2)
        else:
            pen = pg.mkPen(pg.intColor(self.PlottedSpectrumNumber-1,6), width = 1)
            if self.PlottedSpectrumNumber == 6:
                self.PlottedSpectrumNumber = 0
        self.PlottedSpectrumNumber += 1

        # Shear spectral array
        if self.SpectralShear_rollbox.value() != 0.0:
            zoom = 10
            shear = int(self.SpectralShear_rollbox.value())

            spectral_array_zoomed = scipy.ndimage.zoom(self.spectral_array, (1,zoom))
            roi_w = np.shape(spectral_array_zoomed)[1]

            spectral_array_sheared = np.zeros((np.shape(spectral_array_zoomed)[0], np.shape(spectral_array_zoomed)[1]))

            for i in range(np.shape(spectral_array_zoomed)[1]):
                spectral_array_sheared[:,i] = np.roll(spectral_array_zoomed[:,i],int((shear-i*shear/roi_w)))

            spectral_array_sheared = np.roll(spectral_array_sheared, int(np.round(-shear/2)), axis=0)
            spectral_array_unzoomed = scipy.ndimage.zoom(spectral_array_sheared, (1,1/zoom))
        else:
            spectral_array_unzoomed = self.spectral_array

        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        # plt.subplots(figsize=(676*px, 100*px))
        # plt.imshow(spectral_array_unzoomed.T)
        # plt.show()

        # Set spectral profile
        spectral_profile = (np.mean(spectral_array_unzoomed, axis=1) - self.BiasLevel_rollbox.value()) * self.YScaler_rollbox.value()

        # spectral_profile = np.mean(self.spectral_array, axis=1)

        # Init array for the scaled profile
        global scaled_spectral_profile
        scaled_spectral_profile = np.zeros(len(spectral_profile))

        # Scaling parameters
        # s = 1/2.15 # px/nm
        s = 1/self.SpectralScale_rollbox.value() # px/nm
        nm0 = float(self.nm0_edit.text()) # nm 

        # Calculate wavelength values as they correspond to each pixel
        for i in range(len(scaled_spectral_profile)):
            nmt = (((i - self.x) / s) + nm0)
            scaled_spectral_profile = np.append(scaled_spectral_profile, nmt)

        

        # Take the part of the array with the desired values
        length = len(scaled_spectral_profile)
        middle_index = length//2
        
        # Reset array
        scaled_spectral_profile = scaled_spectral_profile[middle_index:]

        # plt.figure(figsize=(16,8))
        # plt.plot(scaled_spectral_profile,spectral_profile)
        # plt.show()

        scaled_spectral_profile_short = []
        spectral_profile_short = []
        for i in range(len(scaled_spectral_profile)):
            if scaled_spectral_profile[i] > 300 and scaled_spectral_profile[i] < 1100:
                scaled_spectral_profile_short.append(scaled_spectral_profile[i])
                spectral_profile_short.append(spectral_profile[i])

        # print(min(scaled_spectral_profile))

        # Interpolate responsivity curve to match scaled_spectral_profile
        xM = self.responsivityDefault[:,0]
        # xM = np.linspace(300,1100,800)
        # print(self.responsivityDefault[:,1])
        yM = self.responsivityDefault[:,1]/np.max(self.responsivityDefault[:,1])
        # yM = np.ones(800)
        fM = interpolate.interp1d(xM,yM)
        yMnew = fM(scaled_spectral_profile_short)

        # print(len(yMnew))
        # print(np.min(scaled_spectral_profile), np.max(scaled_spectral_profile), len(scaled_spectral_profile))

        self.plotMax = np.max(np.divide(spectral_profile_short,yMnew))

        # Set axis titles 
        self.Plot.setLabel('left', 'Intensity')
        self.Plot.setLabel('bottom', 'Wavelength (nm)')

        # Create the plot
        self.Plot.plot([518,518],[0,self.plotMax], pen=pg.mkPen(color=(0,0,0), style=QtCore.Qt.DotLine, width=2))
        self.Plot.plot([589,589],[0,self.plotMax], pen=pg.mkPen(color=(0,0,0), style=QtCore.Qt.DotLine, width=2))
        self.Plot.plot([777,777],[0,self.plotMax], pen=pg.mkPen(color=(0,0,0), style=QtCore.Qt.DotLine, width=2))

        self.Plot.plot(scaled_spectral_profile_short, spectral_profile_short, pen = pg.mkPen(color=(50,50,50), width = 2))
        self.Plot.plot(scaled_spectral_profile_short, np.divide(spectral_profile_short,yMnew), pen=pg.mkPen(color=(0,0,255), width=2)) # Uncomment for responsivity

        # line = pg.InfiniteLine(pos=589, angle=90, pen=pen)

        self.t1 = pg.TextItem('518.2 nm', anchor=(1,1), angle=90)
        self.t2 = pg.TextItem('589.1 nm', anchor=(1,1), angle=90)
        self.t3 = pg.TextItem('777.1 nm', anchor=(1,1), angle=90)
        self.Plot.addItem(self.t1)
        self.Plot.addItem(self.t2)
        self.Plot.addItem(self.t3)
        self.t1.setPos(518,self.plotMax)
        self.t2.setPos(589,self.plotMax)
        self.t3.setPos(777,self.plotMax)
        
        self.Plot.setXRange(np.min(scaled_spectral_profile_short),np.max(scaled_spectral_profile_short))
        self.Plot.setYRange(0,self.plotMax)
        self.Plot.setBackground('w')
        # self.Plot.setYRange(0,30000)
        self.CalibrateSpectrum_button.setEnabled(True)
        self.SetReference_button.setEnabled(True)

        

    def setReference(self):
        self.SetReference_button.setText('Remove Reference')

    # clear the spectrum
    def clearSpec(self):
        self.PlottedSpectrumNumber = 0
        self.Plot.clear()

    def savePlot(self):
        print('Saving plot...')

        plt.figure(figsize=(14,10))
        plt.show()

    # def updateFlatName(self):
    #     """ update flat structure when FlatName_linedit is selected and enter pressed """
    #     flat_path = self.FlatPath_linedit.text()
    #     flat_name = self.FlatName_linedit.text()

    #     if os.path.exists(os.path.join(flat_path, flat_name)):
    #         self.flat_structure = loadFlat(flat_path, flat_name)



    ##### Command Button Helper Functions  #####

    # def addExtinction(self, increment=0.1):
    #     """ add increment to the current exctinction coefficient value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.AtExtinct_linedit.text()
    #     if not current_value:
    #         current_value = 0.9

    #     # add, while handling variable types and floating point wandering
    #     added_value = float(current_value) + increment
    #     added_value = np.around(added_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.AtExtinct_linedit.setText(str(added_value))
    #     self.updateExtinctionValue()

    # def subExtinction(self, increment=0.1):
    #     """ subtract increment from the current exctinction coefficient value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.AtExtinct_linedit.text()
    #     if not current_value:
    #         current_value = 1.1

    #     # subtract, while handling variable types and floating point wandering
    #     subbed_value = float(current_value) - increment
    #     subbed_value = np.around(subbed_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.AtExtinct_linedit.setText(str(subbed_value))
    #     self.updateExtinctionValue

    # def updateExtinctionValue(self):

    #     extinctionValue = self.AtExtinct_linedit.text()
    #     if extinctionValue:
    #         extinctionValue = float(extinctionValue)
    #     else:
    #         return

        # update GuralSpectral object

    # Modified for spinner - MJM
    def updateExtinctionValue(self):
        self.extinctionValue = self.Extinction_rollbox.value()
        print('Updated extinction vale: %f' % self.extinctionValue)

    # def addRoll(self, increment=0.1):
    #     """ add increment to the current roll value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.Roll_linedit.text()
    #     if not current_value:
    #         current_value = 0.9

    #     # add, while handling variable types and floating point wandering
    #     added_value = float(current_value) + increment
    #     added_value = np.around(added_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.Roll_linedit.setText(str(added_value))
    #     self.updateRollValue()

    # def subRoll(self, increment=0.1):
    #     """ subtract increment from the current roll value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.Roll_linedit.text()
    #     if not current_value:
    #         current_value = 1.1

    #     # subtract, while handling variable types and floating point wandering
    #     subbed_value = float(current_value) - increment
    #     subbed_value = np.around(subbed_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.Roll_linedit.setText(str(subbed_value))
    #     self.updateRollValue()

    # def updateRollValue(self):

    #     rollValue = self.Roll_linedit.text()
    #     if rollValue:
    #         rollValue = float(rollValue)
    #     else:
    #         return

        # update GuralSpectral object

    # Modified for spinner - MJM
    def updateRollValue(self):
        rollValue = self.Roll_rollbox.value()


    # def addLmm(self, increment=0.1):
    #     """ add increment to the current L/mm value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.Lmm_linedit.text()
    #     if not current_value:
    #         current_value = 0.9

    #     # add, while handling variable types and floating point wandering
    #     added_value = float(current_value) + increment
    #     added_value = np.around(added_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.Lmm_linedit.setText(str(added_value))
    #     self.updateLmmValue()

    # def subLmm(self, increment=0.1):
    #     """ subtract increment from the currentL/mm value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.Lmm_linedit.text()
    #     if not current_value:
    #         current_value = 1.1

    #     # subtract, while handling variable types and floating point wandering
    #     subbed_value = float(current_value) - increment
    #     subbed_value = np.around(subbed_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.Lmm_linedit.setText(str(subbed_value))
    #     self.updateLmmValue()

    # def updateLmmValue(self):

    #     LmmValue = self.Lmm_linedit.text()
    #     if LmmValue:
    #         LmmValue = float(LmmValue)
    #     else:
    #         return

    #     # update GuralSpectral object

    # Modified for spinner - MJM
    def updateLmmValue(self):
        lmmValue = self.Lmm_rollbox.value()


    # def addHighTemp(self, increment=100):
    #     """ add increment to the current Temp high value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.HighTemp_linedit.text()
    #     if not current_value:
    #         current_value = 9900

    #     # add, while handling variable types and floating point wandering
    #     added_value = float(current_value) + increment
    #     added_value = np.around(added_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.HighTemp_linedit.setText(str(added_value))
    #     self.updateHighTempValue()

    # def subHighTemp(self, increment=100):
    #     """ subtract increment from the current Temp high value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.HighTemp_linedit.text()
    #     if not current_value:
    #         current_value = 10100

    #     # subtract, while handling variable types and floating point wandering
    #     subbed_value = float(current_value) - increment
    #     subbed_value = np.around(subbed_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.HighTemp_linedit.setText(str(subbed_value))
    #     self.updateHighTempValue()

    # def updateHighTempValue(self):

    #     highTempValue = self.HighTemp_linedit.text()
    #     if highTempValue:
    #         highTempValue = float(highTempValue)
    #     else:
    #         return

    #     # update GuralSpectral object

    # Modified for spinner - MJM
    def updateHighTempValue(self):
        highTempValue = self.HighTemp_rollbox.value()


    # def addLowTemp(self, increment=100):
    #     """ add increment to the current Temp low value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.LowTemp_linedit.text()
    #     if not current_value:
    #         current_value = 4400

    #     # add, while handling variable types and floating point wandering
    #     added_value = float(current_value) + increment
    #     added_value = np.around(added_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.LowTemp_linedit.setText(str(added_value))
    #     self.updateLowTempValue()

    # def subLowTemp(self, increment=100):
    #     """ subtract increment from the current Temp low value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.LowTemp_linedit.text()
    #     if not current_value:
    #         current_value = 4600

    #     # subtract, while handling variable types and floating point wandering
    #     subbed_value = float(current_value) - increment
    #     subbed_value = np.around(subbed_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.LowTemp_linedit.setText(str(subbed_value))
    #     self.updateLowTempValue()

    # def updateLowTempValue(self):

    #     lowTempValue = self.LowTemp_linedit.text()
    #     if lowTempValue:
    #         lowTempValue = float(lowTempValue)
    #     else:
    #         return

    #     # update GuralSpectral object

    # Modified for spinner - MJM
    def updateLowTempValue(self):
        lowTempValue = self.LowTemp_rollbox.value()


    # def addSigma(self, increment=0.1):
    #     """ add increment to the current sigma value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.Sigma_linedit.text()
    #     if not current_value:
    #         current_value = 0.9

    #     # add, while handling variable types and floating point wandering
    #     added_value = float(current_value) + increment
    #     added_value = np.around(added_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.Sigma_linedit.setText(str(added_value))
    #     self.updateSigmaValue()

    # def subSigma(self, increment=0.1):
    #     """ subtract increment from the current Sigma value """

    #     # get the current value, and if it is not set, assign it a value
    #     current_value = self.Sigma_linedit.text()
    #     if not current_value:
    #         current_value = 1.1

    #     # subtract, while handling variable types and floating point wandering
    #     subbed_value = float(current_value) - increment
    #     subbed_value = np.around(subbed_value, decimals=6)

    #     # set value and update in GuralSpectral object
    #     self.Sigma_linedit.setText(str(subbed_value))
    #     self.updateSigmaValue()

    # def updateSigmaValue(self):

    #     sigmaValue = self.Sigma_linedit.text()
    #     if sigmaValue:
    #         sigmaValue = float(sigmaValue)
    #     else:
    #         return

    #     # update GuralSpectral object

    def updateSigmaValue(self):
        # self.sigmaValue = self.Sigma_rollbox.value()
        self.spectral.changeBroadening(self.Sigma_rollbox.value())
        # self.refreshPlot()

    def updateHot2WarmRatio(self):
        print('Changing hot2warm ratio')
        self.spectral.changeHot2WarmRatio(self.Hot2WarmRatio_rollbox.value())
        # self.refreshPlot()
        # self.hot2Warm = self.Hot2WarmRatio_rollbox.value()
        # self.spectral.elemdata.hot2warm = self.Hot2WarmRatio_rollbox.value()


    def hotTempToggle(self):
        """ handle the toggling of the Hot button """
        if self.HotTempOn_button.isChecked():
            pass # do focus
        else:
            pass # do unfocus

    def warmTempToggle(self):
        """ handle the toggling of the Warm button """
        if self.WarmTempOn_button.isChecked():
            pass # do focus
        else:
            pass # do unfocus

    def ionsToggle(self):
        """ handle the toggling of the Ions button """
        if self.Ions_button.isChecked():
            pass # do focus
        else:
            pass # do unfocus

    def neutralToggle(self):
        """ handle the toggling of the Neutral button """
        if self.Neutral_button.isChecked():
            pass # do focus
        else:
            pass # do unfocus

    def responsivityToggle(self):
        """ handle the toggling of the Responsivity button """
        if self.Responsivity_check.isChecked():
            pass # do focus
        else:
            pass # do unfocus

    def extinctionToggle(self):
        """ handle the toggling of the Extinction button """
        if self.Extinction_check.isChecked():
            pass # do focus
        else:
            pass # do unfocus

###############################################################################################
################################ /// OPEN THE APPLICATION /// #################################
###############################################################################################

app = QtWidgets.QApplication(sys.argv) # create instance of QtWidgets.QApplication
window = Ui()                          # create instance of class
app.exec_()                            # start the application
