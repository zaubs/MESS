##testing file

###################################################################################################
################################## IMPORTS AND SUPPORT FUNCTIONS ##################################
###################################################################################################

#################### STANDARD IMPORTS ####################

import numpy as np
import pyqtgraph as pg
import sys 
sys.path.append('../RMS/RMS/Routines')

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

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
from astropy.modeling import models
from astropy import units as u
from astropy.visualization import quantity_support

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
from RMS.Routines.BinImageCy import binImage as binImageCy

spectral_library = __import__("CAMO-Spectral_Library")

########TESTING OUT THE SHOW SPECTRUM#############



spec = np.genfromtxt('testspectrum.txt', delimiter = ' ', skip_header = 7)
resp = np.genfromtxt('DefaultResponsivity-EMCCD.csv', delimiter = ',', skip_header = 2)
#resp = np.genfromtxt('DefaultResponsivity-EMCCD2.csv', delimiter = ',', )
x = spec[:,0]
y = spec[:,1]
xM = resp[:,0]
yMsg = savgol_filter(resp[:,1],101,2)
yM = yMsg/np.max(yMsg)
fM = interpolate.interp1d(xM,yM)
yMnew = -1/fM(x)
y_corr = np.divide(y*1000,yMnew- np.min(yMnew))

plt.figure()
plt.plot(x,y)
plt.plot(x,y_corr,color = 'orange')
#plt.plot(yMnew)
plt.show()

#######INTENSITY PROFILE CHECK FUNCTION############
"""
import sys                                                                                          
from PySide2.QtWidgets import QApplication, QPushButton                                                                                                                                                                                                                                                                  

#hight, frame number, speed, list of values/dictionary for frames, uncompresss it to a /tmp directory 
#Giving path from the stream and everything from the weblog 
#20240905_025347A, get the name from the directory (after the _)

#use unix time to convert to standard time and convert to day time 
#Use values to the RIGHT of the variable 
#Beg and End has the begin and end value 

#tag is the station number 
# t - time 
# th - zenith angle (90 - elevation angle ) 
# phi - azmith angle (no clue which direction(either north or east) )
# tau - he forgo 
# cx, cy - Coordinates on the mirror system(not important)
# lsp - log on pixel values, 
# mag - magnitude, transformation*lsp + offset value 
# L,R - he forgo 
# vel - velocity 
# ht - height 
# lat - latitude 
# lon- longitude 

# vel_g is the final orbital solution(velocity)


import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# load Image as Grayscale
i = Image.open("afterFlat.png").convert("L")
 convert to numpy array
n = np.array(i)

# average columns and rows
# left to right
cols = n[12]
# bottom to top
rows = n.mean(axis=1)

# plot histograms
f, ax = plt.subplots(2, 1)
ax[0].plot(cols)
ax[1].plot(rows)

#histogram,bin_edges = np.histogram(i,bins = 256,range = [0,1])
plt.show()
# configure and draw the histogram figure
fig, ax = plt.subplots()
ax.set_title("Grayscale Histogram")
ax.set_xlabel("grayscale value")
ax.set_ylabel("pixel count")
ax.set_xlim([0.0, 1.0])  # <- named arguments do not work here

ax.plot(bin_edges[0:-1], histogram)  # <- or here

plt.show()
#from skimage.measure import profile_line

import matplotlib.pyplot as plt
import numpy as np
import skimage.measure 
from skimage import data 
import skimage 

#raw = imageio.v2.imread("RawFrame.png")
#tform = skimage.transform.SimilarityTransform(scale = 1.4,rotation = 0.968657735,translation=(raw.shape[0]/2 +20 ,- raw.shape[0]/2))
#image = skimage.transform.warp(raw,tform)

image = imageio.v2.imread("afterFlat_median.png")
start = (300,0) #Start of the profile line row=100, col=0
end = (300,image.shape[1] - 1) #End of the profile line row=100, col=last

profile = skimage.measure.profile_line(image, start, end)

fig, ax = plt.subplots(1, 2)
ax[0].set_title('The image fater np.median flat ')
ax[0].imshow(image)
ax[0].plot([start[1], end[1]], [start[0], end[0]], 'r')
ax[1].set_title('Intensity profile of the line')
ax[1].plot(profile/255)
plt.show()

"""