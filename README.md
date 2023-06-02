# MESS
 The **M**eteor **E**lemental **S**pectra **S**oftware is a package based on the CAMO-S code for reducing meteor spectra from a static camera.
 
 # Installation
 To successfully run the software, Python and QT are required along with the following Python modules...
 - Cython
 - Numpy
 - Matplotlib
 - pyqtgraph
 - scipy
 - imageio
 - sklearn

 Installation in a virtualenv or a Conda environment is recommended.
 
 SpectralTest.so.0 in spectral_library should be soft-linked to SpectralTest.so
 ln -s SpectralTest.so.0 SpectralTest.so
 
 **Input data types:** VID files, ev* txt files, PNG images<br>
 **Output data types:** CSV and PNG images
 
# Using MESS
