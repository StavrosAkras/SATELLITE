# testing_unequal_arrays_script.py:
# Determines the lenght of 2 arrays 
# (C) Stavros Akras

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import array
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
from astropy.io import fits
from astropy.io.fits import getdata
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from scipy.interpolate import spline
import matplotlib.patches as patches
import matplotlib.font_manager as font_manager
import random
import math


def testing_unequal_arrays(test,test_error,test_radial,pixscale,slitlength):
    
    for j in range(0,slitlength):
        for i in range(0,len(test_radial)):
            if test_radial[i]!=i*pixscale:
                test_radial[i:i] = [i*pixscale] 
                test[i:i] = [0] 
                test_error[i:i] = [0]            

    
    return test,test_error,test_radial





