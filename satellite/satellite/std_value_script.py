# std_value_script.py:
# Computes and returns the standard deviation of a 1D array (e.g. emission line ratio table).
# The "value" parameter defines the minimum value of the table than will be considered in the 
# estimation of the standard deviation, excluding the pixels wth a value of -100000. 
# The mean has been calculated in the mean_value_script.
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
import matplotlib.patches as patches
import matplotlib.font_manager as font_manager
import random
import math



def stdvalue(a1,mean,value):
        sum=0.0
        l=0
        for j in range(0,len(a1)):
                if a1[j] > value:        
                    sum=sum+(a1[j]-mean)*(a1[j]-mean)
                    l=l+1
                    
        return np.sqrt(sum/(l-1))
