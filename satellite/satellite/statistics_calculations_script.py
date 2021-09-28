# statistics_calculations_script.py:
# Computes and returns the 5%, 25%, 50%, 75%, and 95% percentiles
# as well as the mean value and standard deviation from a 2D maps
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




def statistic_numbers(data_array,sx,sy):
    
    data=[]
    k=0
    for j in range(0,sy):
        for i in range(0,sx):
            if data_array[i,j]!=0:
                data.append(0)
                data[k]=data_array[i,j]
                k=k+1
    
    a1=k
    a2=np.percentile(data, 5)
    a3=np.percentile(data, 25)
    a4=np.percentile(data, 50)
    a5=np.percentile(data, 70)
    a6=np.percentile(data, 95)
    a7=np.mean(data)
    a8=np.std(data)
    
    return a1,a2,a3,a4,a5,a6,a7,a8
