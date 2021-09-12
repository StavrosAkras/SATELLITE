# define_Te_Ne_for_null_spaxels_script.py: 
# Defines the Te and Ne classes of arrays with zero values.
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
import numpy as np
from scipy.stats import norm
import seaborn as sns
import sys


def nullTeNe(param_mod_name,param_model_values):
    
    index_step_ang=param_mod_name.index("position_angle_step")
    index_star_ang=param_mod_name.index("minimum_position_angle")
    index_end_ang=param_mod_name.index("maximum_position_angle")
    start_angle=param_model_values[index_star_ang]
    end_angle=param_model_values[index_end_ang]
    step_angle=param_model_values[index_step_ang]
    
    numbersize=int((end_angle-start_angle)/step_angle)+1
    if numbersize>0:
        class Te:
            NIISII=np.zeros(numbersize)
            OISII=np.zeros(numbersize)
            OIISII=np.zeros(numbersize)
            OIIISII=np.zeros(numbersize)
            SIIISII=np.zeros(numbersize)
            OIIOII=np.zeros(numbersize)
            NIIOII=np.zeros(numbersize)
            OIOII=np.zeros(numbersize)
            OIIIClIII=np.zeros(numbersize)
            SIIIClIII=np.zeros(numbersize)
            OIIIArVI=np.zeros(numbersize)
            SIIIArVI=np.zeros(numbersize)
            NIIClIII=np.zeros(numbersize)

        class Ne:
            NIISII=np.zeros(numbersize)
            OISII=np.zeros(numbersize)
            OIISII=np.zeros(numbersize)
            OIIISII=np.zeros(numbersize)
            SIIISII=np.zeros(numbersize)
            OIIOII=np.zeros(numbersize)
            NIIOII=np.zeros(numbersize)
            OIOII=np.zeros(numbersize)
            OIIIClIII=np.zeros(numbersize)
            SIIIClIII=np.zeros(numbersize)
            OIIIArVI=np.zeros(numbersize)
            SIIIArVI=np.zeros(numbersize)
            NIIClIII=np.zeros(numbersize)
        
        class Te_error:
            NIISII=np.zeros(numbersize)
            OISII=np.zeros(numbersize)
            OIISII=np.zeros(numbersize)
            OIIISII=np.zeros(numbersize)
            SIIISII=np.zeros(numbersize)
            OIIOII=np.zeros(numbersize)
            NIIOII=np.zeros(numbersize)
            OIOII=np.zeros(numbersize)
            OIIIClIII=np.zeros(numbersize)
            SIIIClIII=np.zeros(numbersize)
            OIIIArVI=np.zeros(numbersize)
            SIIIArVI=np.zeros(numbersize)
            NIIClIII=np.zeros(numbersize)
        
        class Ne_error:
            NIISII=np.zeros(numbersize)
            OISII=np.zeros(numbersize)
            OIISII=np.zeros(numbersize)
            OIIISII=np.zeros(numbersize)
            SIIISII=np.zeros(numbersize)
            OIIOII=np.zeros(numbersize)
            NIIOII=np.zeros(numbersize)
            OIOII=np.zeros(numbersize)
            OIIIClIII=np.zeros(numbersize)
            SIIIClIII=np.zeros(numbersize)
            OIIIArVI=np.zeros(numbersize)
            SIIIArVI=np.zeros(numbersize)
            NIIClIII=np.zeros(numbersize)
    
    return Te, Ne


def nullTeNe_specific_slits(param_mod_name,param_model_values):
    
    number_specific_slits=10
    class Te:
        NIISII=np.zeros(number_specific_slits)
        OISII=np.zeros(number_specific_slits)
        OIISII=np.zeros(number_specific_slits)
        OIIISII=np.zeros(number_specific_slits)
        SIIISII=np.zeros(number_specific_slits)
        OIIOII=np.zeros(number_specific_slits)
        NIIOII=np.zeros(number_specific_slits)
        OIOII=np.zeros(number_specific_slits)
        OIIIClIII=np.zeros(number_specific_slits)
        SIIIClIII=np.zeros(number_specific_slits)
        OIIIArVI=np.zeros(number_specific_slits)
        SIIIArVI=np.zeros(number_specific_slits)
        NIIClIII=np.zeros(number_specific_slits)

    class Ne:
        NIISII=np.zeros(number_specific_slits)
        OISII=np.zeros(number_specific_slits)
        OIISII=np.zeros(number_specific_slits)
        OIIISII=np.zeros(number_specific_slits)
        SIIISII=np.zeros(number_specific_slits)
        OIIOII=np.zeros(number_specific_slits)
        NIIOII=np.zeros(number_specific_slits)
        OIOII=np.zeros(number_specific_slits)
        OIIIClIII=np.zeros(number_specific_slits)
        SIIIClIII=np.zeros(number_specific_slits)
        OIIIArVI=np.zeros(number_specific_slits)
        SIIIArVI=np.zeros(number_specific_slits)
        NIIClIII=np.zeros(number_specific_slits)
        
    class Te_error:
        NIISII=np.zeros(number_specific_slits)
        OISII=np.zeros(number_specific_slits)
        OIISII=np.zeros(number_specific_slits)
        OIIISII=np.zeros(number_specific_slits)
        SIIISII=np.zeros(number_specific_slits)
        OIIOII=np.zeros(number_specific_slits)
        NIIOII=np.zeros(number_specific_slits)
        OIOII=np.zeros(number_specific_slits)
        OIIIClIII=np.zeros(number_specific_slits)
        SIIIClIII=np.zeros(number_specific_slits)
        OIIIArVI=np.zeros(number_specific_slits)
        SIIIArVI=np.zeros(number_specific_slits)
        NIIClIII=np.zeros(number_specific_slits)
    
    class Ne_error:
        NIISII=np.zeros(number_specific_slits)
        OISII=np.zeros(number_specific_slits)
        OIISII=np.zeros(number_specific_slits)
        OIIISII=np.zeros(number_specific_slits)
        SIIISII=np.zeros(number_specific_slits)
        OIIOII=np.zeros(number_specific_slits)
        NIIOII=np.zeros(number_specific_slits)
        OIOII=np.zeros(number_specific_slits)
        OIIIClIII=np.zeros(number_specific_slits)
        SIIIClIII=np.zeros(number_specific_slits)
        OIIIArVI=np.zeros(number_specific_slits)
        SIIIArVI=np.zeros(number_specific_slits)
        NIIClIII=np.zeros(number_specific_slits)
    
    return Te, Ne
