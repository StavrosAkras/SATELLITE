# physical_parameters_script.py:
# This routine is an intermedaite step in which the TeNe routine is used
# for each module to compute the Te and Ne parameters.
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
import seaborn as sns
import pyneb as pn

from satellite import diagnotic_diagrams_script as dds
from satellite import TeNe_angles_script as TeNeangles
from satellite import TeNe_specific_slits_script as TeNeslits
from satellite import TeNe_2D_script as TeNe2Ds

global file10

        
def analysiswithPyNeb(flux2D,flux2D_error,flux_angles,flux_angles_error,flux_spec_slit,flux_spec_slit_error,ang,line_names,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,hdr):

    
    ####################
    ### for all the angle in the rotation angle task
    ####################
    index=param_estimated.index('angular_analysis_task')
    if param_requered[index]=="yes":
        Te_PA, Ne_PA = TeNeangles.TeNe(flux_angles,flux_angles_error,ang,line_names,lines_available,param_estimated,param_requered)
    

    ####################
    ### for the 2D analysis
    ####################
    index=param_estimated.index('2D_anaysis_task')
    if param_requered[index]=="yes":
        Te_2D, Ne_2D = TeNe2Ds.TeNe(flux2D,flux2D_error,line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered,hdr)
    
    ####################
    ### for all the specific slit regions defined by the user 
    ####################
    index=param_estimated.index('specific_slits_analysis_task')
    if param_requered[index]=="yes":
        Te_slits, Ne_slits = TeNeslits.TeNe(flux_spec_slit,flux_spec_slit_error,line_names,lines_available,param_estimated,param_requered)
    
    

    return      

        
