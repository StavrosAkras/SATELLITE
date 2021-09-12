# my_analysis2D_script.py:
# The user can develop his/her own script for 2D analysis module ( future release)
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


from satellite import diagnotic_diagrams_script as dds

        
def analysis2D(flux2D,flux_angles_norm,ang,param_mod_name,param_model_values):
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]
    
### define a class for the normalized fluxes to Hb=100  (intensities)
    class flux_norm_2D:
        Ha_6563=np.zeros((sizex, sizey))
        Hb_4861=np.zeros((sizex, sizey))
        Hg_4340=np.zeros((sizex, sizey))
        Hd_4101=np.zeros((sizex, sizey))
        HeIa_5876=np.zeros((sizex, sizey))
        HeIb_6678=np.zeros((sizex, sizey))
        HeIIa_4686=np.zeros((sizex, sizey))
        HeIIb_5412=np.zeros((sizex, sizey))
        NIIa_5755=np.zeros((sizex, sizey))
        NIIb_6548=np.zeros((sizex, sizey))
        NIIc_6584=np.zeros((sizex, sizey))
        NI_5199=np.zeros((sizex, sizey))
        OIIIa_4363=np.zeros((sizex, sizey))
        OIIIb_4959=np.zeros((sizex, sizey))
        OIIIc_5007=np.zeros((sizex, sizey))
        OIIa_3727=np.zeros((sizex, sizey))
        OIIb_3729=np.zeros((sizex, sizey))
        OIIc_7320=np.zeros((sizex, sizey))
        OIId_7330=np.zeros((sizex, sizey))
        OIa_5577=np.zeros((sizex, sizey))
        OIb_6300=np.zeros((sizex, sizey))
        OIc_6363=np.zeros((sizex, sizey))
        SIIa_6716=np.zeros((sizex, sizey))
        SIIb_6731=np.zeros((sizex, sizey))
        SIIIa_6312=np.zeros((sizex, sizey))
        SIIIb_9069=np.zeros((sizex, sizey))
        ClIIIa_5517=np.zeros((sizex, sizey))
        ClIIIb_5538=np.zeros((sizex, sizey))
        ArIII_7136=np.zeros((sizex, sizey))
        ArIVa_4712=np.zeros((sizex, sizey))
        ArIVb_4740=np.zeros((sizex, sizey))
        CI_8727=np.zeros((sizex, sizey))
        CII_6461=np.zeros((sizex, sizey))
        NeIIIa_3868=np.zeros((sizex, sizey))
        NeIIIb_3967=np.zeros((sizex, sizey))


### define a class for the error of the normalized fluxes to Hb=100  (intensities)
    class flux_norm_2D_error:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]




    return      

        
