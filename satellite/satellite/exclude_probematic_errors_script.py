# exclude_probematic_errors_script.py:
# Finds and excludes the low signal-to-noise pixels. 
# The minimum S/N is given by the parameter "value"
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


def exclude_problematic_error_values(flux2D,flux2D_error,param_mod_name,param_model_values):
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]
    
    
    class flux2D_new:
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

    class flux2D_new_error:
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
    
######################################################################    
###this value defined the limit (SNR) of the acceptanced measurements 
######################################################################

    value=2.0 
    for j in range(0,sizey):
        for i in range(0,sizex):
            if flux2D_error.Ha_6563[i,j]/flux2D.Ha_6563[i,j]<value and flux2D.Ha_6563[i,j]!=0:
                flux2D_new.Ha_6563[i,j]=flux2D.Ha_6563[i,j]
                flux2D_new_error.Ha_6563[i,j]=flux2D_error.Ha_6563[i,j]
            else:
                flux2D_new.Ha_6563[i,j]=0.0
                flux2D_new_error.Ha_6563[i,j]=0.0
            if flux2D_error.Hb_4861[i,j]/flux2D.Hb_4861[i,j]<value and flux2D.Hb_4861[i,j]!=0:
                flux2D_new.Hb_4861[i,j]=flux2D.Hb_4861[i,j]
                flux2D_new_error.Hb_4861[i,j]=flux2D_error.Hb_4861[i,j]
            else:
                flux2D_new.Hb_4861[i,j]=0.0
                flux2D_new_error.Hb_4861[i,j]=0.0

            if flux2D_error.Hg_4340[i,j]/flux2D.Hg_4340[i,j]<value and flux2D.Hg_4340[i,j]!=0:
                flux2D_new.Hg_4340[i,j]=flux2D.Hg_4340[i,j]
                flux2D_new_error.Hg_4340[i,j]=flux2D_error.Hg_4340[i,j]
            else:
                flux2D_new.Hg_4340[i,j]=0.0
                flux2D_new_error.Hg_4340[i,j]=0.0    
                
            if flux2D_error.Hd_4101[i,j]/flux2D.Hd_4101[i,j]<value and flux2D.Hd_4101[i,j]!=0:
                flux2D_new.Hd_4101[i,j]=flux2D.Hd_4101[i,j]
                flux2D_new_error.Hd_4101[i,j]=flux2D_error.Hd_4101[i,j]
            else:
                flux2D_new.Hd_4101[i,j]=0.0
                flux2D_new_error.Hd_4101[i,j]=0.0    
                
            if flux2D_error.HeIa_5876[i,j]/flux2D.HeIa_5876[i,j]<value and flux2D.HeIa_5876[i,j]!=0:
                flux2D_new.HeIa_5876[i,j]=flux2D.HeIa_5876[i,j]
                flux2D_new_error.HeIa_5876[i,j]=flux2D_error.HeIa_5876[i,j]
            else:
                flux2D_new.HeIa_5876[i,j]=0.0
                flux2D_new_error.HeIa_5876[i,j]=0.0  
                
            if flux2D_error.HeIb_6678[i,j]/flux2D.HeIb_6678[i,j]<value and flux2D.HeIb_6678[i,j]!=0:
                flux2D_new.HeIb_6678[i,j]=flux2D.HeIb_6678[i,j]
                flux2D_new_error.HeIb_6678[i,j]=flux2D_error.HeIb_6678[i,j]
            else:
                flux2D_new.HeIb_6678[i,j]=0.0
                flux2D_new_error.HeIb_6678[i,j]=0.0 
                
            if flux2D_error.HeIIa_4686[i,j]/flux2D.HeIIa_4686[i,j]<value and flux2D.HeIIa_4686[i,j]!=0:
                flux2D_new.HeIIa_4686[i,j]=flux2D.HeIIa_4686[i,j]
                flux2D_new_error.HeIIa_4686[i,j]=flux2D_error.HeIIa_4686[i,j]
            else:
                flux2D_new.HeIIa_4686[i,j]=0.0
                flux2D_new_error.HeIIa_4686[i,j]=0.0
           
            if flux2D_error.HeIIb_5412[i,j]/flux2D.HeIIb_5412[i,j]<value and flux2D.HeIIb_5412[i,j]!=0:
                flux2D_new.HeIIb_5412[i,j]=flux2D.HeIIb_5412[i,j]
                flux2D_new_error.HeIIb_5412[i,j]=flux2D_error.HeIIb_5412[i,j]
            else:
                flux2D_new.HeIIb_5412[i,j]=0.0
                flux2D_new_error.HeIIb_5412[i,j]=0.0 
            
                
            if flux2D_error.NIIa_5755[i,j]/flux2D.NIIa_5755[i,j]<value and flux2D.NIIa_5755[i,j]!=0:
                flux2D_new.NIIa_5755[i,j]=flux2D.NIIa_5755[i,j]
                flux2D_new_error.NIIa_5755[i,j]=flux2D_error.NIIa_5755[i,j]
            else:
                flux2D_new.NIIa_5755[i,j]=0.0
                flux2D_new_error.NIIa_5755[i,j]=0.0      
            
            if flux2D_error.NIIb_6548[i,j]/flux2D.NIIb_6548[i,j]<value and flux2D.NIIb_6548[i,j]!=0:
                flux2D_new.NIIb_6548[i,j]=flux2D.NIIb_6548[i,j]
                flux2D_new_error.NIIb_6548[i,j]=flux2D_error.NIIb_6548[i,j]
            else:
                flux2D_new.NIIb_6548[i,j]=0.0
                flux2D_new_error.NIIb_6548[i,j]=0.0  
                
            if flux2D_error.NIIc_6584[i,j]/flux2D.NIIc_6584[i,j]<value and flux2D.NIIb_6548[i,j]!=0:
                flux2D_new.NIIc_6584[i,j]=flux2D.NIIc_6584[i,j]
                flux2D_new_error.NIIc_6584[i,j]=flux2D_error.NIIc_6584[i,j]
            else:
                flux2D_new.NIIc_6584[i,j]=0.0
                flux2D_new_error.NIIc_6584[i,j]=0.0  
                
            if flux2D_error.NI_5199[i,j]/flux2D.NI_5199[i,j]<value and flux2D.NIIc_6584[i,j]!=0:
                flux2D_new.NI_5199[i,j]=flux2D.NI_5199[i,j]
                flux2D_new_error.NI_5199[i,j]=flux2D_error.NI_5199[i,j]
            else:
                flux2D_new.NI_5199[i,j]=0.0
                flux2D_new_error.NI_5199[i,j]=0.0      
                
            if flux2D_error.OIIIa_4363[i,j]/flux2D.OIIIa_4363[i,j]<value and flux2D.OIIIa_4363[i,j]!=0:
                flux2D_new.OIIIa_4363[i,j]=flux2D.OIIIa_4363[i,j]
                flux2D_new_error.OIIIa_4363[i,j]=flux2D_error.OIIIa_4363[i,j]
            else:
                flux2D_new.OIIIa_4363[i,j]=0.0
                flux2D_new_error.OIIIa_4363[i,j]=0.0       
                
            if flux2D_error.OIIIb_4959[i,j]/flux2D.OIIIb_4959[i,j]<value and flux2D.OIIIb_4959[i,j]!=0:
                flux2D_new.OIIIb_4959[i,j]=flux2D.OIIIb_4959[i,j]
                flux2D_new_error.OIIIb_4959[i,j]=flux2D_error.OIIIb_4959[i,j]
            else:
                flux2D_new.OIIIb_4959[i,j]=0.0
                flux2D_new_error.OIIIb_4959[i,j]=0.0    
                
            if flux2D_error.OIIIc_5007[i,j]/flux2D.OIIIc_5007[i,j]<value and flux2D.OIIIc_5007[i,j]!=0:
                flux2D_new.OIIIc_5007[i,j]=flux2D.OIIIc_5007[i,j]
                flux2D_new_error.OIIIc_5007[i,j]=flux2D_error.OIIIc_5007[i,j]
            else:
                flux2D_new.OIIIa_4363[i,j]=0.0
                flux2D_new_error.OIIIa_4363[i,j]=0.0    
                
            if flux2D_error.OIIa_3727[i,j]/flux2D.OIIa_3727[i,j]<value and flux2D.OIIa_3727[i,j]!=0:
                flux2D_new.OIIa_3727[i,j]=flux2D.OIIa_3727[i,j]
                flux2D_new_error.OIIa_3727[i,j]=flux2D_error.OIIa_3727[i,j]
            else:
                flux2D_new.OIIa_3727[i,j]=0.0
                flux2D_new_error.OIIa_3727[i,j]=0.0    
                
            if flux2D_error.OIIb_3729[i,j]/flux2D.OIIb_3729[i,j]<value and flux2D.OIIb_3729[i,j]!=0:
                flux2D_new.OIIb_3729[i,j]=flux2D.OIIb_3729[i,j]
                flux2D_new_error.OIIb_3729[i,j]=flux2D_error.OIIb_3729[i,j]
            else:
                flux2D_new.OIIb_3729[i,j]=0.0
                flux2D_new_error.OIIb_3729[i,j]=0.0    
                
            if flux2D_error.OIIc_7320[i,j]/flux2D.OIIc_7320[i,j]<value and flux2D.OIIc_7320[i,j]!=0:
                flux2D_new.OIIc_7320[i,j]=flux2D.OIIc_7320[i,j]
                flux2D_new_error.OIIc_7320[i,j]=flux2D_error.OIIc_7320[i,j]
            else:
                flux2D_new.OIIc_7320[i,j]=0.0
                flux2D_new_error.OIIc_7320[i,j]=0.0        
                
            if flux2D_error.OIId_7330[i,j]/flux2D.OIId_7330[i,j]<value and flux2D.OIId_7330[i,j]!=0:
                flux2D_new.OIId_7330[i,j]=flux2D.OIId_7330[i,j]
                flux2D_new_error.OIId_7330[i,j]=flux2D_error.OIId_7330[i,j]
            else:
                flux2D_new.OIId_7330[i,j]=0.0
                flux2D_new_error.OIId_7330[i,j]=0.0    
                
            if flux2D_error.OIa_5577[i,j]/flux2D.OIa_5577[i,j]<value and flux2D.OIa_5577[i,j]!=0:
                flux2D_new.OIa_5577[i,j]=flux2D.OIa_5577[i,j]
                flux2D_new_error.OIa_5577[i,j]=flux2D_error.OIa_5577[i,j]
            else:
                flux2D_new.OIa_5577[i,j]=0.0
                flux2D_new_error.OIa_5577[i,j]=0.0    
                
            if flux2D_error.OIb_6300[i,j]/flux2D.OIb_6300[i,j]<value and flux2D.OIb_6300[i,j]!=0:
                flux2D_new.OIb_6300[i,j]=flux2D.OIb_6300[i,j]
                flux2D_new_error.OIb_6300[i,j]=flux2D_error.OIb_6300[i,j]
            else:
                flux2D_new.OIb_6300[i,j]=0.0
                flux2D_new_error.OIb_6300[i,j]=0.0    
                
            if flux2D_error.OIc_6363[i,j]/flux2D.OIc_6363[i,j]<value and flux2D.OIc_6363[i,j]!=0:
                flux2D_new.OIc_6363[i,j]=flux2D.OIc_6363[i,j]
                flux2D_new_error.OIc_6363[i,j]=flux2D_error.OIc_6363[i,j]
            else:
                flux2D_new.OIc_6363[i,j]=0.0
                flux2D_new_error.OIc_6363[i,j]=0.0        
             
            if flux2D_error.SIIa_6716[i,j]/flux2D.SIIa_6716[i,j]<value and flux2D.SIIa_6716[i,j]!=0:
                flux2D_new.SIIa_6716[i,j]=flux2D.SIIa_6716[i,j]
                flux2D_new_error.SIIa_6716[i,j]=flux2D_error.SIIa_6716[i,j]
            else:
                flux2D_new.SIIa_6716[i,j]=0.0
                flux2D_new_error.SIIa_6716[i,j]=0.0 
                
            if flux2D_error.SIIb_6731[i,j]/flux2D.SIIb_6731[i,j]<value and flux2D.SIIb_6731[i,j]!=0:
                flux2D_new.SIIb_6731[i,j]=flux2D.SIIb_6731[i,j]
                flux2D_new_error.SIIb_6731[i,j]=flux2D_error.SIIb_6731[i,j]
            else:
                flux2D_new.SIIb_6731[i,j]=0.0
                flux2D_new_error.SIIb_6731[i,j]=0.0 
                
            if flux2D_error.SIIIa_6312[i,j]/flux2D.SIIIa_6312[i,j]<value and flux2D.SIIIa_6312[i,j]!=0:
                flux2D_new.SIIIa_6312[i,j]=flux2D.SIIIa_6312[i,j]
                flux2D_new_error.SIIIa_6312[i,j]=flux2D_error.SIIIa_6312[i,j]
            else:
                flux2D_new.SIIIa_6312[i,j]=0.0
                flux2D_new_error.SIIIa_6312[i,j]=0.0 
                
            if flux2D_error.SIIIb_9069[i,j]/flux2D.SIIIb_9069[i,j]<value and flux2D.SIIIb_9069[i,j]!=0:
                flux2D_new.SIIIb_9069[i,j]=flux2D.SIIIb_9069[i,j]
                flux2D_new_error.SIIIb_9069[i,j]=flux2D_error.SIIIb_9069[i,j]
            else:
                flux2D_new.SIIIb_9069[i,j]=0.0
                flux2D_new_error.SIIIb_9069[i,j]=0.0     
                
            if flux2D_error.ClIIIa_5517[i,j]/flux2D.ClIIIa_5517[i,j]<value and flux2D.ClIIIa_5517[i,j]!=0:
                flux2D_new.ClIIIa_5517[i,j]=flux2D.ClIIIa_5517[i,j]
                flux2D_new_error.ClIIIa_5517[i,j]=flux2D_error.ClIIIa_5517[i,j]
            else:
                flux2D_new.ClIIIa_5517[i,j]=0.0
                flux2D_new_error.ClIIIa_5517[i,j]=0.0
                
            if flux2D_error.ClIIIb_5538[i,j]/flux2D.ClIIIb_5538[i,j]<value and flux2D.ClIIIb_5538[i,j]!=0:
                flux2D_new.ClIIIb_5538[i,j]=flux2D.ClIIIb_5538[i,j]
                flux2D_new_error.ClIIIb_5538[i,j]=flux2D_error.ClIIIb_5538[i,j]
            else:
                flux2D_new.ClIIIb_5538[i,j]=0.0
                flux2D_new_error.ClIIIb_5538[i,j]=0.0
                
            if flux2D_error.NeIIIa_3868[i,j]/flux2D.NeIIIa_3868[i,j]<value and flux2D.NeIIIa_3868[i,j]!=0:
                flux2D_new.NeIIIa_3868[i,j]=flux2D.NeIIIa_3868[i,j]
                flux2D_new_error.NeIIIa_3868[i,j]=flux2D_error.NeIIIa_3868[i,j]
            else:
                flux2D_new.NeIIIa_3868[i,j]=0.0
                flux2D_new_error.NeIIIa_3868[i,j]=0.0
                
            if flux2D_error.NeIIIb_3967[i,j]/flux2D.NeIIIb_3967[i,j]<value and flux2D.NeIIIb_3967[i,j]!=0:
                flux2D_new.NeIIIb_3967[i,j]=flux2D.NeIIIb_3967[i,j]
                flux2D_new_error.NeIIIb_3967[i,j]=flux2D_error.NeIIIb_3967[i,j]
            else:
                flux2D_new.NeIIIb_3967[i,j]=0.0
                flux2D_new_error.NeIIIb_3967[i,j]=0.0    
             
            if flux2D_error.ArIII_7136[i,j]/flux2D.ArIII_7136[i,j]<value and flux2D.ArIII_7136[i,j]!=0:
                flux2D_new.ArIII_7136[i,j]=flux2D.ArIII_7136[i,j]
                flux2D_new_error.ArIII_7136[i,j]=flux2D_error.ArIII_7136[i,j]
            else:
                flux2D_new.ArIII_7136[i,j]=0.0
                flux2D_new_error.ArIII_7136[i,j]=0.0 
            
            if flux2D_error.ArIVa_4712[i,j]/flux2D.ArIVa_4712[i,j]<value and flux2D.ArIVa_4712[i,j]!=0:
                flux2D_new.ArIVa_4712[i,j]=flux2D.ArIVa_4712[i,j]
                flux2D_new_error.ArIVa_4712[i,j]=flux2D_error.ArIVa_4712[i,j]
            else:
                flux2D_new.ArIVa_4712[i,j]=0.0
                flux2D_new_error.ArIVa_4712[i,j]=0.0 
                
            if flux2D_error.ArIVb_4740[i,j]/flux2D.ArIVb_4740[i,j]<value and flux2D.ArIVb_4740[i,j]!=0:
                flux2D_new.ArIVb_4740[i,j]=flux2D.ArIVb_4740[i,j]
                flux2D_new_error.ArIVb_4740[i,j]=flux2D_error.ArIVb_4740[i,j]
            else:
                flux2D_new.ArIVb_4740[i,j]=0.0
                flux2D_new_error.ArIVb_4740[i,j]=0.0     
            
            if flux2D_error.CI_8727[i,j]/flux2D.CI_8727[i,j]<value and flux2D.CI_8727[i,j]!=0:
                flux2D_new.CI_8727[i,j]=flux2D.CI_8727[i,j]
                flux2D_new_error.CI_8727[i,j]=flux2D_error.CI_8727[i,j]
            else:
                flux2D_new.CI_8727[i,j]=0.0
                flux2D_new_error.CI_8727[i,j]=0.0 
            
            if flux2D_error.CII_6461[i,j]/flux2D.CII_6461[i,j]<value and flux2D.CII_6461[i,j]!=0:
                flux2D_new.CII_6461[i,j]=flux2D.CII_6461[i,j]
                flux2D_new_error.CII_6461[i,j]=flux2D_error.CII_6461[i,j]
            else:
                flux2D_new.CII_6461[i,j]=0.0
                flux2D_new_error.CII_6461[i,j]=0.0
            
                
    return flux2D_new,flux2D_new_error





